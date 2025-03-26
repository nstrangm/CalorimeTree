import yaml
import subprocess
import os
import re
import pathlib
import sys
import argparse
from collections import defaultdict
from difflib import SequenceMatcher
import time
import threading
import importlib
from multiprocessing import Process

try:
    importlib.import_module('rich')
except ImportError:
    activate_script = os.path.join("CaloEnv", "bin", "activate")
    if not os.path.exists(activate_script):
        print(f"Virtual environment activation script not found at {activate_script}")
        sys.exit(1)

    # Use the appropriate shell to activate the virtual environment
    # The current process is replaced by the new shell in a subshell.
    subprocess.run(["bash", "-c", f"source {activate_script} && exec python {' '.join(sys.argv)}"], shell=False)
    sys.exit()  # Exit after running the subprocess to avoid continuing with the current script execution

from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn, TimeElapsedColumn, SpinnerColumn
from rich.console import Console
import logging
from rich.logging import RichHandler
from rich.style import Style
from itertools import cycle
from rich.traceback import install
install()
from rich import print
from rich.panel import Panel

# Define a list of colors to cycle through
COLORS = ["red", "green", "yellow", "blue", "magenta", "cyan", "white"]
settings_colors = []

console = Console()
FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
log = logging.getLogger("rich")

def source_environment(script_path):
    """Source the script and return the updated environment variables."""
    command = f"source {script_path} && env"
    proc = subprocess.Popen(command, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    
    if proc.returncode != 0:
        raise RuntimeError(f"Failed to source environment: {err.decode().strip()}")
    
    env = {}
    for line in out.decode().splitlines():
        key, _, value = line.partition("=")
        env[key] = value
    return env

def load_root_environment():
    """Ensure the ROOT environment is loaded."""
    try:
        # Check if ROOT is already loaded
        subprocess.run(['root-config', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log.info(f"ROOT environment is already loaded.")
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Load ROOT environment
        root_setup_script = "/software/nstrangmann/root/bin/thisroot.sh"  # Update this path
        if os.path.exists(root_setup_script):
            new_env = source_environment(root_setup_script)
            os.environ.update(new_env)
            log.info(f"ROOT environment loaded successfully.")
        else:
            raise FileNotFoundError("ROOT environment setup script not found.")

# Do Event mixing (takes MB and triggered trees, mixes them using  Gale-Shapley stable matching algorithm as root files)
def do_event_mixing():
    print("nothing")

# Read Cuts.yaml file to create folders in IsoGammaJet: In all folders (datasets) create folder for each yaml config
def create_folders():
    print("nothing")

# Function to read YAML files
def read_yaml(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

# def find_files(directory, pattern, output_file):
#     with open(output_file, 'w') as f:
#         for root, dirs, files in os.walk(directory):
#             print(f"Searching in: {root}")  # Debug statement
#             for file in files:
#                 if file.endswith(pattern):
#                     file_path = os.path.join(root, file)
#                     print(f"Found: {file_path}")  # Debug statement
#                     f.write(file_path + '\n')

# if __name__ == "__main__":
#     directory = "/alf/data/nstrangmann/CalorimeTree/Input"  # Replace with the path to your directory
#     pattern = "_101_Charged.root"
#     output_file = "output.txt"  # Output file name
#     find_files(directory, pattern, output_file)


def clear_logs():
    # Delete all existing log files before starting new jobs
    for root, dirs, files in os.walk("./"):
        for file in files:
            if file.endswith("CutsAnalysis.log"):
                os.remove(os.path.join(root, file))

def hadd_root_files(input_folder, output_file):
    # Collect all ROOT files in the given folder
    root_files = [f for f in os.listdir(input_folder) if (f.startswith('HistosFromTree_') and f.endswith('.root'))]
    if not root_files:
        raise RuntimeError("No ROOT files found in the specified folder.")
    
    # Create the full paths for the ROOT files
    root_files = [os.path.join(input_folder, f) for f in root_files]
    
    # Construct the hadd command
    command = ['hadd', '-f', output_file ] + root_files
    
    # Execute the hadd command
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        log.error("Error running hadd:")
        log.error(result.stderr.decode('utf-8'))
        raise RuntimeError("hadd failed")
    else:
        print(f"Successfully combined ROOT files into {output_file}")

        # Delete the original ROOT files
        for root_file in root_files:
            try:
                os.remove(root_file)
                # log.info(f"Want to delete file: {root_file}")
            except OSError as e:
                log.error(f"Error deleting file {root_file}: {e}")

def run_macro(dataset, setting, cut, nSplit, task_id, progress, runOptions):
    processes = []
    progress_values = [0] * nSplit  # Initialize progress list for all jobs
    
    doIsoGamma = runOptions['doIsoGamma']
    doJets = runOptions['doJets']
    doGGPi0 = runOptions['doGGPi0']
    domPi0 = runOptions['domPi0']
    doPlotting = runOptions['doPlotting']
    doAnalysisExclGammaJet = runOptions['doAnalysisExclGammaJet']
    doPlottingExclGammaJet = runOptions['doPlottingExclGammaJet']

    def monitor_progress():
        while any(process.poll() is None for _, process in processes):
            for iJob, process in processes:
                log_file = f"{dataset}/{setting}/{cut}/log_{iJob}_CutsAnalysis.log"
                if os.path.exists(log_file):
                    with open(log_file, 'r') as f:
                        lines = f.readlines()
                        for line in reversed(lines):
                            match = re.search(r'\[(\d+(\.\d+)?)%\]', line)
                            if match:
                                progress_value = float(match.group(1))
                                progress_values[iJob - 1] = progress_value
                                break
            min_progress = min(progress_values)
            progress.update(task_id, completed=min_progress)
            time.sleep(0.5)  # Adjust the sleep time as needed

    for iJob in range(1, nSplit + 1):
        command = f'srun --partition=short --job-name=ct_{iJob} --output={dataset}/{setting}/{cut}/log_{iJob}_CutsAnalysis.log root -b -q ./Analysis/makeHistosFromTree\_C.so\(\\"{dataset}/{setting}/{cut}\\"\,\{iJob}\)'
        # if none of the analysis options are enables, just run a dummy job
        if not (doIsoGamma or doJets or doGGPi0 or domPi0):
            command = f'echo "Dummy job {iJob}. Nothing to see here!" >> {dataset}/{setting}/{cut}/log_{iJob}_CutsAnalysis.log'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        processes.append((iJob, process))
        if (doIsoGamma or doJets or doGGPi0 or domPi0):
            time.sleep(0.2)

    progress_thread = threading.Thread(target=monitor_progress)
    progress_thread.start()
    # Wait for all subprocesses to complete
    for _, process in processes:
        process.wait()

    # Ensure the progress is marked as 100% upon completion
    progress.update(task_id, completed=100)
    progress_thread.join()  # Wait for the progress monitoring thread to complete
    if doIsoGamma or doJets or doGGPi0 or domPi0:
        hadd_root_files(f'{dataset}/{setting}/{cut}', f'{dataset}/{setting}/{cut}/HistosFromTree.root')
    if doPlotting:
        command = f'srun --partition=vip --job-name=ctp --output={dataset}/{setting}/{cut}/log_plot_CutsAnalysis.log root -b -q ./Analysis/plotHistosFromTree\_C.so\(\\"{dataset}/{setting}/{cut}\\"\)'
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            log.error(result.stderr)
            raise RuntimeError("plotHistosFromTree.C failed")
    if doAnalysisExclGammaJet:
        command = f'srun --partition=vip --job-name=ctp --output={dataset}/{setting}/{cut}/log_AnalysisExclGammaJet.log root -b -q ./Analysis/analyseExclGammaJet\_C.so\(\\"{dataset}/{setting}/{cut}\\"\)'
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            log.error(result.stderr)
            raise RuntimeError("analyseExclGammaJet.C failed")
    if doPlottingExclGammaJet:
        command = f'srun --partition=vip --job-name=ctp --output={dataset}/{setting}/{cut}/log_PlottingExclGammaJet.log root -b -q ./Analysis/plotExclGammaJet\_C.so\(\\"{dataset}/{setting}/{cut}\\"\)'
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            log.error(result.stderr)
            raise RuntimeError("plotExclGammaJet.C failed")

def run_multiple_macros(jobs , runOptions):
    color_cycle = cycle(COLORS)  # Create a cycle iterator for colors
    trainconfig_colors = {}
    console = Console()
    bar_columns = {}

    # Initialize a single Progress object
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        TimeRemainingColumn(),
        console=console
    ) as progress:
        tasks = {}
        processes = []

        for job in jobs:
            dataset, trainconfig, cut, nSplit = job[:4]  # Ensure we only unpack the expected number of elements
            nSplit = int(nSplit)  # Ensure nSplit is an integer
            description = f"{dataset}, {trainconfig}, {cut}"

            # Temporarily update progress columns to include the custom BarColumn
            progress.columns = (
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TimeRemainingColumn(),
            )

            # Assign a new color if not already assigned
            if trainconfig not in trainconfig_colors:
                color = next(color_cycle)
                style = f"[{color}]"
                trainconfig_colors.update({trainconfig : style})

            # Create a new task with a styled description
            task_id = progress.add_task(f"{trainconfig_colors[trainconfig]}{description}[/]", total=100)
            tasks[description] = task_id

            process = Process(target=run_macro, args=(dataset, trainconfig, cut, nSplit, task_id, progress, runOptions))
            processes.append(process)
            process.start()

        for process in processes:
            process.join()


def check_files_need_distributing(subdirectory, nSplit):
    max_i = -1
    pattern = re.compile(r'_group_(\d+)')

    # Iterate through the files in the given subdirectory
    for filename in os.listdir(subdirectory):
        match = pattern.search(filename)
        if match:
            i = int(match.group(1))
            if i > max_i:
                max_i = i

    return max_i != nSplit

def delete_existing_grouped_files(subdirectory):
    pattern = "_group_"

    # Iterate through the files in the given subdirectory
    for filename in os.listdir(subdirectory):
        if pattern in filename:
            file_path = os.path.join(subdirectory, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
                log.info(f"Deleted file: {file_path}")

def get_file_size(file_path):
    if os.path.isfile(file_path):
        return os.path.getsize(file_path)
    else:
        return 0

def distribute_files_evenly(file_groups, num_output_files):
    # Sort groups by the size of the main files (excluding _histos_)
    file_groups.sort(key=lambda x: get_file_size(x[0]), reverse=True)
    
    output_files = [[] for _ in range(min(num_output_files, len(file_groups)))]
    output_sizes = [0] * len(output_files)
    
    for main_file, histos_file in file_groups:
        smallest_file_idx = output_sizes.index(min(output_sizes))
        output_files[smallest_file_idx].append(main_file)
        if histos_file:
            output_files[smallest_file_idx].append(histos_file)
        output_sizes[smallest_file_idx] += get_file_size(main_file)
    
    return output_files

def distribute_files(analysisdirectory, inputdatapath, trainconfig, num_output_files):
    analysisdirectory+='/InputFiles'
    if check_files_need_distributing(analysisdirectory, num_output_files):
        log.info(f"Distributing files for {analysisdirectory}...")
        delete_existing_grouped_files(analysisdirectory)
    else:
        log.info(f"Files in {analysisdirectory} do not need to be redistributed")
        return

    input_file=inputdatapath + '/InputFiles_GammaIsoTree_' + trainconfig + '.txt'

    input_path = pathlib.Path(input_file).resolve()
    input_directory = input_path.parent
    
    try:
        with open(input_path, 'r') as f:
            file_paths = [line.strip() for line in f.readlines()]
    except FileNotFoundError:
        log.fatal(f"The file {input_file} does not exist.")
    
    # Group files by their main part (excluding _histos_)
    file_dict = defaultdict(lambda: [None, None])
    for path in file_paths:
        if "_histos_" in path:
            main_path = path.replace("_histos_", "_")
            file_dict[main_path][1] = path
        else:
            file_dict[path][0] = path
    
    # Create list of pairs (main file, histos file or None)
    file_groups = [(main_file, histos_file) for main_file, histos_file in file_dict.values() if main_file]
    
    if len(file_groups) < num_output_files:
        log.warning(f"Warning: Number of input file groups ({len(file_groups)}) is less than the number of requested output files ({num_output_files}). Creating {len(file_groups)} output files instead.")
        num_output_files = len(file_groups)
    
    distributed_files = distribute_files_evenly(file_groups, num_output_files)
    
    for i, file_group in enumerate(distributed_files):
        output_file_name = f"InputFiles_group_{i+1}{input_path.suffix}"
        output_file_path = f"{analysisdirectory}/{output_file_name}"
        with open(output_file_path, 'w') as f:
            for path in file_group:
                f.write(f"{path}\n")
    
    log.info(f"Distributed paths into {num_output_files} files in the directory: {analysisdirectory}")

def compile_makeHistosFromTree():
    command = 'root -q -b -x ./Analysis/makeHistosFromTree.C+\(\\"\\"\,\-1\)'
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # if compilation fails, raise an error and print the whole compilation output
    if result.returncode != 0:
        log.error("Error compiling makeHistosFromTree.C:")
        log.error("Error message reads:")
        log.exception(result.stderr)
        raise RuntimeError("Compilation failed")
def compile_plotHistosFromTree():
    command = 'root -q -b -x ./Analysis/plotHistosFromTree.C+\(\\"\\"\)'
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # if compilation fails, raise an error and print the whole compilation output
    if result.returncode != 0:
        log.error("Error compiling plotHistosFromTree.C:")
        log.error("Error message reads:")
        log.exception(result.stderr)
        raise RuntimeError("Compilation failed")
def compile_analyseExclGammaJet():
    command = 'root -q -b -x ./Analysis/analyseExclGammaJet.C+\(\\"\\"\)'
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # if compilation fails, raise an error and print the whole compilation output
    if result.returncode != 0:
        log.error("Error compiling analyseExclGammaJet.C:")
        log.error("Error message reads:")
        log.exception(result.stderr)
        raise RuntimeError("Compilation failed")
def compile_plotExclGammaJet():
    command = 'root -q -b -x ./Analysis/plotExclGammaJet.C+\(\\"\\"\)'
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # if compilation fails, raise an error and print the whole compilation output
    if result.returncode != 0:
        log.error("Error compiling plotExclGammaJet.C:")
        log.error("Error message reads:")
        log.exception(result.stderr)
        raise RuntimeError("Compilation failed")

def check_and_create_folder(folder_path):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        # If the folder does not exist, create it
        os.makedirs(folder_path)
        log.info(f"Folder '{folder_path}' created.")
    else:
        log.info(f"Folder '{folder_path}' already exists.")

def run_debug(doPlotting, doAnalysisExclGammaJet, doPlottingExclGammaJet):
    process = subprocess.run('root -b -q -x ./Analysis/makeHistosFromTree.C\(\\"DummyDataSet/DummyTrainConfig/Standard\\"\,\\0\\)', shell=True)
    process = subprocess.run('mv DummyDataSet/DummyTrainConfig/Standard/HistosFromTree_0.root DummyDataSet/DummyTrainConfig/Standard/HistosFromTree.root', shell=True)
    if doPlotting:
        process = subprocess.run('root -b -q -l ./Analysis/plotHistosFromTree.C\(\\"DummyDataSet/DummyTrainConfig/Standard\\"\\)', shell=True, check=True)
    if doAnalysisExclGammaJet:
        process = subprocess.run('root -b -q -l ./Analysis/analyseExclGammaJet.C\(\\"DummyDataSet/DummyTrainConfig/Standard\\"\\)', shell=True, check=True)
    if doPlottingExclGammaJet:
        process = subprocess.run('root -b -q -l ./Analysis/plotExclGammaJet.C\(\\"DummyDataSet/DummyTrainConfig/Standard\\"\\)', shell=True, check=True)

# Function to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Script for processing datasets.')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    parser.add_argument('--config', type=str, help='Path to the configuration file, defaults to RunConfig.yaml', default='RunConfig.yaml')
    return parser.parse_args()

# Main function
def main():

    args = parse_args()
    doDebug = args.debug

    # Ensure ROOT environment is loaded
    load_root_environment()

    analysis_config = read_yaml(args.config)

    doPlotting = analysis_config.get('doPlotting')
    doIsoGamma = analysis_config.get('doIsoGamma')
    doJets = analysis_config.get('doJets')
    doGGPi0 = analysis_config.get('doGGPi0')
    domPi0 = analysis_config.get('domPi0')
    doAnalysisExclGammaJet = analysis_config.get('doAnalysisExclGammaJet')
    doPlottingExclGammaJet = analysis_config.get('doPlottingExclGammaJet')

    log.info("Starting the analysis...")
    log.info("Detected options:")
    log.info(f"doIsoGamma: {doIsoGamma} | doJets: {doJets} | doGGPi0: {doGGPi0} | domPi0: {domPi0} | doPlotting: {doPlotting} | doAnalysisExclGammaJet: {doAnalysisExclGammaJet} | doPlottingExclGammaJet: {doPlottingExclGammaJet}")

    
    
    if doDebug:
        run_debug(doPlotting, doAnalysisExclGammaJet, doPlottingExclGammaJet)
        exit()

    cuts_config = read_yaml('Cuts.yaml')
    # count number of lines in file "./Analysis/makeHistosFromTree.C"
    # lines = 0
    # with open("./Analysis/makeHistosFromTree.C") as f:
    #     for line in f:
    #         lines += 1
    # log.info("Be patient, I am busy compiling {} lines of code ...".format(lines))
    log.info("Compiling makeHistosFromTree.C...")
    compile_makeHistosFromTree()
    log.info("Compiling plotHistosFromTree.C...")
    compile_plotHistosFromTree()
    log.info("Compiling analyseExclGammaJet.C...")
    compile_analyseExclGammaJet()
    log.info("Compiling plotExclGammaJet.C...")
    compile_plotExclGammaJet()

    nSplit = analysis_config.get('nParallelJobsPerVar', 1)
    log.info(f"Running {nSplit} parallel jobs per variation.")

     # Clear all existing log files
    clear_logs()

    # Keep track of processed MC productions to ensure distribute_files runs only once per MC production
    processed_settings = set()
    jobs = []

    for dataset, settings in analysis_config.items():
        if not isinstance(settings, dict):
            continue
        if 'trainconfigs' not in settings:
            continue
        # Now we know this is actually a dataset
        inputdatapath = settings.get('path', None)
        for trainconfig, cut in settings['trainconfigs'].items():
            trainconfigdir = f'{dataset}/{trainconfig}'
            if (dataset, trainconfig) not in processed_settings and cut != 'disabled':
                if not inputdatapath:
                    log.warning(f"No path specified for dataset {dataset}. Skipping.")
                    continue
                check_and_create_folder(f'{trainconfigdir}/InputFiles')
                distribute_files(trainconfigdir, inputdatapath, trainconfig, nSplit)
                processed_settings.add((dataset, trainconfig))
            
            if cut == 'disabled':
                continue
            # elif cut == 'Standard':
            #     check_and_create_folder(f'{trainconfigdir}/{cut}')
            #     jobs.append((dataset, trainconfig, cut, nSplit))
            elif cut == 'full':
                for cut_name in cuts_config.keys():
                    check_and_create_folder(f'{trainconfigdir}/{cut_name}')
                    jobs.append((dataset, trainconfig, cut_name, nSplit))
            else:
                # split string comma separated to extract multiple cuts
                try:
                    cuts = cut.split(',')
                # if it does not work for any reason, just continue
                except:
                    log.warning(f"Could not split cuts for {dataset}/{trainconfig}. Skipping.")
                    continue
                for cut_name in cuts:
                    # remove all white spaces from cut_name
                    cut_name = cut_name.replace(" ", "")
                    # check if cut name is in cuts_config.keys(). If it is not, if one found an entry that is onlye different by one character ask "Did you mean ..."
                    foundCut = False
                    similarCut = ""
                    simScore = 0
                    for cut_name_config in cuts_config.keys():
                        if cut_name_config.lower() == cut_name.lower():
                            cut_name = cut_name_config
                            foundCut = True
                            break
                        # if you found a match that is similar (only one character different) ask if the user meant that
                        else:
                            s = SequenceMatcher(None, cut_name_config.lower(), cut_name.lower()).ratio()
                            if s > simScore:
                                similarCut = cut_name_config
                                simScore = s
                    if not foundCut:
                        log.error(f'Cut "{cut_name}" not found in Cuts.yaml. Skipping.')
                        if similarCut:
                            log.error(f'Did you maybe mean "{similarCut}" (y or n)?')
                            # wait for user interaction
                            usrInp = input("(y or n)")
                            if usrInp == "y":
                                cut_name = similarCut
                                foundCut = True
                            else:
                                continue  
                    if foundCut:
                        check_and_create_folder(f'{trainconfigdir}/{cut_name}')
                        jobs.append((dataset, trainconfig, cut_name, nSplit))
    runOptions = { 'doIsoGamma': doIsoGamma, 'doJets': doJets, 'doGGPi0': doGGPi0, 'domPi0': domPi0, 'doPlotting': doPlotting, 'doAnalysisExclGammaJet': doAnalysisExclGammaJet, 'doPlottingExclGammaJet': doPlottingExclGammaJet }
    run_multiple_macros(jobs,runOptions)


if __name__ == "__main__":
    main()
