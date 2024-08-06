import yaml
import subprocess
import os
import re
import pathlib
import sys
from collections import defaultdict
import time
import threading
from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn, TimeElapsedColumn, SpinnerColumn
from rich.console import Console
import logging
from rich.logging import RichHandler
from rich.style import Style
from itertools import cycle

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



# Function to read YAML files
def read_yaml(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)


def clear_logs():
    # Delete all existing log files before starting new jobs
    for root, dirs, files in os.walk("./"):
        for file in files:
            if file.endswith("_prepML.log"):
                os.remove(os.path.join(root, file))



def run_macro(dataset, setting, cut, nSplit, task_id, progress):
    processes = []
    progress_values = [0] * nSplit  # Initialize progress list for all jobs

    def monitor_progress():
        while any(process.poll() is None for _, process in processes):
            for iJob, process in processes:
                log_file = f"{dataset}/{setting}/{cut}/log_{iJob}_prepML.log"
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
        command = f'srun --partition=short --job-name=ct_{iJob} --output={dataset}/{setting}/{cut}/log_{iJob}_prepML.log root -b -q -l ./ML/PrepareForML_Multiple.C\(\\"{dataset}/{setting}/{cut}\\"\,\{iJob}\)'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        processes.append((iJob, process))

        time.sleep(0.2)

    progress_thread = threading.Thread(target=monitor_progress)
    progress_thread.start()
    # Wait for all subprocesses to complete
    for _, process in processes:
        process.wait()

    # Ensure the progress is marked as 100% upon completion
    progress.update(task_id, completed=100)
    progress_thread.join()  # Wait for the progress monitoring thread to complete
    
    #if result.returncode != 0:
    #    log.error(result.stderr)
    #    raise RuntimeError("plotHistosFromTree.C failed")

def run_multiple_macros(jobs):
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
        threads = []

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

            thread = threading.Thread(target=run_macro, args=(dataset, trainconfig, cut, nSplit, task_id, progress))
            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()


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

def compile_PrepareForML_Multiple():
    log.info("Here")
    command = 'root -q -b -x ML/PrepareForML_Multiple.C+\'(\"\",-1)\''
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)#subprocess.PIPE
    if result.returncode != 0:
        raise RuntimeError("Compilation failed")
    output = result.stdout
    log.info(output)

def check_and_create_folder(folder_path):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        # If the folder does not exist, create it
        os.makedirs(folder_path)
        log.info(f"Folder '{folder_path}' created.")
    else:
        log.info(f"Folder '{folder_path}' already exists.")

# Main function
def main():

    # Ensure ROOT environment is loaded
    load_root_environment()

    analysis_config = read_yaml('PrepMLConfig.yaml')
    cuts_config = read_yaml('Cuts.yaml')

    compile_PrepareForML_Multiple()

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
            elif cut == 'Standard':
                check_and_create_folder(f'{trainconfigdir}/{cut}')
                jobs.append((dataset, trainconfig, cut, nSplit))
            elif cut == 'full':
                for cut_name in cuts_config.keys():
                    check_and_create_folder(f'{trainconfigdir}/{cut_name}')
                    jobs.append((dataset, trainconfig, cut_name, nSplit))
    
    run_multiple_macros(jobs)


if __name__ == "__main__":
    main()
