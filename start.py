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

def load_root_environment():
    """Ensure the ROOT environment is loaded."""
    try:
        # Check if ROOT is already loaded
        subprocess.run(['root-config', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("ROOT environment is already loaded.")
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Load ROOT environment
        root_setup_script = "/software/nstrangmann/root/bin/thisroot.sh"  # Update this path
        if os.path.exists(root_setup_script):
            command = f"source {root_setup_script} && root-config --version"
            # Use a subshell to source the script and check ROOT version
            result = subprocess.run(command, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                raise RuntimeError("Failed to load ROOT environment.")
            print("ROOT environment loaded successfully.")
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


def clear_logs(directory):
    # Delete all existing log files before starting new jobs
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".log"):
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
        print("Error running hadd:")
        print(result.stderr)
        raise RuntimeError("hadd failed")
    else:
        print(f"Successfully combined ROOT files into {output_file}")

        # Delete the original ROOT files
        for root_file in root_files:
            try:
                os.remove(root_file)
                print(f"Deleted file: {root_file}")
            except OSError as e:
                print(f"Error deleting file {root_file}: {e}")

def run_macro(directory, dataset, setting, cut, nSplit, task_id, progress):
    processes = []
    progress_values = [0] * nSplit  # Initialize progress list for all jobs

    def monitor_progress():
        while any(process.poll() is None for _, process in processes):
            for iJob, process in processes:
                log_file = f"{directory}/{dataset}/{setting}/{cut}/log_{iJob}.log"
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
        # command = f'srun --job-name=ct_{iJob} --output={directory}/{dataset}/{setting}/{cut}/log_{iJob}.log "root -q -b -x gSystem->Load("./Analysis/makeHistosFromTree_C.so"); makeHistosFromTree_C.so\(\\"{directory}/{dataset}/{setting}/{cut}\\"\,\{iJob}\)"'
        command = f'srun --job-name=ct_{iJob} --output={directory}/{dataset}/{setting}/{cut}/log_{iJob}.log root -b -q -l ./Analysis/makeHistosFromTree.C\(\\"{directory}/{dataset}/{setting}/{cut}\\"\,\{iJob}\)'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        processes.append((iJob, process))

        time.sleep(0.2)

    progress_thread = threading.Thread(target=monitor_progress)
    progress_thread.start()
    progress_thread.join()  # Wait for the progress monitoring thread to complete

    hadd_root_files(f'{directory}/{dataset}/{setting}/{cut}', f'{directory}/{dataset}/{setting}/{cut}/HistosFromTree.root')

def run_multiple_macros(directory, jobs):
    color_cycle = cycle(COLORS)  # Create a cycle iterator for colors

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
            dataset, setting, cut, nSplit = job[:4]  # Ensure we only unpack the expected number of elements
            nSplit = int(nSplit)  # Ensure nSplit is an integer
            description = f"{dataset}, {setting}, {cut}"

            # Assign a new color if not already assigned
            color = next(color_cycle)
            style = f"[{color}]"

            # Temporarily update progress columns to include the custom BarColumn
            progress.columns = (
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TimeRemainingColumn(),
            )

            # Create a new task with a styled description
            task_id = progress.add_task(f"{style}{description}[/]", total=100)
            tasks[description] = task_id

            thread = threading.Thread(target=run_macro, args=(directory, dataset, setting, cut, nSplit, task_id, progress))
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
                print(f"Deleted file: {file_path}")

def get_file_size(file_path):
    return os.path.getsize(file_path)

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

def distribute_files(subdirectory, num_output_files):
    subdirectory+='/InputFiles'
    if check_files_need_distributing(subdirectory, num_output_files):
        log.info(f"Distributing files for {subdirectory}...")
        delete_existing_grouped_files(subdirectory)
    else:
        log.info(f"Files in {subdirectory} do not need to be redistributed")
        return

    input_file=subdirectory + '/InputFiles.txt'

    input_path = pathlib.Path(input_file).resolve()
    input_directory = input_path.parent
    
    try:
        with open(input_path, 'r') as f:
            file_paths = [line.strip() for line in f.readlines()]
    except FileNotFoundError:
        print(f"Error: The file {input_file} does not exist.")
        return
    
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
        print(f"Warning: Number of input file groups ({len(file_groups)}) is less than the number of requested output files ({num_output_files}). Creating {len(file_groups)} output files instead.")
        num_output_files = len(file_groups)
    
    distributed_files = distribute_files_evenly(file_groups, num_output_files)
    
    for i, file_group in enumerate(distributed_files):
        output_file_name = f"{input_path.stem}_group_{i+1}{input_path.suffix}"
        output_file_path = input_directory / output_file_name
        with open(output_file_path, 'w') as f:
            for path in file_group:
                f.write(f"{path}\n")
    
    print(f"Distributed paths into {num_output_files} files in the directory: {input_directory}")

def compile_makeHistosFromTree():
    command = 'root -q -b -x ./Analysis/makeHistosFromTree.C+\(\\"\\"\,\-1\)'
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print("Error compiling makeHistosFromTree.C:")
        print(result.stderr)
        raise RuntimeError("Compilation failed")

def check_and_create_folder(folder_path):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        # If the folder does not exist, create it
        os.makedirs(folder_path)
        log.info(f"Folder '{folder_path}' created.")
    else:
        log.info(f"Folder '{folder_path}' already exists.")


# Main function
def main(analysis_directory):

    # Ensure ROOT environment is loaded
    load_root_environment()

    compile_makeHistosFromTree()

    # Load configurations
    analysis_config_path = os.path.join(analysis_directory, 'RunConfig.yaml')
    cuts_config_path = os.path.join(analysis_directory, 'Cuts.yaml')

    analysis_config = read_yaml(analysis_config_path)
    cuts_config = read_yaml(cuts_config_path)

    nSplit = analysis_config.get('nParallelJobsPerVar', 1)
    print(f"Running {nSplit} parallel jobs per variation")

     # Clear all existing log files
    clear_logs(analysis_directory)

    # Keep track of processed MC productions to ensure distribute_files runs only once per MC production
    processed_settings = set()
    jobs = []

    for dataset, settings in analysis_config.items():
        if not isinstance(settings, dict):
            continue
        if 'trainconfigs' not in settings:
            continue
        # Now we know this is actually a dataset
        for trainconfig, cut in settings['trainconfigs'].items():
            trainconfigdir = f'{analysis_directory}/{dataset}/{trainconfig}'
            if (dataset, trainconfig) not in processed_settings and cut != 'disabled':
                distribute_files(trainconfigdir, nSplit)
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

    run_multiple_macros(analysis_directory, jobs)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py <AnalysisDirectory>")
        sys.exit(1)

    analysis_directory = sys.argv[1]
    main(analysis_directory)
