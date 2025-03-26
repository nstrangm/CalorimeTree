import sys
import os
import argparse
import subprocess
import threading
import logging
import time
from pathlib import Path
from rich.logging import RichHandler
import glob

# setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler()]
)

log = logging.getLogger("rich")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Finds stuck jobs and resubmits them according to the runcommand.sh")
    parser.add_argument("--inputfolder", help="Folder containing the input files")
    args = parser.parse_args()

    inputfolder = Path(args.inputfolder)

    # Check if input folder exists
    if not inputfolder.exists():
        log.error(f"Input folder {inputfolder} does not exist")
        sys.exit(1)
    
    log.info(f"Checking for stuck jobs in {inputfolder}")
    
    # List to store stuck job folders
    stuck_jobs = []
    
    # Loop over all folders in the input folder
    for folder in inputfolder.iterdir():
        if folder.is_dir():
            # Check if dpl-config.json exists in the folder
            dpl_config_file = folder / "dpl-config.json"
            if not dpl_config_file.exists():
                # Job is considered stuck
                stuck_jobs.append(folder)
                log.info(f"Found stuck job in folder: {folder}")
    
    # Print summary of stuck jobs
    if stuck_jobs:
        log.info(f"Found {len(stuck_jobs)} stuck jobs:")
        for job in stuck_jobs:
            log.info(f"  {job}")
    else:
        log.info("No stuck jobs found.")

    # Loop over all folders of stuck jobs and cd into directory and delete the AO2D.root file
    for job in stuck_jobs:
        slurm_log = job / "slurm_retry.log"
        log.info(f"Deleting AO2D.root file in {job}")
        os.system(f"cd {job} && rm AO2D.root")
        # Run runcommand.sh
        
        # Now go through input_data.txt in {job} folder and create one subfolder per file using numbers 0 ... N
        with open(job / "input_data.txt", "r") as f:
            for i, line in enumerate(f):
                file = line.strip()
                log.info(f"Creating subfolder for {file}")
                os.makedirs(job / f"{i}", exist_ok=True)
                # copy from job dir, runcommand.sh, configuration.json OutputDirector.json to subfolder
                os.system(f"cp {job}/runcommand.sh {job}/{i}/")
                os.system(f"cp {job}/configuration.json {job}/{i}/")
                os.system(f"cp {job}/OutputDirector.json {job}/{i}/")
                # in subfolder create a file called input_data.txt with the file
                with open(job / f"{i}/input_data.txt", "w") as f:
                    f.write(file)
                # submit job
                log.info(f"Submitting job {i} in {job}")
                slurm_log = job / f"{i}/slurm_retry.log"
                os.system(f"cd {job}/{i} && sbatch --partition short --output {slurm_log} runcommand.sh")
                
