# EXAMPLE: python3 runGammaJetTreeProducer.py --input=/alf/data/calorimetrees/je_derived_data/JE_LHC23zzm_pass4_Cl_8/fullInputList_AO2D.txt --output=/alf/data/calorimetrees/JE_PbPb_AO2D/000003_LHC23zzm_PbPb --configuration=/software/flo/CalorimeTree/tools/LocalRunConfigs/PbPb_JEDerived
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
import stat

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
SIF_PATH = SCRIPT_DIR.parent / "el9.sif"
SIF_IMAGE = "docker://alisw/slc9-builder:latest"

# setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler()]
)

log = logging.getLogger("rich")

def ensure_sif() -> None:
    if SIF_PATH.exists():
        log.info(f"SIF image found: {SIF_PATH}")
        return
    log.warning(f"SIF image not found at {SIF_PATH} – pulling {SIF_IMAGE} …")
    result = subprocess.run(["apptainer", "pull", str(SIF_PATH), SIF_IMAGE], check=False)
    if result.returncode != 0:
        log.error("Failed to pull SIF image. Is apptainer installed?")
        sys.exit(1)
    log.info("SIF pulled successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts a directory of files to a different format")
    parser.add_argument("--inputfiles", help="txt files with the input files")
    parser.add_argument("--output", help="The output directory. If empty, then inputdir/converted will be used", default="")
    parser.add_argument("--configuration", help="Folder containing the configuration files (configuration.json, OutputDirector.json and runcommand.sh)", default="")
    parser.add_argument("--nFilesPerJob", help="The number of files to convert per job", default=1)
    parser.add_argument("--partition", help="The slurm partition to use", default="short")
    args = parser.parse_args()

    input_files = Path(args.inputfiles)
    output_dir = args.output
    configuration_dir = args.configuration
    partition = args.partition

    # check that nothing empty was provided
    if not input_files.exists():
        log.error(f"Input files {input_files} does not exist")
        sys.exit(1)
    if not output_dir:
        log.error(f"Output directory not provided")
        sys.exit(1)
    if not configuration_dir or not Path(configuration_dir).exists():
        log.error(f"Configuration directory not provided")
        sys.exit(1)
    env_sh = Path(configuration_dir) / "env.sh"
    if not env_sh.exists():
        log.error(f"env.sh not found in configuration directory {configuration_dir}")
        sys.exit(1)
    log.info(f"Using input files {input_files}")
    ensure_sif()

    # create output directory if it does not exist
    output_dir = Path(output_dir)
    if not output_dir.exists():
        log.info(f"Creating output directory {output_dir}")
        os.makedirs(output_dir)

    # read all lines in input_files txt list and split them into chunks each containing nFilesPerJob
    with open(input_files, "r") as f:
        files = f.readlines()
    files = [f.strip() for f in files]
    log.info(f"Found {len(files)} files to process")
    
    # split the files into chunks each containing nFilesPerJob
    nFilesPerJob = int(args.nFilesPerJob)
    filechunks = [files[i:i + nFilesPerJob] for i in range(0, len(files), nFilesPerJob)]
    log.info(f"Splitting files into {len(filechunks)} chunks with {nFilesPerJob} files each")
    log.info(f"Using configuration directory {configuration_dir}")

    # print the content of {configuration_dir}/runcommand.sh
    with open(Path(configuration_dir) / "runcommand.sh", "r") as f:
        log.info(f"Using runcommand.sh: {f.read()}")
    
    # loop over chunks
    for i, filechunk in enumerate(filechunks):
        log.info(f"Submitting chunk {i+1}/{len(filechunks)}: ")
        for f in filechunk:
            log.info(f"{f}")

        # create folder {i} in output_dir
        output_dir_chunk = output_dir / str(i)
        if not output_dir_chunk.exists():
            os.makedirs(output_dir_chunk)
        # copy all files in configuration_dir to output_dir_chunk
        os.system(f"cp -r {configuration_dir}/* {output_dir_chunk}")

        # write files in filechunk to {output_dir_chunk}/input_data.txt
        with open(output_dir_chunk / "input_data.txt", "w") as f:
            f.write("\n".join(filechunk))

        # determine top-level directories of input files for apptainer bind mounts
        bind_dirs = set()
        for f in filechunk:
            parts = Path(f).parts
            if len(parts) >= 2:
                bind_dirs.add(f"/{parts[1]}")
        bind_args = " \\\n  ".join(
            [f"--bind /cvmfs:/cvmfs"]
            + [f"--bind {d}:{d}" for d in sorted(bind_dirs)]
            + [f"--bind {output_dir_chunk}:/workdir"]
        )

        # generate apptainer wrapper script for SLURM
        submit_sh_content = f"""#!/bin/bash
source /etc/profile.d/modules.sh
module use /software/centralsoft/Modules
module load apptainer

apptainer exec \\
  {bind_args} \\
  --cleanenv \\
  {SIF_PATH} \\
  bash -c "cd /workdir && source env.sh && bash runcommand.sh"
"""
        submit_sh = output_dir_chunk / "submit.sh"
        submit_sh.write_text(submit_sh_content)
        submit_sh.chmod(submit_sh.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

        slurm_log = output_dir_chunk / "slurm.log"
        os.system(f"sbatch --partition {partition} --output {slurm_log} {submit_sh}")


        




    
    

    

