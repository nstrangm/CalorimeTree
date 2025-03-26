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

# setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler()]
)

log = logging.getLogger("rich")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts a directory of files to a different format")
    parser.add_argument("--inputfiles", help="txt files with the input files")
    parser.add_argument("--output", help="The output directory", default="")
    parser.add_argument("--centralities", help="The centralities to process", default="0-10,10-30,30-50,50-90,0-90")
    args = parser.parse_args()

    input_files = Path(args.inputfiles)
    output_dir = args.output

    # check that output_dir string is not empty
    if not output_dir:
        log.error(f"Output directory not provided")
        sys.exit(1)

    # decode centralities. Each comma separated string contains a centrality range centMin and centMax
    centralities = args.centralities.split(",")
    centralities = [tuple(map(int, c.split("-"))) for c in centralities]
    log.info(f"Using centralities {centralities}")



    # create output dir
    output_dir = Path(output_dir)
    if not output_dir.exists():
        log.info(f"Creating output directory {output_dir}")
        os.makedirs(output_dir)
    log.info(f"Using input files {input_files}")

    # read all lines in input_files txt list
    with open(input_files, "r") as f:
        files = f.readlines()
    files = [f.strip() for f in files]
    log.info(f"Found {len(files)} files to process")

    # loop over files and create a folcder 0...files in output folder
    workDir = Path(__file__).resolve().parent
    slurmJobIDs = []
    for i, file in enumerate(files):
        # create folder
        folder = output_dir / str(i)
        if not folder.exists():
            log.info(f"Creating folder {folder}")
            os.makedirs(folder)
        
        # in folder create a job script
        job_script = folder / "job.sh"
        with open(job_script, "w") as f:
            f.write(f"#!/bin/bash\n")
            # write command to run root -x -q -b workDir/eventsorting.cpp
            # loop over centralities
            for centMin, centMax in centralities:
                cutstring = f"event_centrality >= {centMin} && event_centrality < {centMax}"
                outputfile = folder / f"GammaJetTree_Split_{centMin}_{centMax}.root"
                f.write(f"root -q -b \"{workDir}/eventsorting_cpp.so(\\\"{file}\\\", \\\"{outputfile}\\\", \\\"{cutstring}\\\")\"\n")
        # send job to slurm
        cmd = f"sbatch -p long {job_script}"
        slurmJobIDs.append(subprocess.check_output(cmd, shell=True).decode('utf-8').split()[-1])
        log.info(f"Submitted job for {file} with id {slurmJobIDs[-1]}")

    filemergerscript = output_dir / "merger.sh"
    with open(filemergerscript, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"outputdir={output_dir}\n")
        # loop over centralities
        for centMin, centMax in centralities:
            # find all files in outputdir with name GammaJetTree_Split_{centMin}_{centMax}.root
            f.write(f"find {output_dir} -name \"GammaJetTree_Split_{centMin}_{centMax}.root\" > $outputdir/InputFiles_GammaIsoTree_Run3-{centMin}-{centMax}.txt\n")
    cmd = f"sbatch --dependency=afterok:{':'.join(slurmJobIDs)} {filemergerscript}"
    subprocess.check_output(cmd, shell=True)



        
        



    
        




    
    

    