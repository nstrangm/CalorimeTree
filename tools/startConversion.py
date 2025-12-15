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
    parser.add_argument("--input", help="The input directory")
    parser.add_argument("--output", help="The output directory. If empty, then inputdir/converted will be used", default="")
    parser.add_argument("--filename", help="The input filename", default="AO2D.root")
    parser.add_argument("--nFilesPerJob", help="The number of files to convert per job", default=1)
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = args.output

    # check if input dir is specified
    if not input_dir.exists():
        log.error(f"Input directory {input_dir} does not exist")
        sys.exit(1)
    log.info(f"Using input directory {input_dir}")
    
    # if args.output is empty, then use input_dir/converted
    if output_dir == "":
        log.info(f"No output directory specified. Using {input_dir / 'converted'}")
        output_dir = input_dir / "converted"

   
    # create output directory if it does not exist
    if not output_dir.exists():
        log.info(f"Creating output directory {output_dir}")
        os.makedirs(output_dir)

    # find recursively all files in the input directory and all its subfolders
    # Find all files with the specified filename recursively in the input directory
    filename = args.filename  # Default is "AO2D.root" from the arguments
    log.info(f"Searching for files named '{filename}' in {input_dir} and its subdirectories")
    
    # Use glob to find all matching files recursively
    files = glob.glob(str(input_dir / "**" / filename), recursive=True)
    
    # Convert to Path objects for easier handling
    files = [str(Path(f)) for f in files]
    
    log.info(f"Found {len(files)} {filename} files in total")
    
    # make sure none of the paths contains the word "converted"
    files = [f for f in files if "converted" not in f]
    # for each file that is included, the same directory also has to contain a file called dpl-config.json, oderwise, skip the file
    # Filter files to only include those where the directory also contains dpl-config.json
    log.info(f"Found {len(files)} files to convert before checking for stuck jobs")
    filtered_files = []
    for file in files:
        file_path = Path(file)
        file_dir = file_path.parent
        dpl_config_path = file_dir / "dpl-config.json"
        
        filtered_files.append(file)

        # if dpl_config_path.exists():
        #     filtered_files.append(file)
        # else:
        #     log.info(f"Skipping {file} as dpl-config.json not found in the same directory")
    
    # Update the files list with only the valid ones
    files = filtered_files

    
    log.info(f"Found {len(files)} files to convert after checking for stuck jobs")
    # print all files
    for f in files:
        log.info(f"{f}")
    # exit()


    
    # split the files into chunks each containing nFilesPerJob
    nFilesPerJob = int(args.nFilesPerJob)
    filechunks = [files[i:i + nFilesPerJob] for i in range(0, len(files), nFilesPerJob)]
    log.info(f"Splitting files into {len(filechunks)} chunks with {nFilesPerJob} files each")

    # create a tmp directory to store the tmp files
    tmpdir = input_dir / "tmp"

    # create a empty file to store a list of all output files with touch
    outputfileslist = output_dir / "GammaTreeList.txt"
    os.system(f"touch {outputfileslist}")


    # remove the tmp directory if it exists
    if tmpdir.exists():
        os.system(f"rm -rf {tmpdir}")
    log.info(f"Creating tmp directory {tmpdir}")
    os.makedirs(tmpdir)
    # write each chunk to a tmp file and submit a job
    slurmJobIDs = []
    for i, filechunk in enumerate(filechunks):
        # create outputfolder in output_dir
        outputfolder = output_dir / f"{i}"
        if not outputfolder.exists():
            os.makedirs(outputfolder)


        tmpfile = tmpdir / f"tmpfile_{i}.txt"
        with open(tmpfile, "w") as f:
            for file in filechunk:
                f.write(file + "\n")

        tmpscript = tmpdir / f"RunJob_{i}.sh" 
        with open(tmpscript, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"inputfilelist={tmpfile}\n")
            f.write(f"outputfile={outputfolder / 'GammaJetTrees.root'}\n")
            f.write("root -x -q -b \"convertGammaJetRun3Tree.cpp(\\\"$inputfilelist\\\",\\\"$outputfile\\\")\"\n")
            # append the outputfile to the outputfileslist using sem with lock
            f.write("rm $inputfilelist")
        
        log.info(f"Writing chunk {i} to {tmpfile}")
        # submit a job
        log.info(f"Submitting job for chunk {i}")
        # submit a job
        cmd=f"sbatch --partition=short {tmpscript}" 
        # save the job id and run
        slurmJobIDs.append(subprocess.check_output(cmd, shell=True).decode('utf-8').split()[-1])


    filelistmakerscript = tmpdir / "filelistmaker.sh"
    with open(filelistmakerscript, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"outputdir={output_dir}\n")
        f.write(f"find {output_dir} -name \"GammaJetTrees.root\" > $outputdir/InputFiles_GammaIsoTree_Run3.txt\n")    


    # run sbatch with filelistmakerscript after all jobs from slurJobIDs are done
    cmd = f"sbatch --dependency=afterok:{':'.join(slurmJobIDs)} {filelistmakerscript}"
    log.info(f"Running {cmd}")
    os.system(cmd)
    

    

