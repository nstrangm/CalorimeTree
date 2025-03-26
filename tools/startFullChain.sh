#!/bin/bash
outputname=JE_LHC23_PbPb_pass4_calo
outputTreeName=0000014_LHC23_PbPb_pass4_calo
# Optional, download trees directly
# python3 downloadHyperloop.py --inputfilelist=inputfiles.txt --outputfolder=/alf/data/calorimetrees/je_derived_data/${outputname} --filename=AO2D.root,AnalysisResults.root --nThreads=7

# Run Gamma-Jet Tree producer task
# python3 runGammaJetTreeProducer.py --input=/alf/data/calorimetrees/je_derived_data/${outputname}/fullInputList_AO2D.txt --output=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName} --configuration=/software/flo/CalorimeTree/tools/LocalRunConfigs/PbPb_JEDerived --nFilesPerJob=10 --partition=long

# check for broken jobs and resubmit them
# python3 findCrashedJobsAndResubmit.py --inputfolder=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName}

# hadd -f -k /alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName}/AnalysisResults_merged.root $(find /alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName} -type f -name "AnalysisResults.root")

# Run conversion to appropriate format
# python3 startConversion.py --input=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName} --nFilesPerJob=2 

# Split into different centrality classes
python3 RunEventSorting.py --inputfiles=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName}/converted/InputFiles_GammaIsoTree_Run3.txt --output=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName}/sorted

# AO2D merger test
# o2-aod-merger --input /software/flo/CalorimeTree/LocalPbPbTestFlorian_GammaJetTree/inputfiles.txt --output AO2D_merged_GammaJetTree.root  --max-size 1

# cleanup
# Cleanup: Delete everything in outputTreeName folder except /sorted and /converted directories
# echo "Cleaning up ${outputTreeName} directory..."
# find /alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName} -mindepth 1 -not -path "*/sorted*" -not -path "*/converted*" -delete
# echo "Cleanup completed."
