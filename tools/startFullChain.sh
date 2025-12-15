#!/bin/bash
folder=JE_OO_AO2D
outputname=JE_LHC24ar_pass1_PbPbRef
outputTreeName=451418_LHC25ae_OO_cpass0_calo_IsoR02
# python3 downloadHyperloop.py --inputfilelist=inputfiles.txt --outputfolder=/alf/data/calorimetrees/je_derived_data/${outputname} --filename=AO2D.root,AnalysisResults.root --nThreads=7

# optional, download trees directlz
# python3 downloadHyperloop.py --inputfilelist=inputfiles.txt --outputfolder=/alf/data/calorimetrees/${folder}/${outputTreeName} --filename=AO2D.root,AnalysisResults.root --nThreads=5


# Run Gamma-Jet Tree producer task
# python3 runGammaJetTreeProducer.py --input=/alf/data/calorimetrees/je_derived_data/${outputname}/fullInputList_AO2D.txt --output=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName} --configuration=/software/flo/CalorimeTree/tools/LocalRunConfigs/PbPb_JEDerived --nFilesPerJob=10 --partition=long

# check for broken jobs and resubmit them
# python3 findCrashedJobsAndResubmit.py --inputfolder=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName}

hadd -f -k /alf/data/calorimetrees/${folder}/${outputTreeName}/AnalysisResults_merged.root $(find /alf/data/calorimetrees/${folder}/${outputTreeName} -type f -name "AnalysisResults.root")

# Run conversion to appropriate format
# python3 startConversion.py --input=/alf/data/calorimetrees/${folder}/${outputTreeName} --nFilesPerJob=1 

# Split into different centrality classes (skip for pp)
# python3 RunEventSorting.py --inputfiles=/alf/data/calorimetrees/${folder}/${outputTreeName}/converted/InputFiles_GammaIsoTree_Run3.txt --output=/alf/data/calorimetrees/${folder}/${outputTreeName}/sorted --centralities=0-10,0-30,10-30,30-50,50-90,0-100

# AO2D merger test
# o2-aod-merger --input /software/flo/CalorimeTree/LocalPbPbTestFlorian_GammaJetTree/inputfiles.txt --output AO2D_merged_GammaJetTree.root  --max-size 1

# cleanup
# Cleanup: Delete everything in outputTreeName folder except /sorted and /converted directories
# echo "Cleaning up ${outputTreeName} directory..."
# find /alf/data/calorimetrees/JE_PbPb_AO2D/${outputTreeName} -mindepth 1 -not -path "*/sorted*" -not -path "*/converted*" -delete
# echo "Cleanup completed."
