#!/bin/bash
outputname=000008_LHC23zzm_JetPt0_NewTM_Occupancy_PbPb
# Run Gamma-Jet Tree producer task
# python3 runGammaJetTreeProducer.py --input=/alf/data/calorimetrees/je_derived_data/JE_LHC23zzm_pass4_Cl_8/fullInputList_AO2D.txt --output=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputname} --configuration=/software/flo/CalorimeTree/tools/LocalRunConfigs/PbPb_JEDerived --nFilesPerJob=1 --partition=long
# python3 runGammaJetTreeProducer.py --input=/alf/data/calorimetrees/je_derived_data/JE_LHC23zzm_pass4_Cl_8/fullInputList_AO2D.txt --output=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputname} --configuration=/software/flo/CalorimeTree/tools/LocalRunConfigs/PbPb_JEDerived_uniform --nFilesPerJob=1 --partition=long
python3 runGammaJetTreeProducer.py --input=/software/flo/CalorimeTree/tools/LocalRunConfigs/PbPb_AO2D/inputGridLHC23zzm.txt --output=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputname} --configuration=/software/flo/CalorimeTree/tools/LocalRunConfigs/PbPb_AO2D --nFilesPerJob=1 --partition=long

# hadd -f -k /alf/data/calorimetrees/JE_PbPb_AO2D/${outputname}/AnalysisResults_merged.root $(find /alf/data/calorimetrees/JE_PbPb_AO2D/${outputname} -type f -name "AnalysisResults.root")

# Run conversion to appropriate format
# python3 startConversion.py --input=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputname} --nFilesPerJob=1 

# Split into different centrality classes
# python3 RunEventSorting.py --inputfiles=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputname}/converted/InputFiles_GammaIsoTree_Run3.txt --output=/alf/data/calorimetrees/JE_PbPb_AO2D/${outputname}/sorted