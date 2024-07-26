#Run the following command in terminal to source python: "/software/amortensen/venv/bin/activate"

#Import necessary packages:
import sys
import os
import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml import plot_utils

def  main(LoadData, BuildHyperparams, TrainModel, model_file_name, OutDir):

    if(LoadData==True):

        #List of inputfiles:, MANUAL list for now. To Be Fixed...
        InputFiles=["/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_6/294013/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_5/289444/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_1/290459/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_7/289308/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_4/292430/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_6/294529/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_4/289444/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_1/292430/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_4/286731/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_3/286931/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_2/294529/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_2/294925/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_5/292298/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_3/286336/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_3/287654/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_3/292430/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_3/289582/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_2/286592/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_6/289444/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_7/294529/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_4/294241/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_2/294925/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_3/294586/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_6/293588/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_5/289815/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_4/286876/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_1/289370/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_4/287654/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_5/290459/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_6/289582/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_1/287480/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_4/287480/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_6/292298/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_2/294013/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_4/289582/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_2/292430/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_2/294241/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_2/294529/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_2/289582/GammaIsoTree_100_Charged.root",
                    "/alf/data/calorimetrees/ML_pp_MC_AOD/6038_20240523-1434_13TeV_GJDCal18/pTHardBin_6/293831/GammaIsoTree_100_Charged.root",]

        #Import data
        MCData = TreeHandler(InputFiles,'caloclustertree;1',
                              ["Event_Rho","Event_ZVtx","Event_NPrimaryTracks","Cluster_E","Cluster_Pt",
                               "Cluster_M02","Cluster_M20","Cluster_V1SplitMass","Cluster_MinMassDiffToPi0",
                               "Cluster_MinMassDiffToEta","Cluster_EFrac","Cluster_IsoGammaCorrected",
                               "Cluster_MatchTrackdEta","Cluster_MatchTrackdPhi","Cluster_MatchTrackP","Cluster_isSignal"])#

        #Separate in signal and background
        signalHo = MCData.get_subset('Cluster_isSignal==1')
        bkgHo = MCData.get_subset('Cluster_isSignal==0')

        #Apply Isolation cut.
        signalH = signalHo.get_subset() #signalHo.get_subset('Cluster_IsoGammaCorrected < 1.5')
        bkgH = bkgHo.get_subset()#bkgHo.get_subset('Cluster_IsoGammaCorrected < 1.5')

        #Generate train/test datasets:
        train_test_data = train_test_generator([signalH,bkgH], [1,0], test_size=0.5, random_state=42)

        #Do correlation plots between signal and data:
        vars_to_draw = signalH.get_var_names()
        vars_to_draw.remove('Cluster_isSignal')
        leg_labels = ['background','signal']

        plot_utils.plot_distr([bkgH, signalH], vars_to_draw, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plot_utils.plot_corr([bkgH, signalH],vars_to_draw, leg_labels)

        if(TrainModel==False):
            save_all_figures(OutDir, 'png')

        if(TrainModel==False and BuildHyperparams==False):
            return

        #Remove quantities of interest (pT and E):
        features_for_train = vars_to_draw.copy()
        features_for_train.remove('Cluster_E')
        features_for_train.remove('Cluster_Pt')
        #features_for_train.remove('Cluster_IsoGammaCorrected')

        if(BuildHyperparams==True):
            #Set classification algorithm:
            #XGBoost is a gradient boosted decision tree.
            model_clf = xgb.XGBClassifier()
            model_hdl = ModelHandler(model_clf, features_for_train)

            #Optimize hyper parameters with optuna:
            hyper_pars_ranges = {'n_estimators': (200, 1000), 'max_depth': (2, 4), 'learning_rate': (0.01, 0.1)}
            model_hdl.optimize_params_optuna(train_test_data, hyper_pars_ranges, cross_val_scoring='roc_auc', n_trials=100,n_jobs=-1, timeout=120, direction='maximize', show_progress_bar=True)

            #Save model:
            model_hdl.dump_model_handler(str(OutDir)+'/'+str(model_file_name))
        if(BuildHyperparams==False):
            model_clf = xgb.XGBClassifier()
            model_hdl = ModelHandler(model_clf, features_for_train)
            model_hdl.load_model_handler(str(OutDir)+'/'+str(model_file_name))

        if(TrainModel==True):
            #Train model:
            model_hdl.train_test_model(train_test_data)
            model_hdl.dump_model_handler(str(OutDir)+'/'+str(model_file_name)+str("_Trained"))
            #plot model output from training and test datasets.
            y_pred_train = model_hdl.predict(train_test_data[0], False)
            y_pred_test = model_hdl.predict(train_test_data[2], False)
            plot_utils.plot_output_train_test(model_hdl, train_test_data)
            save_all_figures(OutDir, 'png')
            plt.close()

def save_all_figures(directory='figures', file_format='png'):
    if not os.path.exists(directory):
        os.makedirs(directory)
    figures = [plt.figure(n) for n in plt.get_fignums()]
    for i, fig in enumerate(figures):
        fig.savefig(f"{directory}/figure_{i + 1}.{file_format}")


if __name__ == "__main__":
    if(sys.argv[1]=="True"):
        LoadData=True
    elif(sys.argv[1]=="False"):
        LoadData=False
    else:
        print("Wrong 1st argument!")

    if(sys.argv[2]=="True"):
        BuildHyperparams=True
    elif(sys.argv[2]=="False"):
        BuildHyperparams=False
    else:
        print("Wrong 2nd argument!")

    if(sys.argv[3]=="True"):
        TrainModel=True
    elif(sys.argv[3]=="False"):
        TrainModel=False
    else:
        print("Wrong 3rd argument!")

    model_file_name = str(sys.argv[4])

    OutDir = str(sys.argv[5])

    main(LoadData,BuildHyperparams,TrainModel,model_file_name,OutDir)
    print("LoadData")
