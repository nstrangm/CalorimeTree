#Run the following command in terminal to source python: "/software/amortensen/venv/bin/activate"

#Import necessary packages:
import sys
import os
import numpy as np
import pandas as pd
import xgboost as xgb
import glob
import matplotlib.pyplot as plt
import timeit as t
from sklearn.model_selection import train_test_split
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml import plot_utils

def  main(LoadData, BuildHyperparams, TrainModel, model_file_name, OutDir, RemoveIsoVar, ApplyIsoCut, Test):

    #Set values:
    #Hyperparameter tuning:
    n_estimator_range=(200, 1000)
    max_depth_range=(2, 4)
    learning_rate_range=(0.01, 0.1)
    MaxTrails=100
    TimeOut=120
    cross_val_scoring='roc_auc'
    direction_hyperparm='maximize'
    dataratio_hyperparm=0.01

    #Model training:
    dataratio_train=0.50
    dataratio_test=0.25

    if(dataratio_hyperparm+dataratio_test+dataratio_train>1):
        print("Error! Training/Testing and Hyperparameter optimisation data overlap!")

    #Cuts and pT-binning:
    IsoCut=1.5


    if(LoadData==True):

        #Move working directory to calorimetrees folder
        os.chdir("/alf/data/calorimetrees/ML_pp_MC_AOD")
        #List of inputfiles:, MANUAL list for now. To Be Fixed...
        if(Test==True):
            InputFiles=["/alf/data/calorimetrees/ML_pp_MC_AOD/6039_20240527-1821_13TeV_JJDCal18/BASE/pTHardBin_6/294013/GammaIsoTree_100_Charged.root"]
        else:
            InputFiles = glob.glob('**/*.root',recursive=True)
            InputFiles = [i for i in InputFiles if "294925" not in i]
            InputFiles = [i for i in InputFiles if "GammaIsoTree_100_Charged.root" in i]
            InputFiles = [i for i in InputFiles if "285222" not in i]
            print("Length inputfiles")
            print(len(InputFiles))
            #print("Input Files:")
            #print(InputFiles)

        #Import data
        MCData = TreeHandler(InputFiles,'caloclustertree;1',
                              ["Event_Rho","Event_ZVtx","Event_NPrimaryTracks","Cluster_E","Cluster_Pt",
                               "Cluster_M02","Cluster_M20","Cluster_V1SplitMass","Cluster_MinMassDiffToPi0",
                               "Cluster_MinMassDiffToEta","Cluster_EFrac","Cluster_IsoGammaCorrected",
                               "Cluster_MatchTrackdEta","Cluster_MatchTrackdPhi","Cluster_MatchTrackP",
                               "Cluster_isSignal", "Cluster_SplitFloat"])#
        
        #Move working directory back to amortensen
        os.chdir("/alf/data/amortensen")

        #Separate for tuning hyper parameters and training:
        MCData_HyperParms=MCData.get_subset('Cluster_SplitFloat < '+str(dataratio_hyperparm))
        MCData_TrainTest=MCData.get_subset('Cluster_SplitFloat > '+str(1.00-(dataratio_train+dataratio_train)))
        #For tuning hyperparameters:
        signal_HyperParms=MCData_HyperParms.get_subset('Cluster_isSignal==1')
        bkg_HyperParms=MCData_HyperParms.get_subset('Cluster_isSignal==0')
        #For training and testing model:
        signal_TrainTest=MCData_TrainTest.get_subset('Cluster_isSignal==1')
        bkg_TrainTest=MCData_TrainTest.get_subset('Cluster_isSignal==0')
        #All Signal and background:
        signalH=MCData.get_subset('Cluster_isSignal==1')
        bkgH=MCData.get_subset('Cluster_isSignal==0')

        #Checking data ratios:
        print("HyperParm float cut<"+str(dataratio_hyperparm))
        print("Signal, Hyperparm:")
        print(signal_HyperParms.get_data_frame().shape[0]/signalH.get_data_frame().shape[0])
        print("Bkg, Hyperparm:")
        print(bkg_HyperParms.get_data_frame().shape[0]/bkgH.get_data_frame().shape[0])

        print("TrainTest float cut>"+str(1-(dataratio_train+dataratio_train)))
        print("Signal, TrainTest:")
        print(signal_TrainTest.get_data_frame().shape[0]/signalH.get_data_frame().shape[0])
        print("Bkg, TrainTest:")
        print(bkg_TrainTest.get_data_frame().shape[0]/bkgH.get_data_frame().shape[0])

        #Apply isolation cut if required:
        print("Check that IsoCut did something:")
        print(signal_HyperParms.get_data_frame().shape[0])
        if(ApplyIsoCut==True):
            signal_HyperParms = signal_HyperParms.get_subset('Cluster_IsoGammaCorrected < '+str(IsoCut))
            bkg_HyperParms = bkg_HyperParms.get_subset('Cluster_IsoGammaCorrected < '+str(IsoCut))
            signal_TrainTest = signal_TrainTest.get_subset('Cluster_IsoGammaCorrected < '+str(IsoCut))
            bkg_TrainTest = bkg_TrainTest.get_subset('Cluster_IsoGammaCorrected < '+str(IsoCut))
            print("Isolation cut applied to data!")
            print('Cluster_IsoGammaCorrected < '+str(IsoCut))
        else:
            print("Isolation cut not applied to data!")
        print(signal_HyperParms.get_data_frame().shape[0])


        #Generate train/test datasets:
        hyperparm_optim_data = train_test_generator([signal_HyperParms,bkg_HyperParms], [1,0], test_size=2./3., random_state=123)
        train_test_data = train_test_generator([signalH,bkgH], [1,0], test_size=dataratio_test/(dataratio_test+dataratio_train), random_state=42)

        #Do correlation plots between signal and data:
        vars_to_draw = signalH.get_var_names()
        vars_to_draw.remove('Cluster_isSignal')
        leg_labels = ['signal','background']

        plot_utils.plot_distr([signalH, bkgH], vars_to_draw, bins=100, labels=leg_labels, log=True, density=True, figsize=(12, 7), alpha=0.3, grid=False)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plot_utils.plot_corr([signalH, bkgH],vars_to_draw, leg_labels)

        if(TrainModel==False):
            save_all_figures(OutDir, 'png')

        if(TrainModel==False and BuildHyperparams==False):
            return

        #Remove quantities of interest (pT and E):
        features_for_train = vars_to_draw.copy()
        features_for_train.remove('Cluster_E')
        features_for_train.remove('Cluster_Pt')

        if(RemoveIsoVar==True):
            features_for_train.remove('Cluster_IsoGammaCorrected')
            print("Removing Isolation quantity from variables!")
        else:
            print("Isolation quantity kept in variables!")

        if(BuildHyperparams==True):
            t1=t.default_timer()
            #Set classification algorithm:
            #XGBoost is a gradient boosted decision tree.
            model_clf = xgb.XGBClassifier()
            model_hdl = ModelHandler(model_clf, features_for_train)

            #Optimize hyper parameters with optuna:
            hyper_pars_ranges = {'n_estimators': n_estimator_range, 'max_depth': max_depth_range, 'learning_rate': learning_rate_range}
            model_hdl.optimize_params_optuna(hyperparm_optim_data, hyper_pars_ranges, cross_val_scoring=cross_val_scoring, n_trials=MaxTrails,n_jobs=-1, timeout=TimeOut, direction=direction_hyperparm, show_progress_bar=True)

            #Save model:
            model_hdl.dump_model_handler(str(OutDir)+'/'+str(model_file_name))
            print(str("Time to optimize hyper params:")+str(t.default_timer()-t1))
            print("Optuna settings:")
            print("dataratio_hyperparm="+str(dataratio_hyperparm))
            print("n_estimator_range="+str(n_estimator_range))
            print("depth="+str(max_depth_range))
            print("learning_rate_range="+str(learning_rate_range))
            print("n_trails"+str(MaxTrails))
            print("Timeout="+str(TimeOut))
            print("cross_val_scoring="+str(cross_val_scoring))
            print("direction="+str(direction_hyperparm))
        if(BuildHyperparams==False):
            model_clf = xgb.XGBClassifier()
            model_hdl = ModelHandler(model_clf, features_for_train)
            model_hdl.load_model_handler(str(OutDir)+'/'+str(model_file_name))

        if(TrainModel==True):
            t2=t.default_timer()
            #Train model:
            model_hdl.train_test_model(train_test_data)
            model_hdl.dump_model_handler(str(OutDir)+'/'+str(model_file_name)+str("_Trained"))
            #plot model output from training and test datasets.
            y_pred_train = model_hdl.predict(train_test_data[0], False)
            y_pred_test = model_hdl.predict(train_test_data[2], False)
            plot_utils.plot_output_train_test(model_hdl, train_test_data)
            save_all_figures(OutDir, 'png')
            plt.close()
            print(str("Time to train:")+str(t.default_timer()-t2))
            print(str("dataratio_train=")+str(dataratio_train)) 
            print(str("dataratio_test=")+str(dataratio_test)) 

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

    if(sys.argv[6]=="True"):
        RemoveIsoVar=True
    elif(sys.argv[6]=="False"):
        RemoveIsoVar=False
    else:
        print("Wrong 6th argument!")

    if(sys.argv[7]=="True"):
        ApplyIsoCut=True
    elif(sys.argv[7]=="False"):
        ApplyIsoCut=False
    else:
        print("Wrong 7th argument!")

    if(sys.argv[8]=="True"):
        Test=True
    elif(sys.argv[8]=="False"):
        Test=False
    else:
        print("Wrong 8th argument!")

    model_file_name = str(sys.argv[4])

    OutDir = str(sys.argv[5])

    main(LoadData,BuildHyperparams,TrainModel,model_file_name,OutDir, RemoveIsoVar, ApplyIsoCut, Test)
    print("LoadData")
