import glob
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

#This script takes as input the path to a trained model 
#from  PromptPhotonID_Hipe4ML_Restructure.py, applies it to data 
#and saves the output in a root file and plots Signal, Background and total BDT
#Distributions. Data used is implied to be MC-data.
#NB:
#ModelPath = Path to trained model (No file suffix is necessary)
#OutName = Name of file with dataset+BDT scores pt. entry. (Don't need ".root" suffix.)

def  main(ModelPath,OutName,OutDir,pTintegrated,Test=True,ReRead=False):
    dataratio_train=0.5
    dataratio_test=0.25
    pTbinEdges = [3,8]#[8,40]#[3,40]#[3,4,4.5,5,5.5,6,6.5,7,7.5,8,9,10,12.5,15,17.5,20,40]
    if(pTintegrated==True):
        pTbinEdges = [3,40]#[8,40]#
    #Cut="(Cluster_Pt>8)&(Cluster_IsoGammaCorrected<1.5)"
    #InputFiles=[OutDir+"/TestData.parquet"]

    #Load model
    model_hdl = ModelHandler()
    model_hdl.load_model_handler(ModelPath)

    for i in range(len(pTbinEdges)-1):
        #Get pTbin
        pTmin=pTbinEdges[i]
        pTmax=pTbinEdges[i+1]
    
    
        #Update signal, bkg and data with correct subset
        #if ReRead==False:
        dataH = TreeHandler(OutDir+"/TestData.parquet")
        print(dataH.get_data_frame().keys())
        #else:
        #    print("Reloading data!")
        #    InputFiles = glob.glob('**/*ClusterCuts.root',recursive=True)
        #    InputFiles = [i for i in InputFiles if "294925" not in i]
        #    InputFiles = [i for i in InputFiles if "GammaIsoTree_100_Charged_ClusterCuts.root" in i]
        #    InputFiles = [i for i in InputFiles if "285222" not in i]
        #    InputFiles = [i for i in InputFiles if "287249" not in i]
        #    InputFiles = [i for i in InputFiles if "6097_20240626-1006_13TeV_GJDCal18" not in i ]

        #    dataH = TreeHandler(InputFiles,'caloclustertree',
        #                      ['Cluster_E','Cluster_Pt', 'Cluster_M02', 
        #                      'Cluster_M20', 'Cluster_V1SplitMass','Cluster_MinMassDiffToPi0',
        #                      'Cluster_MinMassDiffToEta', 'Cluster_EFrac','Cluster_IsoGammaCorrected',
        #                      'Cluster_DistanceToBadChannel','Cluster_NLM','Cluster_NCells',
        #                      'Cluster_SplitFloat', 'Cluster_isSignal','Event_Weight', 'Cluster_SplitFloat'],
        #                      array_cache="inherit",
        #                      cut=Cut)
            
        #Apply Pt-cut:
        print("Before Cut")
        print(sys.getsizeof(dataH))
        dataH = dataH.get_subset('Cluster_Pt<'+str(pTmax))
        dataH = dataH.get_subset('Cluster_Pt>'+str(pTmin))
        print("Before Cut")
        print(sys.getsizeof(dataH))
        #Get data and signal:
        signalH = dataH.get_subset('Cluster_isSignal>0')
        bkgH = dataH.get_subset('Cluster_isSignal<1')
        
        # plot feautre importance
        print(model_hdl.get_original_model().feature_importances_)


        #Apply model to data:
        dataH.apply_model_handler(model_hdl, True, "BDT")
        signalH.apply_model_handler(model_hdl,True,"BDT")
        bkgH.apply_model_handler(model_hdl,True,"BDT")

        #Fearture importance:
        booster = model_hdl.get_original_model().get_booster()
        feature_names = booster.feature_names
        new_feature_names =  [j.replace('Cluster_', '') for j in feature_names]
        new_feature_names = [j.replace('Event_', '') for j in new_feature_names]
        feature_name_mapping = dict(zip(feature_names, new_feature_names))
        booster.feature_names = [feature_name_mapping.get(f, f) for f in booster.feature_names]

        #Roc Curve
        plot_utils.plot_roc(dataH["Cluster_isSignal"],dataH["BDT"])

        print("pT min/max")
        print(pTmin)
        print(pTmax)

        #Save data:
        dataH.write_df_to_root_files(base_file_name=str(OutDir)+'/'+str(OutName)+"_pT_"+str(pTmin)+"_"+str(pTmax), tree_name='caloclustertree', path='./', save_slices=False)
        selected_data_hndl = dataH.get_subset()#'model_output<-4'

        labels_list = ["Signal+Background","Signal: "+str(signalH.get_data_frame().shape[0]),"Background: "+str(bkgH.get_data_frame().shape[0])]
        colors_list = ['orangered', 'cornflowerblue','green']

        plot_utils.plot_distr([dataH,signalH,bkgH], column='BDT', bins=200, labels=labels_list, colors=colors_list, density=False,fill=True, histtype='step', alpha=0.5)#, weights=[dataH['Event_Weight'],signalH['Event_Weight'],bkgH['Event_Weight']]
        ax = plt.gca()
        ax.set_xlabel(r'BDT')
        ax.margins(x=0)
        ax.xaxis.set_label_coords(0.9, -0.075)

        fig, ax = plt.subplots()
        xgb.plot_importance(booster)
        plt.yticks(rotation=45)
        plt.subplots_adjust(left=0.2)
        save_all_figures(pTmin,pTmax,OutDir,'pdf')
    
        xgb.plot_tree(model_hdl.get_original_model(), num_trees=4, rankdir='LR')
        fig = plt.gcf()
        fig.set_size_inches(150, 100)
        fig.savefig(OutDir+str('tree.png'))


    #plt.show()

def save_all_figures(pTmin,pTmax,directory='figures', file_format='png'):
    if not os.path.exists(directory):
        os.makedirs(directory)
    figures = [plt.figure(n) for n in plt.get_fignums()]
    for i, fig in enumerate(figures):
        fig.savefig(f"{directory}/figure_LoadPlot_{i + 1}_pT_{pTmin}_{pTmax}.{file_format}")



if __name__ == "__main__":
    ModelPath = str(sys.argv[1])
    OutName = str(sys.argv[2])
    OutDir = str(sys.argv[3])

    if(str(sys.argv[4])=="True"):
        pTintegrated=True
    else:
        pTintegrated=False

    if(str(sys.argv[5])=="True"):
        Test=True
    else:
        Test=False

    if(str(sys.argv[6])=="True"):
        ReRead=True
    else:
        ReRead=False

    if(Test==False):
        OutName="Data_"+OutName

    main(ModelPath,OutName,OutDir,pTintegrated,Test,ReRead)