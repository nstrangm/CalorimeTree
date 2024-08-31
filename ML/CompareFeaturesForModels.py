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

#Slightly shitty script for reading in multiple models and plotting their feature importance.
#Can only compare three at the moment.


def main():
    
    #Bar width:
    barwidth=0.2


    model_hdl1 = ModelHandler()
    model_hdl1.load_model_handler("OutDir100_Charged_main_IsoCut_01hyppar_ClusterCuts_Full/Model100Charged_main_IsoCut_01hyppar_ClusterCuts_Full_Trained")
    model_hdl2 = ModelHandler()
    model_hdl2.load_model_handler("OutDir100_Charged_main_IsoCut_01hyppar_ClusterCuts_PtTrained/Model100Charged_main_IsoCut_01hyppar_ClusterCuts_PtTrained_Trained")
    model_hdl3 = ModelHandler()
    model_hdl3.load_model_handler("OutDir100_Charged_main_IsoCut_01hyppar_HigherDepth_ClusterCuts_PtTrained/Model100Charged_main_IsoCut_01hyppar_HigherDepth_ClusterCuts_PtTrained_Trained")

    feature_important1 = model_hdl1.get_original_model().get_booster().get_score(importance_type='weight')
    #Get values:
    values1 = list(feature_important1.values())
    #Normalise:
    values1max=max(values1)
    values1 = [i/values1max for i in values1]
    values1 = [0 if i==0 else values1[i-1] for i in range(len(values1)+1)]

    feature_important2 = model_hdl2.get_original_model().get_booster().get_score(importance_type='weight')
    #Get feature names and update for plotting
    keys = list(feature_important2.keys())
    keys =  [j.replace('Cluster_', '') for j in keys]
    keys = [j.replace('Event_', '') for j in keys]
    axis = np.arange(len(keys))
    #Get values:
    values2 = list(feature_important2.values())
    #Normalise:
    values2max=max(values2)
    values2 = [i/values2max for i in values2]

    feature_important3 = model_hdl3.get_original_model().get_booster().get_score(importance_type='weight')
    #Get values:
    values3 = list(feature_important3.values())
    #Normalise:
    values3max=max(values3)
    values3 = [i/values3max for i in values3]

    print(keys)

    print(values1)
    print(values2)
    print(values3)

    plt.barh(axis-barwidth, values1, barwidth, label="Base line model")
    plt.barh(axis, values2, barwidth, label="Trained on $p_{T}$")
    plt.barh(axis+barwidth, values3, barwidth, label=r"Trained on $p_{T}$,"
                                                      "\n" 
                                                      "increased depth")
    plt.yticks(axis,keys)
    plt.yticks(rotation=45)
    plt.subplots_adjust(left=0.25)
    plt.xlabel("Feature importance/max(Feature importance)")
    plt.legend()
    plt.show()

    


if __name__ == "__main__":
    #model_file_name1 = str(sys.argv[0])
    #print(model_file_name1)
   #model_file_name2 = str(sys.argv[1])
   #model_file_name3 = str(sys.argv[2])
    main()#, model_file_name2, model_file_name3