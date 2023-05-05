import os
import subprocess
import scipy.io
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import sys
import Building_MI_matrices
import Inference_algo_module as Inf_algo_mod
import Precision_Recall_module as PR_module

#Constants
#min_bins_or_neighbors = 2 ### 1 is for k and 2 is for FB
#max_bins_or_neighbors = 20 # 10 if or kNN and 20 for FB
#binning_rule = "KNN3" #"Sqrt" #or "Sturges"
binning_rule = "Sqrt" #or "Sturges"

#binning_rule = "Sturges"

#mi_est = "Shannon" #or "Miller-Madow" #or "KSG"
mi_est = "Miller-Madow"
#mi_est = "KSG"
#min_bins_or_neighbors = 1 # 2 ### 1 is for k and 2 is for FB
#max_bins_or_neighbors = 15 #20 # 10 if or kNN and 20 for FB
#infer_algo = "RL" # "CLRvMinet_pandas"
#infer_algo = "CMIA_CLR_vKSG" # "CLRvMinet_pandas"
#infer_algo = "CLRvMinet_pandas"
#infer_algo = "ARACNE"
infer_algo = "SA_CLR_v2"

#expression_data_type = "SS"
expression_data_type = "knockouts"

file = "./DATA/" +"insilico_size100_1_knockouts.tsv"
topology = "insilico_size100_1"
replicate = 1
### Importing data and calculating number of bins to use based on some rule (for binning methods)
input1_data_array = Building_MI_matrices.import_and_clean_input_file(file)
file = topology +"_" + expression_data_type + "_all.tsv"
bins_or_neighbors = Building_MI_matrices.calc_bins_num(input1_data_array.shape[1],binning_rule)

### open AUPR output file for summary statistics
output_filename = "AUPR_insilico_size100_1_knockouts_"+ expression_data_type + "_" + mi_est + "_" + infer_algo + ".dat"

# file is empty or doesn't exist
output_file = open(output_filename,"w")
#AUPR_bins_array = np.zeros(len(range(min_bins_or_neighbors,max_bins_or_neighbors + 1)), dtype=float)
AUPR_bins_array = np.zeros(len(range(1)), dtype=float)

#for idx,bins_or_neighbors in enumerate(range(max_bins_or_neighbors,min_bins_or_neighbors - 1,-1)):
idx = 0 # if not using the for loop above
### Function to calc MI2 and save to file
if mi_est == "Shannon":
    mi_est_string = "FB" + str(bins_or_neighbors) + "_Shan"
elif mi_est == "Miller-Madow":
    mi_est_string = "FB" + str(bins_or_neighbors) + "_MM"
elif mi_est == "KSG":
    mi_est_string = "KNN" + str(max_bins_or_neighbors) + "_KSG" ### kNN is calculated for all k=1..max_bins...

MI2_filename = file[:-4] + "_MI2_" + mi_est_string + ".mat"
TC_filename = file[:-4] + "_TC_" + mi_est_string + ".mat"

if MI2_filename in ''.join(os.listdir()):
    print("MI2 mat file in folder")
else:
    #print("need to build MI2 mat file") # Debug
    Building_MI_matrices.MI2_matrix_build(file,input1_data_array,bins_or_neighbors,mi_est)

if TC_filename in ''.join(os.listdir()):
    print("TC mat file in folder")
else:
    if "FB" in mi_est_string:
        Building_MI_matrices.TC_FB_matrix_build_from_entropies(file,input1_data_array,bins_or_neighbors,mi_est)
    elif "KSG" in mi_est_string:
        Building_MI_matrices.TC_KSG_matrix_build(file,input1_data_array,bins_or_neighbors,mi_est)

if mi_est == "KSG":
    mi_est_string = "KNN" + str(bins_or_neighbors) + "_KSG"

MI2_DB_filename = file[:-4] + "_MI2_" + mi_est_string + ".dat"

if MI2_DB_filename in ''.join(os.listdir()):
    print("MI2_DB_filename already exists")
else:
    [MI2_matrix] = Building_MI_matrices.building_DB_wMI2only(MI2_filename,bins_or_neighbors,mi_est)

TC_DB_filename = file[:-4] + "_MI2andTC_" + mi_est_string + ".dat"

if TC_DB_filename in ''.join(os.listdir()):
    print("TC_DB_filename already exists")
else:
    [MI2_matrix, TC_matrix] = Building_MI_matrices.building_DB_wMI2andTC(TC_filename,bins_or_neighbors,mi_est)

### PR & AUPR
#mi_est_full = [bins_or_neighbors, mi_est]
#print(PR_module.PR_per_infer_algo(topology,expression_data_type,mi_est_full,infer_algo)) # Debug
#AUPR_bins_array[idx] = PR_module.AUPR_calc(PR_module.PR_per_infer_algo(topology,expression_data_type,mi_est_full,infer_algo))

## writing output file for summary statistics
#if "-" in topology: #DREAM3 style = InSilicoSize100-Ecoli1
#    output_file.write("%.3f,%d,%s,%s,%s,%d,%s,%s\n" %(AUPR_bins_array[idx],int(topology.split("-")[0].replace("InSilicoSize",'')),topology.split("-")[1],expression_data_type,mi_est,int(bins_or_neighbors),infer_algo,replicate))
#elif "_" in topology: #DREAM4 style = insilico_size100_1
#    output_file.write("%.3f,%d,%s,%s,%s,%d,%s,%s\n" %(AUPR_bins_array[idx],int(topology.split("_")[1].replace("size",'')),topology.split("_")[2],expression_data_type,mi_est,int(bins_or_neighbors),infer_algo,replicate))

#output_file.close()
print("Done ",file)
