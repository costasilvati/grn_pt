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
import sys
argumentList = sys.argv[1:]
filename = ''
output = ''
for i in range(len(argumentList)):
        if argumentList[i] == '-i':
                filename = argumentList[i+1]
        elif argumentList[i] == '-o':
                output = argumentList[i+1]
                if output[-1] != '/':
                    output+='/'

if filename == '' or output == '':
        print('you MUST provide filename (-i <filename>) AND path out (-o <path_out>)')
        exit(1)
if filename[-3:] == 'csv':
        print('Reading ',filename,' and transforming to tsv')
        df = pd.read_csv(filename)
        df.to_csv(filename.replace('.csv','.tsv'),sep='\t')
        filename = filename.replace('.csv','.tsv')
elif filename[-3:] == 'tsv':
        #df = pd.read_csv(filename,sep='tsv')
        pass
else:
        print('The input file MUST be a .csv OR .tsv')
        exit(1)

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

### CONSIDERANDO UM NOME +- NESTE FORMATO "./DATA/" +"insilico_size100_1_knockouts.tsv"

expression_data_type = filename.split('/')[-1].split('_')[-1][:-4]  # "knockouts"

file = filename
topology = '_'.join(file.split('/')[-1].split('_')[:-1]) # "insilico_size100_1"
replicate = 1
### Importing data and calculating number of bins to use based on some rule (for binning methods)
input1_data_array = Building_MI_matrices.import_and_clean_input_file(file)
file = topology +"_" + expression_data_type + "_all.tsv"
bins_or_neighbors = Building_MI_matrices.calc_bins_num(input1_data_array.shape[1],binning_rule)

### open AUPR output file for summary statistics
partialname = filename.split('/')[-1][:-4] # insilico_size100_1_knockouts
output_filename = "AUPR_"+partialname+"_"+ expression_data_type + "_" + mi_est + "_" + infer_algo + ".dat"

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
TC_filename  = file[:-4] + "_TC_" + mi_est_string + ".mat"

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
print('Starting pos processing...')

gene_names = list(pd.read_csv(filename,sep='\t').columns)
outputdata           = pd.DataFrame(0,index = list(gene_names),columns = list(gene_names))
outputdata_threshold = pd.DataFrame(0,index = list(gene_names),columns = list(gene_names))

with open(partialname+'_all_MI2_FB10_MM.dat','r') as f:
        content = f.readlines()
content = [x.strip('\n').split(',')[2:] for x in content]
#print(content)
values_list = []
for item in content:
        outputdata.loc[item[0],item[1]] = float(item[2])
        values_list.append(float(item[2]))
values_list = np.array(values_list)
threshold = np.percentile(values_list, 75)

for item in content:
        if float(item[2])>=threshold:
                outputdata_threshold.loc[item[0],item[1]] = float(item[2])
        else:
                outputdata_threshold.loc[item[0],item[1]] = 0
outputdata.to_csv(output+'GRN_Iference_predicted_original.csv')
outputdata_threshold.to_csv(output+'GRN_Iference_predicted_threshold.csv')

test = os.listdir('.')

for item in test:
        if item.endswith(".dat") or item.endswith(".mat"):
                os.remove(item)
