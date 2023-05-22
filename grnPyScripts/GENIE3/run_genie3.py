import numpy as np
from GENIE3 import *
import pandas as pd
import sys
import os

# Remove 1st argument from the
# list of command line arguments
argumentList = sys.argv[1:]
filename = ''
output = ''
for i in range(len(argumentList)):
	if argumentList[i] == '-i':
		filename = argumentList[i+1]
	elif argumentList[i] == '-o':
		outputpath = argumentList[i+1]
		if outputpath[-1] != '/':
			outputpath+='/'

if filename == '' or outputpath == '':
	print('you MUST provide filename (-i <filename>) AND path out (-o <path_out>)')
	exit(1)

if filename[-3:] == 'csv':
	print('Lendo arquivo ',filename)
	df = pd.read_csv(filename)
elif filename[-3:] == 'tsv':
	df = pd.read_csv(filename,sep='\t')
else:
	print('The input file MUST be a .csv OR .tsv')
	exit(1)

data = df.values
outputdata           = pd.DataFrame(0,index = list(df.columns),columns = list(df.columns))
outputdata_threshold = pd.DataFrame(0,index = list(df.columns),columns = list(df.columns))

gene_names = list(df.columns)

VIM = GENIE3(data)

get_link_list(VIM,gene_names=gene_names,file_name='GENIE3_ranking.txt')

# Pos processamento para deixar no formato de matrix
print('Starting pos processing...')
with open('GENIE3_ranking.txt','r') as f:
        content = f.readlines()
content = [x.strip('\n').split('\t') for x in content]

values_list = []
for item in content:
        outputdata.loc[item[0],item[1]] = float(item[2])
        values_list.append(float(item[2]))

threshold = np.percentile(np.array(values_list), 75)
for item in content:
        if float(item[2])>=threshold:
                outputdata_threshold.loc[item[0],item[1]] = float(item[2])
        else:
                outputdata_threshold.loc[item[0],item[1]] = 0

outputdata.to_csv(outputpath+'GENIE3_predicted_original.csv')
outputdata_threshold.to_csv(outputpath+'GENIE3_predicted_threshold.csv')

os.remove('GENIE3_ranking.txt')
print('Feito')
