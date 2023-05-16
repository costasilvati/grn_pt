import pandas as pd
import os

with open('ranking.txt','r') as f:
	content = f.readlines()
content = [x.strip('\n').split('\t') for x in content]

#df = pd.read_csv('../PREMER/insilico_size100_1_knockouts.csv')

f = open('data.txt')
gene_names = f.readline()
f.close()
gene_names = gene_names.rstrip('\n').split('\t')

outputdata = pd.DataFrame(0,index = list(gene_names),columns = list(gene_names))

for item in content:
	print('index: ',item[0],'columns: ',item[1],'value: ',item[2])
	#outputdata.loc[item[0]][item[1]] = float(item[2])
	outputdata.loc[item[0],item[1]] = float(item[2])
	print(outputdata)
outputdata.to_csv('GENIE3_predicted_teste.csv')
