import sys
import os
import pandas as pd
import subprocess
import pathlib

#pathHere = str(pathlib.Path().resolve())
#os.system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'+pathHere)

argumentList = sys.argv[1:]
filename = ''
output = ''
for i in range(len(argumentList)):
        if argumentList[i] == '-i':
                filename = argumentList[i+1]
        elif argumentList[i] == '-o':
                output = argumentList[i+1]
if filename == '' or output == '':
        print('you MUST provide filename (-i <filename>) AND path out (-o <path_out>)')
        exit(1)
if not (filename[-3:] == 'csv' or filename[-3:] == 'tsv'):
        print('The input file MUST be a .csv OR .tsv file')
        exit(1)

#wasTSV = False
#if filename[-3:] == 'tsv':
#        df = pd.read_csv(filename,sep='\t')
#        df.to_csv(filename.replace('.tsv','.csv'),index=False)
#        filename = filename.replace('.tsv','.csv')
#        wasTSV = True

pathHere = str(pathlib.Path().resolve())
os.system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'+pathHere)

gene_names = list(pd.read_csv(filename).columns)

cmd = 'python3 pyPREMER.py '+ filename
#returned_value = subprocess.call(cmd, shell=True)  # returns the exit code in unix
os.system(cmd)

print('Starting pos processing...')
#if wasTSV:
#        os.remove(filename)

with open('./results/' + filename.split('/')[-1][:-4] + '/out_'+ filename.split('/')[-1][:-4]+'_prediction.txt') as f:
	content = f.readlines()
content = [x.strip('\n').split(',')[1:] for x in content[1:]]

outputdata      = pd.DataFrame(0,index = list(gene_names),columns = list(gene_names))
outputdata_threshold = pd.DataFrame(0,index = list(gene_names),columns = list(gene_names))

with open('./results/' + filename.split('/')[-1][:-4] + '/out_'+ filename.split('/')[-1][:-4]+'_threshold.txt') as f:
         threshold = f.readlines()
threshold = float(threshold[0].strip('\n'))

#print(threshold)

for item in content:
        #print('\n\n#######\n\nindex: ',item[0],'column: ',item[1],'value: ',item[2])
        outputdata.loc[item[0],item[1]] = float(item[2])
        if float(item[2])>=threshold:
                outputdata_threshold.loc[item[0],item[1]] = float(item[2])
        else:
                outputdata_threshold.loc[item[0],item[1]] = 0
if output[-1]!='/':
	output+='/'

outputdata.to_csv(output+'PREMER_predicted_original.csv')
outputdata_threshold.to_csv(output+'PREMER_predicted_threshold.csv')

os.system('rm -r '+ './results/' + filename.split('/')[-1][:-4]+'/')
