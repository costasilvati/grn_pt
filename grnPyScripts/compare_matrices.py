from multiprocessing import Process
import numpy as np
import pandas as pd
import sys
import os
from sklearn.metrics import confusion_matrix
#from sklearn.metrics import f1_score
#from sklearn.metrics import accuracy_score
#from sklearn.metrics import matthews_corrcoef

# Remove 1st argument from the
# list of command line arguments
argumentList = sys.argv[1:]
filenames = ''
output = ''
gold = ''
for i in range(len(argumentList)):
        if argumentList[i] == '-i':
                filenames = argumentList[i+1]
        elif argumentList[i] == '-o':
                output = argumentList[i+1]
                if output[-1] != '/':
                        output+='/'
        elif argumentList[i] == '-g':
                gold = argumentList[i+1]
if gold == '':
        print('MUST provide gold standard matrix')
        exit(1)

if filenames == '':
        print('No file to compare provided... Please provide one or more .csv file(s) (separeted by ",")')
        exit(1)

if output == '':
        print('No output path provided... writing in ./')

def write_Metrics(gold,filename,name,pathout):
        y_pred = pd.read_csv(filename,index_col=0).to_numpy().flatten()
        y_pred = np.where(y_pred>0,1,0)
        #print('GOLD: ',gold)
        #print('Y_PRED: ',y_pred)
        tn, fp, fn, tp = confusion_matrix(gold, y_pred).ravel()
        acc = (tp+tn)/(tp+tn+fp+fn)
        #acc = accuracy_score(y_true, y_pred)
        f1 = (2*tp)/(2*tp+fp+fn)
        #f1 = f1_score(y_true, y_pred)
        mcc = (tp*tn - fp*fn) / ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))**(1/2)
        #mcc = matthews_corrcoef(y_true, y_pred)
        precision = tp/(tp+fp)
        recall = tp /(tp+fn)

        with open(pathout+name+'.txt','w') as f:
                f.write('TP: '+ str(tp)+'\nTN: '+str(tn)+'\nFP: '+str(fp)+'\nFN: '+str(fn)+'\n')
                f.write('ACC: '+str(acc)+'\n')
                f.write('Precision: '+str(precision)+'\n')
                f.write('Recall: '+str(recall)+'\n')
                f.write('F1: '+str(f1)+'\n')
                f.write('MCC: '+str(mcc))

gold = pd.read_csv(gold,index_col=0)
gold = gold.to_numpy().flatten()
gold = np.where(gold>0,1,0)
jobs = []
for filename in filenames.split(','):
        name = filename.split('/')[-1].replace('.csv','')
        p = Process(target=write_Metrics, args=(gold,filename,name,output))
        jobs.append(p)
        p.start()

for proc in jobs:
        proc.join()
            
