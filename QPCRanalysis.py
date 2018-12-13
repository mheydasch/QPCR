#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 16:22:39 2018

@author: max
"""
#Use Following Naming convention
#Protein of Interest: Name= Proteinname, Type= Unknown
#Protein of Interest with GAPDH Primers: Name= Proteinname, Type= Calibrator (RQ)
#Ctrl DNA with primers for protein of Interest: Name= Proteinname, Type= Positive Control
#GAPDH for control DNA: Name= CTRLGAPDH, Type= Calibrator (RQ)
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
path='/Users/max/Desktop/Office/Phd/Data/N1E_115/SiRNA/SiRNA_20/SiRNA_20.xls'
data=pd.read_excel(path, skiprows=8, skip_footer=3, usecols=[1, 2 ,3])
expname=path.split('/')[-1].strip('.xls')
data=data.dropna()
averages=data.groupby(['Name', 'Type'], as_index=False)['Ct'].mean() #creates averages for the CT value of each sample
deltact=pd.DataFrame()
averages['delta Ct']=""
averages['expression level']="" #creates dataframe structure
#%%
d = averages['Name']=='CTRLGAPDH' #boolean mask for Calibrator of Ctrl
indexd= averages[d].index.values[0]
print(averages.groupby('Name').groups.keys())
for group in averages.groupby('Name'):
    a=group[1]['Type']=='Unknown'
    b=group[1]['Type']=='Calibrator (RQ)' 
    c=group[1]['Type']=='Positive Control' #boolean mask for Name groups where type is 'Positive Control'
    try:
        indexc=group[1][c]['Ct'].index.values[0] 
        indexb=group[1][b]['Ct'].index.values[0]
        indexa=group[1][a]['Ct'].index.values[0] #gets the index of value where boolean mask gives the 'Ct' value
        v=group[1][a]['Ct'].subtract(group[1][b]['Ct'][indexb]) #subtracts Calibrator from sample
        w=group[1][c]['Ct'].subtract(averages[d]['Ct'][indexd]) #subracts ctrl Calibrator from Ctrl
        delta=v-w[indexc] #subtracts the two above calculated differences, ctrl from test
        averages.loc[indexa, 'delta Ct']=delta[indexa] #adds the CT value of delta into the 'delta Ct' column and the row corresponding to the test sample
        deltapower=math.pow(2,-delta[indexa])
        averages.loc[indexa, 'expression level']=deltapower*100
    except IndexError:
        print(group[1]['Name'])
#%% figure making
plt.figure()        
toplot=averages[averages['expression level']!='']
bars=plt.bar(toplot['Name'], toplot['expression level'], width=0.3, color='r') #creating bars for expression levels
baseline=np.repeat(100, len(toplot))
#baseline=baseline.tolist()
n=0
for i in baseline: #subtracting expression level from baseline
    print(n)
    baseline[n]-=toplot['expression level'][toplot.index.values[n]]
    n+=1   
plt.bar(toplot['Name'], baseline, width=0.3, color='b', bottom=toplot['expression level']) #creating bars for baseline above expression levels
plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='on') #removing ticks
for spine in plt.gca().spines.values(): #removing grid
    spine.set_visible(False)    
for bar in bars: #adding expression levels on the bars
    plt.gca().text(bar.get_x() + bar.get_width()/2, bar.get_height() - 8, str(int(bar.get_height())) + '%', 
                 ha='center', color='black', fontsize=13)    
plt.title('Expression levels of knocked down RNA') #set title
fig1=plt.gcf()
plt.show()
fig1.savefig('{}{}_expression.jpg'.format(path.strip(('{}{}').format(expname,'xls')), expname),dpi=150) #save figure in folder of excel file
