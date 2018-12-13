#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 11:57:12 2018

@author: max
"""

import pandas as pd
import numpy as np
import itertools
Experimentname='SiRNA_20'
filepath='/Users/max/Desktop/Office/Phd/Data/N1E_115/SiRNA/SiRNA_20/'
samples=['Rac', 'ARHGAP17', 'DOCK9']
samplelist=np.repeat((samples),6)
Ctrllist=np.repeat('Ctrl',len(samplelist)/2+3)
Waterlist=np.repeat('H2O', len(samplelist)/6+1)
GAPDH1=np.repeat('GAPDH', 1)
silist = np.hstack((samplelist, Ctrllist, Waterlist))
primerlist=[]
for i in samples:
    primers=np.repeat(i, 3)
    GAPDH=np.repeat('GAPDH', 3)
    primerlist=np.concatenate((primerlist, primers, GAPDH))
primerlist=np.concatenate((primerlist, np.repeat(samples, 3), np.repeat('GAPDH', 3), samples, GAPDH1))    

df=pd.DataFrame({'Knockdown': silist, 'Primer': primerlist},)
Mastermix=6.25*(len(primerlist)+2)
Water=2.75*(len(primerlist)+2)
primermix=sum(primerlist==samples[1])*9
primer=primermix/18

GAPDHmix=(sum(primerlist=='GAPDH'))*9
GAPDHprimer=GAPDHmix/18
mixindices=['Total Mastermix master', 'Total Mastermix water','Mastermix Sample','2x primer','GAPDH mastermix','2x GAPDH primer']
mixdict=pd.DataFrame({'Knockdown' : ['Total Mastermix master', 'Total Mastermix water', 'Mastermix Sample', '2x primer', 'GAPDH mastermix', '2x GAPDH primer'], 
                     'Primer' : [Mastermix, Water, primermix, primer, GAPDHmix, GAPDHprimer]})
df=df.append(mixdict, ignore_index=True)
df.index +=1
df.to_csv('{}{}QPCR_planning.csv'.format(filepath, Experimentname)) 

       

#Sampleprimermix= 
#df=pd.DataFrame({'SiRNA':np.repeat(['H2O', 'CTRL'],len(samplelist)), 'Primer':np.repeat(['H20', 'CTRl'],len(samplelist))})

#df=pd.DataFrame({'well':[], 'n':[],'mean':[], 'STDV':[], 'normal':[],'Test':[], 'F_value':[], 'P_value':[]},)
