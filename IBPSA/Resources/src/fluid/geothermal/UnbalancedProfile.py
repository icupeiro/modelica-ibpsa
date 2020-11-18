# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:14:43 2018

@author: u0112721
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D
from cycler import cycler
from loadData import *
import os

#plt.style.use('grayscale')
plt.rcParams['axes.prop_cycle'] =  cycler(color= ['#1f77b4', '#2ca02c', '#9467bd', '#d62728', '#ff7f0e'])
plt.rcParams["figure.facecolor"] = 'w'
#plt.rc('xtick', labelsize='x-small')
#plt.rc('ytick', labelsize='x-small')
plt.rcParams['lines.linewidth'] = 1
#plt.rcParams['legend.fontsize'] = 'x-small'
#plt.rcParams['axes.labelsize'] = 12
plt.rc('figure', figsize=(20.69,11.69))
#plt.rc('xtick', labelsize=32)
plt.rc('font', family='serif')

os.chdir("/home/u0112721/Git/borefieldMPC/borefieldMPC/")

reload = True

path = []
path.append("UnbalanceProfile.mat")

stopTime = 31536000*5

varDic = {}

varDic['load1'] = 'product1.y'
varDic['load2'] = 'product2.y'

varDic['T1'] = 'twoUTube.borHol.intHex[1].vol2.T'
varDic['T2'] = 'twoUTube1.borHol.intHex[1].vol2.T'

############# Load data #######################################################
df={}
tr=[]
for i in range(len(path)):
    if reload:
        data = LoadData(path=path[i],varDic=varDic, stop_time=stopTime,year=2018, start_time=dt.datetime(2018,1,8))
        df[path[i]] = data.df 
        df[path[i]].to_pickle(path[i].replace('.mat',''))
    else:
        df[path[i]] = pd.read_pickle(path[i].replace('.mat',''))

print('Data loaded!')

######## Set index

for i in range(len(path)):
     df[path[i]]['Year'] =  (df[path[i]]['time'] - df[path[i]]['time'][0])/31536000
     df[path[i]] = df[path[i]].set_index('Year')




(df[path[0]]['T1']-273.15).plot(ylim=[-3,15], title='Borehole outlet temperature')
(df[path[0]]['T2']-273.15).plot(ylim=[-3,15], linestyle = '--')
plt.margins(0,0)
plt.savefig(fname, bbox_inches='tight', pad_inches=0)


plt.show()

(df[path[0]]['load1']-273.15).plot(ylim=[-3,12], title='Borehole outlet temperature')
(df[path[0]]['load2']-273.15).plot(ylim=[-3,12], linestyle='dashed')
plt.margins(0,0)
plt.savefig(fname, bbox_inches='tight', pad_inches=0)
plt.show()



