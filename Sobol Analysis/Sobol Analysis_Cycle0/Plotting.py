# -*- coding: utf-8 -*-
"""
Created on Mon May 30 19:14:38 2022

@author: Leon Ritchie Salim
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

################################# Main Plotting ###############################################
pd.options.mode.chained_assignment = None  # default='warn'

filename = 'salib65536_gfp_Boundary30Per_Time_6000.csv'
code = filename[11:14].upper()
# code = filename[19:22].upper() # Only for Modifying Plot

data = pd.read_csv(filename)
print(len(data))
data1 = data[0:19]  #modify this value based on the number of constants you have

data1['S1'] = pd.to_numeric(data1["S1"], downcast="float")
data1['S1con'] = pd.to_numeric(data1["S1con"], downcast="float")
data1['ST'] = pd.to_numeric(data1["ST"], downcast="float")
data1['STcon'] = pd.to_numeric(data1["STcon"], downcast="float")

parameters = data1['Param1'].tolist()
s1data = data1['S1'].tolist()
s1devdata = data1['S1con'].tolist()
stdata = data1['ST'].tolist()
stdevdata = data1['STcon'].tolist()

# print(s1data)
plt.rcParams['figure.dpi'] = 500
para_axis = np.arange(len(parameters))
  
plt.bar(para_axis - 0.2, s1data, 0.4, label = '1st Order',yerr = s1devdata)
plt.bar(para_axis + 0.2, stdata, 0.4, label = 'Total',yerr = stdevdata)
  
plt.xticks(para_axis, parameters,fontsize = 5)
plt.xlabel("Parameters")
plt.ylabel("SI Value")
plt.legend(fontsize = 8)
plt.tight_layout()
plt.title(code+" 1st and Total SA",fontsize = 10)   
plt.xlim(left=-1)
# plt.ylim(0,0.7)
plt.ylim(0,0.18)
######################### Labelling Negative Value #########################
for i in range(len(s1data)):
    value = s1data[i]
    if value < 0:
        print(i+1)
        # print(i+1)
        plt.text(i-0.1, 0, '*', fontsize = 6)
plt.show()
################################################################################

data2 = data[22:len(data)]
Secord = data2[['Param1','Param2','S2','S2con']]
# print(Secord)