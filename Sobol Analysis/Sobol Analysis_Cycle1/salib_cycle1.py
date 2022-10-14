# -*- coding: utf-8 -*-
"""
Created on Sun May 29 18:09:42 2022

@author: Leon Ritchie Salim
"""
import time

from datetime import datetime

import pandas as pd

from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import math
import odegfp_cycle1 as gfpmodel
import oderfp_cycle1 as rfpmodel
import timeit
import random
from multiprocessing import Process, Queue, current_process, freeze_support



PROC_NUM = 4
number_of_sample = 65536
name = 'odefull'
odeoutput = 'gfp'


def main(odeoutput1):
    
    name1 = name.lower()
    parameters, bounds= setups(name1)
    problem = {'num_vars': len(parameters),
           'names': parameters,
           'bounds': bounds}
    
    params = saltelli.sample(problem, number_of_sample ,calc_second_order=True)
    # quasi-random low-discrepancy sequencing
    # Number of params = N(2D +2), where D is the number of parameters
    # 3000*(2*22+2) = 138000 Parameters
    # Recomemded for N to be power of 2
    
    print(len(params))
    
    if False:
        params = log_distribute(params)

    done = Queue()
    paramsdiv=[]
    p = []
    remaining = 0
    for i in range(0, PROC_NUM):
        x = params.shape[0]/PROC_NUM
        remaining += x - int(x)
        x = int(x)
        y = params[x*i:x*(i+1), :]

        if i == PROC_NUM-1 and remaining != 0:
            y = np.append(y, params[-int(remaining):, :], axis=0)
        paramsdiv.append(y)
        p.append(Process(target=evaluate, args=(paramsdiv[i], done, i)))
        p[i].start()
    print(np.asarray(paramsdiv).shape)

    Ydiv = np.zeros(PROC_NUM).tolist()
    count = 0
    while True:
        temp, no =done.get()
        Ydiv[no] = temp.tolist()
        count += 1
        if count == PROC_NUM:
            break
    Y = []
    for i in range(0, PROC_NUM):
        Y += Ydiv[i]

    Y = np.asarray(Y)
    print(Y.shape)

    Si = sobol.analyze(problem, Y, print_to_console=True, calc_second_order=True)
    print(Si)
    write_file(Si, parameters, name)


def log_distribute(params):
    
    new = []
    for set in params:
        temp = []
        for i, x in enumerate(set):
            temp.append(math.log(x, 10))
        new.append(temp)
        
    mean = []
    std = []
    for i in range(len(params[0])):
        mean.append(np.mean(new[:][i]))
        std.append(np.std(new[:][i]))

    for i, x in enumerate(new):
        for j, y in enumerate(x):
            Z = random.gauss(0, 1)
            new[i][j] = 10**(mean[j] + std[j]*Z)
    return np.asarray(new)


def evaluate(values, done, no):
    
    Y = np.zeros([values.shape[0]])
    print("Start evaluation method with " + str(len(values)) + " parameter sets")
    start = timeit.default_timer()
    for i, X in enumerate(values):
        if odeoutput.lower() == 'gfp': 
            Y[i] = gfpmodel.main(X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7], X[8], X[9], X[10], X[11], X[12], X[13],\
                                      X[14], X[15], X[16], X[17], X[18],X[19])
        elif odeoutput.lower() == 'rfp':
            Y[i] = rfpmodel.main(X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7], X[8], X[9], X[10], X[11], X[12], X[13],\
                                      X[14], X[15], X[16], X[17], X[18],X[19])
        if i % 100 == 0 and i != 0:
            time_100 = timeit.default_timer() - start
            remaining_time=(int((len(values)-i)/100)*time_100)/60
            print('progress: ' + str((i/float(len(values)))*100) + ' %')
            print('time for 100 calculations: ' + str(time_100) + ' sec')
            print('remaining time: ' + str(remaining_time) + ' min\n')
            start = timeit.default_timer()
    done.put([Y,no])


def write_file(Si, parameters, name):
    
    # change the name of the file based on the boundary and number of samples and gfp or rfp model
    filename = 'salib_Cycle1_'+str(number_of_sample)+'_'+odeoutput+'_Boundary30Per_Time_1000.csv'     
   
    columns=("Param1,Param2,S1,S1con,S2,S2con,ST,STcon")
    file = open(filename, 'w')
    file.write(columns + '\n')
    S1=Si['S1']
    S1con=Si['S1_conf']
    ST=Si['ST']
    STcon=Si['ST_conf']
    S2=Si['S2'].tolist()
    S2con=Si['S2_conf'].tolist()
    for i,x in enumerate(S1):
        file.write(parameters[i] + ', - ,' + str(S1[i]) + ',' + str(S1con[i]) + ', - , - ,' + str(ST[i]) + ',' + str(STcon[i]) + '\n')
    for i in range(0,len(parameters)):
        for j in range(0,len(parameters)):
            if math.isnan(S2[i][j]) == False:
                file.write(parameters[i] + ',' + parameters[j] + ', - , - ,' + str(S2[i][j]) + ',' + str(S2con[i][j]) + ', - , - ' + '\n')
    file.close()



def setups(conf):
    
    data = pd.read_excel('Boundary Data1.xlsx',sheet_name = 'Boundary')
    parameters = data['Names'].tolist()
    lowbound = data['Lower Bound'].tolist()
    uppbound = data['Upper Bound'].tolist()
    b = []
    for i in range(len(lowbound)):
        x = []
        x.append(lowbound[i])
        x.append(uppbound[i])
        b.append(x)
    return parameters, b 



if __name__ == '__main__':
    main(name)


        

    
























