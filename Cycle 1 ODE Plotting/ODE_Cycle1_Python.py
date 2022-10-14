# -*- coding: utf-8 -*-
"""
Created on Wed May  4 12:55:31 2022

@author: Leon Ritchie Salim
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


fig= plt.figure(figsize =(8,6) ,dpi = 300)

#### kn constants all in M/s
k1 = 0.00009627
k2 = 0.00042
k3 = 0.00001
k4 = 0.00002
k5 = 0.000001
k6 = 0.0008
alpha = 0.7 
k7 = 0.0001 

m1 = 0
m2 = 0.01
m3 = 0.01
m4 = 0.09

d1 = 0.00000802
d2 = 0.001
d3 = 0.001
d4 = 0.002
d5 = 0.0021
d6 = 0.0000688
d7 = 0.00000963
d8 = 0.0000642
d9 = 0.000138

#### Initial concentrations

c_h2o2 = 0
c_oxyrn = 0 # A number
c_oxyra = 0
c_pkat = 1
c_mtevp = 0
c_tevp = 0
c_gfpt = 0
c_gfpd = 0
c_rfpt = 0
c_rfpd = 0 


timerx = 1000
timeperiod = 1000

def wrap_full(k1,k2,k3,k4,k5,k6,alpha,k7,m1,m2,m3,m4,d1,d2,d3,d4,d5,d6,d7,d8,d9,
              c_dia,c_h2o2,c_oxyrn,c_oxyra,c_pkat,c_mtevp,c_tevp,c_gfpt,c_gfpd,c_rfpt,c_rfpd):
    
    global timerx, timeperiod
    
    def model(z, t):
        dia = z[0]
        h2o2 = z[1]
        oxyrn = z[2]
        oxyra = z[3]
        pkat = z[4]
        mtevp = z[5]
        tevp = z[6]
        gfpt = z[7]
        gfpd = z[8]
        rfpt = z[9]
        rfpd = z[10]
              

        ddiadt = m1 -k1*dia
        
        dh2o2dt = k1*dia -k2*oxyrn*h2o2 - d1*h2o2
        
        doxyrndt = m2 - k2*oxyrn*h2o2 - d2*oxyrn

        doxyradt = k2*oxyrn*h2o2 - d3*oxyra
        
        dpkatdt = -k3*pkat*oxyra 
        
        dmtevpdt = alpha*(oxyra/k7)/ (1+(oxyra/k7)) - d4*mtevp
        
        dtevpdt = k4*mtevp - d7*tevp
        
        dgfptdt = m3 - k5*gfpt*tevp - d5*gfpt
        
        dgfpddt = k5*gfpt*tevp - d8*gfpd
        
        drfptdt = m4 - k6*rfpt*tevp - d9*rfpt
        
        drfpddt = k6*rfpt*tevp - d6*rfpd
    
        
        dzdt = [ddiadt,dh2o2dt,doxyrndt,doxyradt,
                dpkatdt,dmtevpdt,dtevpdt,dgfptdt,
                dgfpddt,drfptdt,drfpddt]
        
        return dzdt

    z0 = [c_dia,c_h2o2,c_oxyrn,c_oxyra,c_pkat,c_mtevp,c_tevp,c_gfpt,c_gfpd,c_rfpt,c_rfpd]
    
    t = np.linspace(0, timerx, timeperiod)
    return odeint(model, z0, t)

def main_full(k1,k2,k3,k4,k5,k6,alpha,k7,m1,m2,m3,m4,d1,d2,d3,d4,d5,d6,d7,d8,d9,
              c_dia,c_h2o2,c_oxyrn,c_oxyra,c_pkat,c_mtevp,c_tevp,c_gfpt,c_gfpd,c_rfpt,c_rfpd,color1,label1): 
    
    global timerx, timeperiod
    z = wrap_full(k1,k2,k3,k4,k5,k6,alpha,k7,m1,m2,m3,m4,d1,d2,d3,d4,d5,d6,d7,d8,d9,
                  c_dia,c_h2o2,c_oxyrn,c_oxyra,c_pkat,c_mtevp,c_tevp,c_gfpt,c_gfpd,c_rfpt,c_rfpd)
    if True:
        t = np.linspace(0, timerx, timeperiod)
  
        plt.plot(t,z[:,7],str(color1)+'-',label=label1)
        plt.plot(t,z[:,10],str(color1)+'--',label=label1)

        
        plt.ylim(0,30)
        plt.xlim(0,timerx)
        plt.ylabel('Concentration')
        plt.xlabel('Time')
        plt.legend(loc='best')
    
        
        # print(z)
    # return (z[:,4][-1],z[:,5][-1])
    return z[:,7][-1], z[:,10][-1]


###################################

# if __name__ == '__main__':
    # c_dia = 0
#     r= main_full(k1,k2,k3,k4,k5,k6,alpha,k7,m1,m2,m3,m4,d1,d2,d3,d4,d5,d6,d7,d8,d9,
#                   c_dia,c_h2o2,c_oxyrn,c_oxyra,c_pkat,c_mtevp,c_tevp,c_gfpt,c_gfpd,c_rfpt,c_rfpd)
#     print(r)

###############################