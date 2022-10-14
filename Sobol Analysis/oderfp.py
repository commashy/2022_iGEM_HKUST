# -*- coding: utf-8 -*-
"""
Created on Wed May  4 12:55:31 2022

@author: Leon Ritchie Salim
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#### Initial concentrations and constant
#### kn constants all in M/s


m1 = 0
c_dia = 0.01
c_h2o2 = 0
c_oxyrn = 0 # A number
c_oxyra = 0
c_poxy = 1
c_pkat = 1
c_mgfp = 0 
c_gfp = 0
c_mrfp = 0
c_rfp = 0

timerx = 6000
timeperiod = 6000

def wrap(k1,k2,k3,k4,k5,k6,alpha,beta,k,L,m2,d1,d2,d3,d4,d5,d6,d7,d8):
    
    global timerx, timeperiod
    
    def model(z, t):
        dia = z[0]
        h2o2 = z[1]
        oxyrn = z[2]
        oxyra = z[3]
        poxy = z[4]
        pkat = z[5]
        mgfp = z[6]
        gfp = z[7]
        mrfp = z[8]
        rfp = z[9]
              

        ddiadt = m1 -k1*dia - d1*dia
        
        dh2o2dt = k1*dia -k2*oxyrn*h2o2 - d2*h2o2
        
        doxyrndt = m2 - k2*oxyrn*h2o2 - d3*oxyrn

        doxyradt = k2*oxyrn*h2o2 - d4*oxyra
            
        dpoxydt = - k3*poxy*oxyra
        
        dpkatdt = -k4*pkat*oxyra 
        
        dmgfpdt = alpha*(oxyra/k)/ (1+(oxyra/k)) - d5*mgfp
        
        dgfpdt = k5*mgfp - d7*gfp
        
        dmrfpdt = beta*(oxyra/L)/ (1+(oxyra/L))  - d6*mrfp
        
        drfpdt = k6*mrfp - d8*rfp
        
        dzdt = [ddiadt,dh2o2dt,doxyrndt,doxyradt,
                dpoxydt,dpkatdt, dmgfpdt,dgfpdt, 
                dmrfpdt, drfpdt]
        
        return dzdt

    z0 = [c_dia,c_h2o2,c_oxyrn, c_oxyra,c_poxy,c_pkat,c_gfp,c_mgfp,c_mrfp,c_rfp]
    
    t = np.linspace(0, timerx, timeperiod)
    return odeint(model, z0, t)

def main(k1,k2,k3,k4,k5,k6,alpha,beta,k,L,m2,d1,d2,d3,d4,d5,d6,d7,d8): 
    
    global timerx, timeperiod
    z = wrap(k1,k2,k3,k4,k5,k6,alpha,beta,k,L,m2,d1,d2,d3,d4,d5,d6,d7,d8)
    if True:
 
    # return (z[:,4][-1],z[:,5][-1])
        return z[:,9][-1]


###################################
if __name__ == '__main__':
    c_dia = 0.01
    rfp1= main(k1,k2,k3,k4,k5,k6,alpha,beta,k,L,m2,d1,d2,d3,d4,d5,d6,d7,d8)
    print(rfp1)



###############################
