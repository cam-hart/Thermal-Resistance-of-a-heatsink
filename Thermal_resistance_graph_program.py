# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 14:37:44 2021

@author: camer
"""
import matplotlib.pyplot as plt

N_fins=[]
R_sink=[]
with open("R_sink_f(N).txt", "r") as file:
    for N in 10,12,14,16,18:
        for line in file:
            a=line.split(' ')
            if int(a[0])==N:
                N_fins.append(int(a[0]))
                R_sink.append(float(a[1]))
                break
plt.plot(N_fins,R_sink)
plt.xlabel("Number of fins")
plt.ylabel("Thermal Resistance of heatsink (K/W)")
plt.show
            
            
        
        