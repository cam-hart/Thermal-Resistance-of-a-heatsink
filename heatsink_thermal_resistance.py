# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 15:19:54 2021

@author: camer
"""
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
"""
num_fins=N(10,12,14,16,18)
base_thickness=x(10mm)
base_width=w (46.9mm)
gap_thickness=b
fin_thickness=t(1.5mm)
fin_height=H(30mm)
fin_length=L (60mm)
metal_thermal_conductivity=Km (180 W/mk)
fluid_thermal_conductivity=Kf (2.624kW/m*K)
dynamic_viscocity=Vd (1.846*(10**-5))
kinematic_viscosity=Vk (1.586*(10**-5))
freestream_velocity=v_f
Specific_heat_capacity=Cp (1.0049kJ/kg*K)
Q=air flow(m^3/min)
fluid density (1.777kg/m^3)
"""
fan_curve=pd.read_excel('Copy of FanCurveData.xlsx',index_col=None,usecols='A,B')
air_flow=np.array(fan_curve['Air Flow [m^3/minute]'].tolist())
pressure=np.array(fan_curve['Static Pressure [mm-H2O]'].tolist())   
#importing excel data, making arrays for graph
plt.scatter(pressure,air_flow)
def volume_flow_fP(x,a,b,c,d,e,f):
    return (a*(x**5))+(b*(x**4))+(c*(x**3))+(d*(x**2))+(e*x)+f
#making a linear function for volume flow rate as a function of pressure, originally
#did the other way around but realised needs to be flipped
coefficients,other=curve_fit(volume_flow_fP,pressure,air_flow)
a,b,c,d,e,f=coefficients
x=np.arange(pressure[len(pressure)-1],pressure[0],0.001)
plt.plot(x,volume_flow_fP(x,a,b,c,d,e,f))
plt.ylabel("Air flow rate (m^3/s)")
plt.xlabel("Pressure drop (mm(H2O)")
plt.show()
#plotting function to see how accurate it is. Tried function in form a**b, but didn't work. Left order of
#function at 5 because it fits pretty well

class Rsink:
    def __init__(self,N,Q):
        self.N=N
        self.Q=Q
        self.x=0.01
        self.w=0.047
        self.t=0.0015
        self.H=0.03
        self.L=0.06
        self.Km=180
        self.Kf=2624
        self.density=1.777
        self.Vd=1.846*(10**-5) 
        self.Vk=1.586*(10**-5)
        self.Cp=1.0049
        self.b=(self.w-(self.N*self.t))/(self.N-1)
        self.v_f=(self.Q*60)/0.041*0.012
        self.mu=1-((self.N*self.t)/self.w)
        self.Kc=0.42*(1-(self.mu**2))
        self.Ke=(1-(self.mu**2))**2
        self.v_ch=self.v_f*(1+(self.t/self.b))
        self.Rech=(self.v_ch*self.L)/self.Vk
        self.l=self.L/(2*self.b*self.Rech)
        self.prandtl_number=(self.Cp*self.Vd)/self.Kf
        self.Reb=(self.b*self.v_ch)/self.Vk
        self.nusselt_number=(((self.Reb*self.prandtl_number*0.5)**-3)+((0.664*math.sqrt(self.Reb*(self.prandtl_number**1/3)*(1+(3.65/math.sqrt(self.Reb)))))**-3))**(-1/3)
        self.h=(self.Kf*self.nusselt_number)/self.b
        self.fRech=24-(32.527*(self.b/self.H))+(46.721*((self.b/self.h)**2))
        self.Fapp=((((3.44/math.sqrt(self.l))**2)+(self.fRech**2))**0.5)/self.Rech
        self.P_drop_pa=(((self.Fapp*self.N*((2*self.H*self.L)+(self.b*self.L)))/(self.H*self.w))+self.Kc+self.Ke)*0.5*self.density*(self.v_ch**2)
        self.P_drop_mmh20=self.P_drop_pa/9.80665
        self.P=2*(self.t+self.L)
        self.Ac=self.t*self.L
        self.m=math.sqrt((self.h*self.P)/(self.Km*self.Ac))
        self.R_fin=1/math.sqrt(self.h*self.P*self.Km*self.Ac)*math.tanh(self.m*self.H)
        self.R_sink=1/(self.N/self.R_fin)+(self.h*self.N*self.b*self.L)
        #making a class to contain all Rsink data. Starting with the variables, the constants @300K,
        #then coding the equations, ensuring theyre in the right order.
        self.Q2=volume_flow_fP(self.P_drop_mmh20, a, b, c, d, e, f)

Q=0.1
N=int(input("Enter number of fins:"))
Previous=0
while True:
    Q=Rsink(N,Q).Q2
    Q=Rsink(N,Q).Q2
    if Q==Previous:
        break
    Previous=Q
#finding operating point using a while loop until volume flow rate converges, ie equal to last
#for N=10 and N=18, it oscillates between two very close points (less than 1*10^-15), so changed it to return Q when 
#its equal to the second previous, by repeating line

print("R-Sink for "+str(N)+" fins is: "+str(Rsink(N,Q).R_sink)+" (K/W)")
#printing R-sink using operating point
    
with open("R_sink_f(N).txt", "a") as file:
    file.write(str(N)+" "+str(Rsink(N,Q).R_sink)+'\n')
#adds result to file
   
    
    



    




