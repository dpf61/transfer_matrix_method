# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 18:35:19 2023

@author: dpf61
"""
import numpy as np
import matplotlib.pyplot as plt
from typing import Final
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import CubicSpline

plt.style.use(['science', 'notebook', 'grid'])

import colour

from colour import XYZ_to_sRGB, XYZ_to_Lab
from colour import xy_to_XYZ, XYZ_to_xyY, xyY_to_XYZ, XYZ_to_RGB
import TransferMatrixMethod as tf

###############################################################################
data = np.loadtxt('SiO2.csv', skiprows=1, delimiter=';').transpose()
landas = data[0]
n = data[1] 
k = data[2]
nombre = "$SiO_2$"

vidrio = tf.Material(n, k, landas, nombre)

aire=tf.Material(np.ones(2000), np.zeros(2000), np.arange(0,2000), "aire")

data = np.loadtxt('Crystalline_Sb2S3_Bari.csv', skiprows=1, delimiter=';').transpose()
landa = data[1]
n = data[4]
k = data[5]

cristalino = tf.Material(n, k, landa, "$Sb_2S_3$ crist")

data = np.loadtxt('Amorphous_Sb2S3_Bari.csv', skiprows=1, delimiter=';').transpose()
landa = data[1]
n = data[4]
k = data[5]

amorfo = tf.Material(n, k, landa, "$Sb_2S_3$ amorfo")

###############################################################################



landa = 633
theta = np.arange(0,90.1,0.1) #grados
pol = 'p'


d_amorfo = [290*1e-9]

d_cristalino = [255*1e-9]

R_p_cristalino = [[] for j in range(np.size(d_cristalino))]

R_p_amorfo = [[] for j in range(np.size(d_amorfo))]

R_s_cristalino = [[] for j in range(np.size(d_cristalino))]

R_s_amorfo = [[] for j in range(np.size(d_amorfo))]


for k in range(0,np.size(d_amorfo)):
    
    R = [] 
    
    for i in range(0,np.size(theta)):
        
        stack = tf.Multicapa()
        
        stack.añade_capa(amorfo.get_n(landa), amorfo.get_k(landa), d_amorfo[k], amorfo.get_nombre)

        stack.añade_capa(vidrio.get_n(landa), vidrio.get_k(landa), d_amorfo[k], vidrio.get_nombre)
        
        R.append(stack.reflectancia(np.radians(theta[i]), landa*1e-9, pol))

    # guardo el array de R's en el pack
    R_p_amorfo[k] = R


for k in range(0,np.size(d_cristalino)):
    
    R = [] 
    
    for i in range(0,np.size(theta)):
        
        stack = tf.Multicapa()
        
        stack.añade_capa(cristalino.get_n(landa), cristalino.get_k(landa), d_cristalino[k], cristalino.get_nombre)

        stack.añade_capa(vidrio.get_n(landa), vidrio.get_k(landa), d_cristalino[k], vidrio.get_nombre)
        
        R.append(stack.reflectancia(np.radians(theta[i]), landa*1e-9, pol))

    # guardo el array de R's en el pack
    R_p_cristalino[k] = R

for k in range(0,np.size(d_amorfo)):
    
    R = [] 
    
    for i in range(0,np.size(theta)):
        
        stack = tf.Multicapa()
        
        stack.añade_capa(amorfo.get_n(landa), amorfo.get_k(landa), d_amorfo[k], amorfo.get_nombre)

        stack.añade_capa(vidrio.get_n(landa), vidrio.get_k(landa), d_amorfo[k], vidrio.get_nombre)
        
        R.append(stack.reflectancia(np.radians(theta[i]), landa*1e-9, 's'))

    # guardo el array de R's en el pack
    R_s_amorfo[k] = R


for k in range(0,np.size(d_cristalino)):
    
    R = [] 
    
    for i in range(0,np.size(theta)):
        
        stack = tf.Multicapa()
        
        stack.añade_capa(cristalino.get_n(landa), cristalino.get_k(landa), d_cristalino[k], cristalino.get_nombre)

        stack.añade_capa(vidrio.get_n(landa), vidrio.get_k(landa), d_cristalino[k], vidrio.get_nombre)
        
        R.append(stack.reflectancia(np.radians(theta[i]), landa*1e-9, 's'))

    # guardo el array de R's en el pack
    R_s_cristalino[k] = R


### DATOS EXPERIMENTALES#######################################################
Imax = 1000

data = np.genfromtxt('amorfo_rojo.txt', skip_header=1, delimiter='\t').transpose()
theta_1 = data[0]
power_1 = data[1]
theta_2 = data[2]
power_2 = data[3]
theta_3 = data[4]
power_3 = data[5]

#PROMEDIO AMORFO
power_1_f = interp1d(theta_1, power_1)
power_2_f = interp1d(theta_2, power_2)
power_3_f = interp1d(theta_3, power_3)

th = np.linspace(5,80,1000)

promedio_amorfo = []
std_amorfo = []
    
# Iterar sobre los elementos correspondientes de los arrays
for valor1, valor2, valor3 in zip(power_1_f(th), power_2_f(th), power_3_f(th)):
    valor_medio = np.mean([valor1,valor2,valor3])
    std = np.std([valor1,valor2,valor3])
    promedio_amorfo.append(valor_medio)
    std_amorfo.append(std)



data = np.genfromtxt('cristalino_rojo.txt', skip_header=1, delimiter='\t').transpose()
thetac_1 = data[0]
powerc_1 = data[1]
thetac_2 = data[2]
powerc_2 = data[3]
thetac_3 = data[4]
powerc_3 = data[5]

#PROMEDIO CRISTALINO
powerc_1_f = interp1d(thetac_1, powerc_1)
powerc_2_f = interp1d(thetac_2, powerc_2)
powerc_3_f = interp1d(thetac_3, powerc_3)

th = np.linspace(5,80,1000)

promedio_cristalino = []
std_cristalino = []

# Iterar sobre los elementos correspondientes de los arrays
for valor1, valor2, valor3 in zip(powerc_1_f(th), powerc_2_f(th), powerc_3_f(th)):
    valor_medio = np.mean([valor1,valor2,valor3])
    std = np.std([valor1,valor2,valor3])
    promedio_cristalino.append(valor_medio)
    std_cristalino.append(std)



###############################################################################



R_p_amorfo = np.array(R_p_amorfo[0])
R_p_cristalino = np.array(R_p_cristalino[0])
R_s_cristalino = np.array(R_s_cristalino[0])
R_s_amorfo = np.array(R_s_amorfo[0])


s = [0.0,0.1,0.2,0.3,0.4,0.5]
p = [1.0,0.9,0.8,0.7,0.6,0.5]



#%% GRAFICA ADICIONALES
plt.figure()
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('$P$ / $\mu$W', fontsize=18)
plt.title(r'$Sb_2S_3$ Amorfo ($\lambda=633$ nm)', fontsize=18)
plt.plot(theta_1, power_1, 'o-', label='Medida 1')
plt.plot(theta_2, power_2, 'o-', label='Medida  2')
plt.plot(theta_3, power_3,'o-', label='Medida  3')
plt.grid(None)
plt.legend(fontsize=8)
plt.xlim([0,90])
# plt.ylim([0,1])

plt.figure()
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('$P$ / $\mu$W', fontsize=18)
plt.title(r'$Sb_2S_3$ Cristalino ($\lambda=633$ nm)', fontsize=18)
plt.plot(thetac_1, powerc_1,'o-', label='Medida  1')
plt.plot(thetac_2, powerc_2,'o-', label='Medida  2')
plt.plot(thetac_3, powerc_3,'o-', label='Medida  3')
plt.grid(None)
plt.legend(fontsize=8)
plt.xlim([0,90])
# plt.ylim([0,1])



