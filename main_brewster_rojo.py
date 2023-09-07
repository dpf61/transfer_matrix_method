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
power_1 = data[1]/Imax
theta_2 = data[2]
power_2 = data[3]/Imax
theta_3 = data[4]
power_3 = data[5]/Imax

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
powerc_1 = data[1]/Imax
thetac_2 = data[2]
powerc_2 = data[3]/Imax
thetac_3 = data[4]
powerc_3 = data[5]/Imax

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



#%% GRAFICA EXPERIMENTAL AMORFO
plt.figure()
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia ($\lambda=633$ nm)', fontsize=18)
plt.title(r'$Sb_2S_3$ Amorfo', fontsize=18)
plt.plot(theta_1, power_1, 'o-', label='Medida 1')
plt.plot(theta_2, power_2, 'o-', label='Medida  2')
plt.plot(theta_3, power_3,'o-', label='Medida  3')
plt.grid(None)
plt.legend(fontsize=8)
plt.xlim([0,90])
# plt.ylim([0,1])


## teorica amorfo
plt.figure()
plt.title(f'$Sb_2S_3$ Amorfo ($d={d_amorfo[0]*1e9:.0f}$ nm)', fontsize=18)
for k in range(np.size(s)):  
    plt.plot(theta, p[k]*R_p_amorfo+s[k]*R_s_amorfo, label=f'{p[k]*100:.0f} %')
plt.grid(None)
plt.legend(title='Polarización $p$', fontsize=10)
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia ($\lambda=633$ nm)', fontsize=18)
plt.xlim([0,90])
plt.ylim([0,1])


#%% GRAFICA EXPERIMENTAL CRISTRALINO
plt.figure()
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia ($\lambda=633$ nm)', fontsize=18)
plt.title(r'$Sb_2S_3$ Cristalino', fontsize=18)
plt.plot(thetac_1, powerc_1,'o-', label='Medida  1')
plt.plot(thetac_2, powerc_2,'o-', label='Medida  2')
plt.plot(thetac_3, powerc_3,'o-', label='Medida  3')
plt.grid(None)
plt.legend(fontsize=8)
plt.xlim([0,90])
# plt.ylim([0,1])

## teorica cristalino
plt.figure()
plt.title(f'$Sb_2S_3$ Cristalino ($d={d_cristalino[0]*1e9:.0f}$ nm)', fontsize=18)
for k in range(np.size(s)):  
    plt.plot(theta, p[k]*R_p_cristalino+s[k]*R_s_cristalino, label=f'{p[k]*100:.0f} %')
plt.grid(None)
plt.legend(title='Polarización $p$', fontsize=10)
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia ($\lambda=633$ nm)', fontsize=18)
plt.xlim([0,90])
plt.ylim([0,1])
    



#%% COMPARATIVA TEÓRICA VS POLARIZACIÓN

# AMORFO
plt.figure()
plt.plot(theta, 0.9*R_p_amorfo+0.1*R_s_amorfo, label = 'Simulación (90% polarización $p$)')
plt.errorbar(th, promedio_amorfo, yerr=std_amorfo, capsize=0, 
              fmt='none', 
              # label=r'$\langle R_{exp} \rangle \pm \sigma_{\langle R_{exp} \rangle}$ Valor medio exp.',
              color='g',
              alpha=0.1)
plt.plot(th, promedio_amorfo, 'g')
prom_a_f = interp1d(th, promedio_amorfo)
plt.plot(theta_3, prom_a_f(theta_3),'og', label='Experimental')
plt.grid(None)
plt.xlim([0,90])
plt.ylim([0,1])
plt.title(f'$Sb_2S_3$ Amorfo ($d={d_amorfo[0]*1e9:.0f}$ nm)', fontsize=18)
plt.legend(fontsize=14)
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia ($\lambda=633$ nm)', fontsize=18)


# CRISTALINO
plt.figure()
plt.plot(theta, 1.0*R_p_cristalino + 0.0*R_s_cristalino, label = 'Simulación (100% polarización $p$)')
# xerr = np.ones_like(std_cristalino)*2
plt.errorbar(th, promedio_cristalino, yerr=std_cristalino, capsize=0, 
              fmt='none', 
              # label=r'$\langle R_{exp} \rangle \pm \sigma_{\langle R_{exp} \rangle}$ Valor medio exp.',
              color='g',
              alpha=0.1)
plt.plot(th, promedio_cristalino, 'g')
prom_c_f = interp1d(th, promedio_cristalino)
plt.plot(thetac_3, prom_c_f(thetac_3),'og', label='Experimental')
plt.grid(None)
plt.xlim([0,90])
plt.ylim([0,1])
plt.title(f'$Sb_2S_3$ Cristalino ($d={d_cristalino[0]*1e9:.0f}$ nm)', fontsize=18)
plt.legend(fontsize=14)
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia ($\lambda=633$ nm)', fontsize=18)


#%% GRAFICAS TEORICAS ESPESORES

color_crist = 'brown'

color_amorfo = 'darkblue'


landa = 633
theta = np.arange(0,90.1,1) #grados
pol = 'p'


d_amorfo = np.linspace(280,300,100)*1e-9

d_cristalino = np.linspace(245,265,100)*1e-9

R_p_cristalino = [[] for j in range(np.size(d_cristalino))]

R_p_amorfo = [[] for j in range(np.size(d_amorfo))]

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
        
        R.append(stack.reflectancia(np.radians(theta[i]), landa*1e-9, pol))

    # guardo el array de R's en el pack
    R_p_cristalino[k] = R


plt.figure()

for k in range(np.size(d_amorfo)):

    # DESOMENTAR SI SE QUIEREN LINEAS CON LABEL DE CADA ESPESOR
    # if k==0:
    #     plt.plot(theta, R_p_amorfo[k], color='red', label=f'($d={d_amorfo[k]*1e9:.0f}$ nm)')
    
    # elif k==np.size(d_amorfo)/2:
    #     plt.plot(theta, R_p_amorfo[k], color=color_amorfo, label=f'($d={d_amorfo[k]*1e9:.0f}$ nm)')
        
    # elif k == np.size(d_amorfo)-1:
    #     plt.plot(theta, R_p_amorfo[k], color='yellow', label=f'($d={d_amorfo[k]*1e9:.0f}$ nm)')
    
    # else:      
    #     plt.plot(theta, R_p_amorfo[k], color=color_amorfo, alpha=0.1)
        
    plt.plot(theta, R_p_amorfo[k], color=color_amorfo, alpha=0.1)

plt.plot(theta, R_p_amorfo[int(np.size(d_amorfo)/2)], color=color_amorfo, label='Simulación')
    
plt.grid(None)
plt.xlim([0,81])
plt.ylim([0,0.3])
plt.title(f'$Sb_2S_3$ Amorfo ($\lambda=633$ nm)', fontsize=18)
plt.legend()
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia', fontsize=18)
    


plt.figure()

for k in range(np.size(d_cristalino)):
       
    plt.plot(theta, R_p_cristalino[k], color=color_crist, alpha=0.1)

plt.plot(theta, R_p_cristalino[int(np.size(d_cristalino)/2)], color = color_crist, label='Simulación')
    
plt.grid(None)
plt.xlim([0,81])
plt.ylim([0,0.3])
plt.title(f'$Sb_2S_3$ Cristalino ($\lambda=633$ nm)', fontsize=18)
plt.legend()
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia', fontsize=18)
    
    
    
    
# GRAFICAS EXPERIMENTALES VALOR MEDIO 
    


# AMORFO
plt.figure()
plt.errorbar(th, promedio_amorfo, yerr=std_amorfo, capsize=0, 
              fmt='none', 
              # label=r'$\langle R_{exp} \rangle \pm \sigma_{\langle R_{exp} \rangle}$ Valor medio exp.',
              color=color_amorfo,
              alpha=0.1)
# plt.plot(th, promedio_amorfo, 'g')
prom_a_f = interp1d(th, promedio_amorfo)
plt.plot(theta_3, prom_a_f(theta_3),'-o',color=color_amorfo, label='Experimental')
plt.grid(None)
plt.xlim([0,81])
plt.ylim([0,0.3])
plt.title(f'$Sb_2S_3$ Amorfo ($\lambda=633$ nm)', fontsize=18)
plt.legend()
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia', fontsize=18)


# CRISTALINO
plt.figure()
plt.errorbar(th, promedio_cristalino, yerr=std_cristalino, capsize=0, 
              fmt='none', 
              # label=r'$\langle R_{exp} \rangle \pm \sigma_{\langle R_{exp} \rangle}$ Valor medio exp.',
              color=color_crist,
              alpha=0.1)
# plt.plot(th, promedio_cristalino, 'g')
prom_c_f = interp1d(th, promedio_cristalino)
plt.plot(thetac_3, prom_c_f(thetac_3),'-o',color=color_crist, label='Experimental')
plt.grid(None)
plt.xlim([0,81])
plt.ylim([0,0.3])
plt.title(f'$Sb_2S_3$ Cristalino ($\lambda=633$ nm)', fontsize=18)
plt.legend()
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia', fontsize=18)

# grafica contraste teorica
plt.figure()

R_log = [[] for j in range(np.size(d_amorfo))]
for k in range(np.size(d_amorfo)):
    R_log[k] = 10*np.log10(np.array(R_p_cristalino[k]) / np.array(R_p_amorfo[k]))
    

for C in R_log:
    plt.plot(theta, C, 'red', alpha=0.05)

plt.plot(theta, R_log[int(np.size(d_amorfo)/2)], 'red', label='Simulación')

# GRAFICA DE CONTRASTE 

C = 10*np.log10(np.array(prom_c_f(thetac_3))/np.array(prom_a_f(thetac_3)))


plt.plot(thetac_3, C, '-ok', label='Experimental')
plt.title(f'Contraste $Sb_2S_3$ ($\lambda=633$ nm)', fontsize=18)
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('$C$ / dB', fontsize=18)
plt.ylim([-40,40])
plt.xlim([0,81])
plt.grid(None)
plt.legend()
plt.plot(theta, np.zeros_like(theta), '--k', linewidth=0.5)


#%% BREWSTERS


plt.figure(figsize=(10,5))

sim = [63.0,64.0]
sim_y = ['Amorfo', 'Cristalino']
err_sim = [2.0,2.0]

exp = [68.0]
exp_y = 'Cristalino'
err_exp = [1.0]

snell = [65.0,70.0]

plt.plot(sim,sim_y,'o')

plt.errorbar(sim, sim_y, xerr=err_sim, capsize=5, 
              fmt='none', 
              label= 'Simulación',
              color='blue')

plt.plot(exp, exp_y, 'og')
plt.errorbar(exp, exp_y, xerr=err_exp, capsize=5, 
              fmt='none', 
              label= 'Experimental',
              color='green')

plt.plot(snell,sim_y, 'or' ,label='Snell')

plt.legend(loc='center right')

plt.title('Ángulo de Brewster', fontsize=18)
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.legend(fontsize=14, loc='center right')

plt.grid(None)

#%% GRAFICA POLARIZACION CON ESPESORES

landa = 633
theta = np.arange(0,90.1) #grados
pol = 'p'


d_amorfo = np.linspace(280,300,100)*1e-9

d_cristalino = np.linspace(245,265,100)*1e-9

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

###### COMPARATIVA TEÓRICA VS POLARIZACIÓN

# AMORFO
plt.figure()

for k in range(np.size(d_amorfo)):
    if k == int(np.size(d_cristalino)/2):
        plt.plot(theta, 0.9*np.array(R_p_amorfo[k])+0.1*np.array(R_s_amorfo[k]),
             'darkblue',
             label = 'Simulación (90% polarización $p$)')
    else:
        plt.plot(theta, 0.9*np.array(R_p_amorfo[k])+0.1*np.array(R_s_amorfo[k]),'darkblue', alpha=0.01)       

plt.errorbar(th, promedio_amorfo, yerr=std_amorfo, capsize=0, 
              fmt='none', 
              # label=r'$\langle R_{exp} \rangle \pm \sigma_{\langle R_{exp} \rangle}$ Valor medio exp.',
              color='g',
              alpha=0.1)
plt.plot(th, promedio_amorfo, 'g')
prom_a_f = interp1d(th, promedio_amorfo)
plt.plot(theta_3, prom_a_f(theta_3),'og', label='Experimental')
plt.grid(None)
plt.xlim([0,90])
plt.ylim([0,1])
plt.title(f'$Sb_2S_3$ Amorfo', fontsize=18)
plt.legend(fontsize=14)
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia ($\lambda=633$ nm)', fontsize=18)


# CRISTALINO
plt.figure()

for k in range(np.size(d_cristalino)):
    if k == int(np.size(d_cristalino)/2):
        plt.plot(theta, R_p_cristalino[k],'brown', label = 'Simulación (100% polarización $p$)')
    else:
        plt.plot(theta, R_p_cristalino[k],'brown', alpha=0.05)       
# xerr = np.ones_like(std_cristalino)*2
plt.errorbar(th, promedio_cristalino, yerr=std_cristalino, capsize=0, 
              fmt='none', 
              # label=r'$\langle R_{exp} \rangle \pm \sigma_{\langle R_{exp} \rangle}$ Valor medio exp.',
              color='g',
              alpha=0.1)
plt.plot(th, promedio_cristalino, 'g')
prom_c_f = interp1d(th, promedio_cristalino)
plt.plot(thetac_3, prom_c_f(thetac_3),'og', label='Experimental')
plt.grid(None)
plt.xlim([0,90])
plt.ylim([0,1])
plt.title(f'$Sb_2S_3$ Cristalino', fontsize=18)
plt.legend(fontsize=14)
plt.xlabel(r"$\theta$ / grados", fontsize=18)
plt.ylabel('Reflectancia ($\lambda=633$ nm)', fontsize=18)
