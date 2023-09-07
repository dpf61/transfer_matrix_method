# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 19:26:12 2023

@author: dpf61
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from typing import Final
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors

import colour

from colour import XYZ_to_sRGB, XYZ_to_Lab, XYZ_to_Hunter_Lab, XYZ_to_RGB
from colour.models import RGB_COLOURSPACE_sRGB


import TransferMatrixMethod as tf

#%% COLECCIÓN DE MATERIALES

aire = tf.Material(np.ones(2000), np.zeros(2000), np.arange(0,2000), "aire")

data = np.loadtxt('SiO2.csv', skiprows=1, delimiter=';').transpose()
landas = data[0]
n = data[1] 
k = data[2]
nombre = "$SiO_2$"

vidrio = tf.Material(n, k, landas, nombre)

data = np.loadtxt('Crystalline_Sb2S3_Bari.csv', skiprows=1, delimiter=';').transpose()
landa = data[1]
n = data[4]
k = data[5]

cristalino = tf.Material(n, k, landa, "$Sb_2S_3$ Cristalino")

data = np.loadtxt('Amorphous_Sb2S3_Bari.csv', skiprows=1, delimiter=';').transpose()
landa = data[1]
n = data[4]
k = data[5]

amorfo = tf.Material(n, k, landa, "$Sb_2S_3$ Amorfo")

#%% ESPECTRAL

def grafica2d_espesores(x, y, R, etiqueta, vmin, vmax):
    fig, ax = plt.subplots()
    
    # Escalar los datos de R en el rango de 0 a 1
    norm = colors.Normalize(vmin, vmax)
    im = ax.imshow(R, cmap='inferno', norm=norm, extent=[min(x), max(x), min(y), max(y)], aspect='auto', origin='lower')

    # Configurar etiquetas de los ejes
    ax.set_xlabel('$\lambda$ / nm', fontsize=18)
    ax.set_ylabel('Espesor, $d$ / nm', fontsize=18)

    # Crear una barra de color para el mapa de escala de grises
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'{etiqueta}', fontsize=18)

    # Mostrar la gráfica
    plt.grid(None)


def main_espectral_cristalino():
    
    n = 300 # pixeles
    landa = np.linspace(420,800,n) #nm
    theta = 0 #grados
    pol = 'p'
    
    # d = np.linspace(245,265,3)*1e-9
    d = np.linspace(200,300,n)*1e-9

    R_pack = [[] for j in range(np.size(d))]
    T_pack = [[] for j in range(np.size(d))]
    A_pack = [[] for j in range(np.size(d))]
    
    for k in range(0,np.size(d)):
        
        R = []
        
        T = []
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(cristalino.get_n(landa[i]), cristalino.get_k(landa[i]), d[k], cristalino.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))
            
            T.append(stack.transmitancia(np.radians(theta), landa[i]*1e-9, pol))
            
        R_pack[k] = R
        T_pack[k] = T
        A_pack[k] = 1 - np.array(R) - np.array(T)
    
   
    plt.figure()
    grafica2d_espesores(landa, d*1e9, R_pack, 'Reflectancia', vmin=0, vmax=1)
    plt.title(r'$\rm Sb_2S_3$ Cristalino', fontsize=18)
    
    plt.figure()
    grafica2d_espesores(landa, d*1e9, T_pack, 'Transmitancia', vmin=0, vmax=1)
    plt.title(r'$\rm Sb_2S_3$ Cristalino', fontsize=18)
    
    plt.figure()
    grafica2d_espesores(landa, d*1e9, A_pack, 'Absorción', vmin=0, vmax=1)
    plt.title(r'$\rm Sb_2S_3$ Cristalino', fontsize=18)

    landa_borde = []
    for A in A_pack:
        landa_borde.append(landa[np.where(A<0.01)[0][0]])
    
    plt.plot(landa_borde,d*1e9, '--w')
    print(landa_borde)



def main_espectral_amorfo():

    landa = np.arange(420,800) #nm
    theta = 0 #grados
    pol = 'p'
    
    # d = np.linspace(280,300,3)*1e-9
    d = np.linspace(200,300,200)*1e-9

    R_pack = [[] for j in range(np.size(d))]
    T_pack = [[] for j in range(np.size(d))]
    A_pack = [[] for j in range(np.size(d))]
    
    for k in range(0,np.size(d)):
        
        R = []
        
        T = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(amorfo.get_n(landa[i]), amorfo.get_k(landa[i]), d[k], amorfo.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))
            
            T.append(stack.transmitancia(np.radians(theta), landa[i]*1e-9, pol))

        R_pack[k] = R
        T_pack[k] = T
        A_pack[k] = 1 - np.array(R) - np.array(T)
     
    plt.figure()
    grafica2d_espesores(landa, d*1e9, R_pack, 'Reflectancia', vmin=0, vmax=1)
    plt.title(r'$\rm Sb_2S_3$ Amorfo', fontsize=18)
    
    plt.figure()
    grafica2d_espesores(landa, d*1e9, T_pack, 'Transmitancia', vmin=0, vmax=1)
    plt.title(r'$\rm Sb_2S_3$ Amorfo', fontsize=18)
    
    plt.figure()
    grafica2d_espesores(landa, d*1e9, A_pack, 'Absorción', vmin=0, vmax=1)
    plt.title(r'$\rm Sb_2S_3$ Amorfo', fontsize=18)
    
    landa_borde = []
    for A in A_pack:
        landa_borde.append(landa[np.where(A<0.01)[0][0]])
    
    plt.plot(landa_borde,d*1e9, '--w')
    print(landa_borde)
    
    

main_espectral_amorfo()
main_espectral_cristalino()


#%% ANGULAR

def grafica2d_angulos(landa, angulos, R, magnitud, vmin, vmax):
    fig, ax = plt.subplots()

    # Escalar los datos de R en el rango de 0 a 1
    norm = colors.Normalize(vmin, vmax)
    im = ax.imshow(R, cmap='inferno', norm=norm, extent=[min(landa), max(landa), min(angulos), max(angulos)], aspect='auto', origin='lower')

    # Configurar etiquetas de los ejes
    ax.set_xlabel(r'$\theta$ / grados', fontsize=18)
    ax.set_ylabel('$\lambda$ / nm', fontsize=18)

    # Crear una barra de color para el mapa de escala de grises
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'{magnitud}', fontsize=18)

    # Mostrar la gráfica
    plt.grid(None)
    
    
def grafica2d_bwr(landa, angulos, R, magnitud, vmin, vmax):
    fig, ax = plt.subplots()

    # Escalar los datos de R en el rango de 0 a 1
    norm = colors.Normalize(vmin, vmax)
    im = ax.imshow(R, cmap='bwr', norm=norm, extent=[min(landa), max(landa), min(angulos), max(angulos)], aspect='auto', origin='lower')

    # Configurar etiquetas de los ejes
    ax.set_xlabel(r'$\theta$ / grados', fontsize=18)
    ax.set_ylabel('$\lambda$ / nm', fontsize=18)

    # Crear una barra de color para el mapa de escala de grises
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'{magnitud}', fontsize=18)

    # Mostrar la gráfica
    plt.grid(None)
    
def grafica2d(landa, angulos, R, etiqueta):
    fig, ax = plt.subplots()

    # Crear la gráfica 2D con el mapa de escala de grises
    im = ax.imshow(R, cmap='bwr', extent=[min(landa), max(landa), min(angulos), max(angulos)], aspect='auto')
    
    # Configurar etiquetas de los ejes
    ax.set_xlabel(r'$\theta$ / grados', fontsize=18)
    ax.set_ylabel('$\lambda$ / nm', fontsize=18)

    # Crear una barra de color para el mapa de escala de grises
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'{etiqueta}', fontsize=18)

    # Mostrar la gráfica
    plt.grid(None)
    

def main_angular_ambos():
    '''
    Método que representa R como colormap en función de lambda y theta.

    '''
    n = 300 #numero pixeles colormap
    landa = np.linspace(600,1500,n)
    theta = np.linspace(49.5,80,n) #grados
    pol = 'p'


    d_amorfo = 290*1e-9

    d_cristalino = 255*1e-9

    R_pack_cristalino = [[] for j in range(np.size(landa))]

    R_pack_amorfo = [[] for j in range(np.size(landa))]


    for k in range(0,np.size(landa)):
        
        R = [] 
        
        for i in range(0,np.size(theta)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(amorfo.get_n(landa[k]), amorfo.get_k(landa[k]), d_amorfo, amorfo.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[k]), vidrio.get_k(landa[k]), d_amorfo, vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta[i]), landa[k]*1e-9, pol))

        # guardo el array de R's en el pack
        R_pack_amorfo[k] = R


    for k in range(0,np.size(landa)):
        
        R = [] 
        
        for i in range(0,np.size(theta)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(cristalino.get_n(landa[k]), cristalino.get_k(landa[k]), d_cristalino, cristalino.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[k]), vidrio.get_k(landa[k]), d_cristalino, vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta[i]), landa[k]*1e-9, pol))

        # guardo el array de R's en el pack
        R_pack_cristalino[k] = R

    #Recogemos los Brewsters
    
    brewsters_cristalino = [[] for j in range(np.size(landa))]
    brewsters_amorfo = [[] for j in range(np.size(landa))]

    for k in range(np.size(landa)):
        brewsters_cristalino[k] = theta[np.argmin(R_pack_cristalino[k])]
        
    for k in range(np.size(landa)):
        brewsters_amorfo[k] = theta[np.argmin(R_pack_amorfo[k])]

    #Reflectancia amorfa
    plt.figure()
    grafica2d_angulos(theta, landa, R_pack_amorfo, 'Reflectancia', vmin=0, vmax=0.25)
    # plt.title(f'$Sb_2S_3$ Amorfo ($d$ = {d_amorfo*1e9:.0f}  nm)', fontsize=18)
    plt.title(f'$Sb_2S_3$ Amorfo', fontsize=18)
    plt.plot(brewsters_amorfo, landa,'--w')
    # plt.plot(brewsters_cristalino, landa,'-w')
    
    plt.figure()
    grafica2d_angulos(theta, landa, R_pack_cristalino, 'Reflectancia', vmin=0, vmax=0.25)
    # plt.title(f'$Sb_2S_3$ Cristalino ($d$ = {d_cristalino*1e9:.0f} nm)', fontsize=18)
    plt.title(f'$Sb_2S_3$ Cristalino', fontsize=18)
    # plt.plot(brewsters_amorfo, landa,'--w')
    plt.plot(brewsters_cristalino, landa,'-w')
                   
    
    R_log = [[] for j in range(np.size(landa))]
    for k in range(np.size(landa)):
        R_log[k] = 10*np.log10(np.array(R_pack_cristalino[k])/np.array(R_pack_amorfo[k]))
        
    plt.figure()
    grafica2d_bwr(theta, landa, R_log, r'$C$ / dB', vmin=-40, vmax=40)
    # plt.title(f'Contraste $Sb_2S_3$ ($d$ = {d_cristalino*1e9:.0f} nm)', fontsize=18)
    plt.title(f'Contraste $Sb_2S_3$', fontsize=18)
    plt.plot(brewsters_amorfo, landa,'--k')
    plt.plot(brewsters_cristalino, landa,'-k')
    


    '''
    Grafica del contraste a 55 grados
    '''
    
    
    # experimental 55 grados
    
    exp = np.loadtxt('amorfo_55.csv', skiprows=0, delimiter=';').transpose() 
    landa_exp_amorfo = exp[0]
    r55_exp_amorfo = exp[1]

    exp = np.loadtxt('cristalino_55.csv', skiprows=0, delimiter=';').transpose() 
    landa_exp_crist = exp[0]
    r55_exp_crist = exp[1]
    
    R_log_55_exp = 10*np.log10(np.array(r55_exp_crist)/np.array(r55_exp_amorfo))
    
    # calculado
    
    R_log = np.array(R_log)
    
    indice_55 = np.where((theta >= 55.0) & (theta <= 56))[0][0]
    
    
    R_log_55 = []
    
    for R in R_log:
        R_log_55.append(R[indice_55])
    
    print(theta[indice_55])
    
    plt.figure()
    plt.plot(landa, R_log_55, '-g' , label = 'Simulación')   
    plt.grid(None)
    plt.title(r'Contraste $Sb_2S_3$ ($\theta$ = 55$\degree$)', fontsize=18)
    plt.xlim(min(landa), max(landa))
    plt.ylim(-41,41)
    plt.plot(landa_exp_amorfo, R_log_55_exp, '-k', label = 'Experimental')
    plt.plot(landa, np.zeros_like(landa), '--k', linewidth=0.5)
    
    plt.ylabel(r'$C$ / dB', fontsize=18)
    plt.xlabel('$\lambda$ / nm', fontsize=18)
    plt.legend()
    

main_angular_ambos()

#%% EXPERIMENTAL ESPECTRAL

import numpy as np

color_crist = 'brown'

color_amorfo = 'darkblue'

def calcular_rmse(R_experimental, R_teorico, landa_experimental, landa_teorico, landa_0, landa_f):
    # Interpolar las funciones
    f_experimental = interp1d(landa_experimental, R_experimental)
    f_teorico = interp1d(landa_teorico, R_teorico)
    
    # Evaluar las funciones en el rango deseado
    landa_eval = np.linspace(landa_0, landa_f, num=1000)
    R_eval_experimental = f_experimental(landa_eval)
    R_eval_teorico = f_teorico(landa_eval)
    
    # Calcular el RMSE
    rmse = np.sqrt(np.mean((R_eval_experimental - R_eval_teorico)**2))
    
    return rmse



def main_exp_cristalino():
    
    landa = np.arange(380,800) #nm
    theta = 0 #grados
    pol = 'p'
    
    d = np.linspace(245,265,50)*1e-9

    R_pack = [[] for j in range(np.size(d))]
    
    for k in range(0,np.size(d)):
        
        R = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(cristalino.get_n(landa[i]), cristalino.get_k(landa[i]), d[k], cristalino.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))

        R_pack[k] = R

    ## FIGURA EXPERIMENTAL MEDIO
    plt.figure()
       
    exp = np.loadtxt('cristalino_exp_final.txt', skiprows=1, delimiter='\t').transpose() 
    landa_plot = exp[0]
    R_plot = exp[1]
    R_ste_plot = exp[2]*np.sqrt(9)
    
    plt.figure()
    plt.errorbar(landa_plot, R_plot, yerr=R_ste_plot, capsize=0, 
                  fmt='none', 
                  # label=r'$\langle R_{exp} \rangle \pm \sigma_{\langle R_{exp} \rangle}$ Valor medio exp.',
                  color=color_crist,
                  alpha=0.1)
    plt.plot(landa_plot, R_plot, color_crist, label='Experimental')
    plt.xlim([420,800])

    plt.xlabel(r'$\lambda$ / nm', fontsize='20')
    plt.ylabel('Reflectancia', fontsize='20')
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'$Sb_2S_3$ Cristalino', fontsize=18)
    plt.legend()
    plt.grid(None)
    plt.xlim([420,800])


    #FIGURA EXPERIMENTALES BRUTOS
    exp = np.loadtxt('R_cristalino_exp.txt', skiprows=1, delimiter='\t').transpose() 
    landa_exp = exp[0]
    
    R_cristalino = [[] for k in range(10)]
    for k in range(10):
        R_cristalino[k] = exp[k+1]
    
    plt.figure()
    for k in range(9):
        plt.plot(landa_exp, R_cristalino[k], alpha=0.5)
       
    plt.xlabel(r'$\lambda$ / nm', fontsize='20')
    plt.ylabel('Reflectancia', fontsize='20')
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'$Sb_2S_3$ Cristalino (medidas)', fontsize=18)
    plt.grid(None)
    plt.xlim([420,800])
    
    
    # FIGURA TEORICO 
    plt.figure()
    for R in R_pack:
        plt.plot(landa, R, color_crist, alpha=0.1)
    
    plt.plot(landa, R_pack[25], color_crist, label='Simulación')
    plt.xlabel(r'$\lambda$ / nm', fontsize='20')
    plt.ylabel('Reflectancia', fontsize='20')
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'$Sb_2S_3$ Cristalino', fontsize=18)       
    plt.legend()
    plt.grid(None)
    plt.xlim([420,800])
    
    # color_exp = tf.calcula_color_xyz(R_pack[25], landa, "D65")
    
    # color_teo = tf.calcula_color_xyz(R_plot, landa_plot, "D65")
    
    # CALCULO DEL RMSE
    
    RMSE = []
    for R in R_pack:
        RMSE.append(calcular_rmse(R_plot, R, landa_plot, landa, 650, 700))
        
    print(d[np.argmin(RMSE)])
    

# main_exp_cristalino()


def main_exp_amorfo():
    
    landa = np.arange(380,800) #nm
    theta = 0 #grados
    pol = 'p'
    
    d = np.linspace(280,300,50)*1e-9

    R_pack = [[] for j in range(np.size(d))]
    
    for k in range(0,np.size(d)):
        
        R = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(amorfo.get_n(landa[i]), amorfo.get_k(landa[i]), d[k], amorfo.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))

        R_pack[k] = R

    ## FIGURA EXPERIMENTAL MEDIO
    plt.figure()
       
    exp = np.loadtxt('amorfo_exp_final.txt', skiprows=1, delimiter='\t').transpose() 
    landa_plot = exp[0]
    R_plot = exp[1]
    R_ste_plot = exp[2]*np.sqrt(13)
    
    plt.figure()
    plt.errorbar(landa_plot, R_plot, yerr=R_ste_plot, capsize=0, 
                  fmt='none', 
                  # label=r'$\langle R_{exp} \rangle \pm \sigma_{\langle R_{exp} \rangle}$ Valor medio exp.',
                  color=color_amorfo,
                  alpha=0.1)
    plt.plot(landa_plot, R_plot, color_amorfo, label='Experimental')
    plt.xlim([420,800])

    plt.xlabel(r'$\lambda$ / nm', fontsize='20')
    plt.ylabel('Reflectancia', fontsize='20')
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'$Sb_2S_3$ Amorfo', fontsize=18)
    plt.legend()
    plt.grid(None)
    plt.xlim([420,800])


    #FIGURA EXPERIMENTALES BRUTOS
    exp = np.loadtxt('R_amorfo_exp.txt', skiprows=1, delimiter='\t').transpose() 
    landa_exp = exp[0]
    
    R_amorfo = [[] for k in range(10)]
    for k in range(10):
        R_amorfo[k] = exp[k+1]
    
    plt.figure()
    for k in range(9):
        plt.plot(landa_exp, R_amorfo[k], alpha=0.5)
       
    plt.xlabel(r'$\lambda$ / nm', fontsize='20')
    plt.ylabel('Reflectancia', fontsize='20')
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'$Sb_2S_3$ Amorfo (medidas)', fontsize=18)
    plt.grid(None)
    plt.xlim([420,800])
    
    
    # FIGURA TEORICO 
    plt.figure()
    for R in R_pack:
        plt.plot(landa, R, color_amorfo, alpha=0.1)
    
    plt.plot(landa, R_pack[25], color_amorfo, label='Simulación')
    plt.xlabel(r'$\lambda$ / nm', fontsize='20')
    plt.ylabel('Reflectancia', fontsize='20')
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'$Sb_2S_3$ Amorfo', fontsize=18)       
    plt.legend()
    plt.grid(None)
    plt.xlim([420,800])
    
    # color_exp = tf.calcula_color_xyz(R_pack[25], landa, "D65")
    
    # color_teo = tf.calcula_color_xyz(R_plot, landa_plot, "D65")
    
    # CALCULO DEL RMSE
    
    RMSE = []
    for R in R_pack:
        RMSE.append(calcular_rmse(R_plot, R, landa_plot, landa, 650, 700))
        
    print(d[np.argmin(RMSE)])
   
   
main_exp_cristalino()

main_exp_amorfo()


#%% EXPERIMENTAL INGLÉS

def main_exp_cristalino():
    
    #FIGURA EXPERIMENTALES BRUTOS
    exp = np.loadtxt('R_cristalino_exp.txt', skiprows=1, delimiter='\t').transpose() 
    landa_exp = exp[0]
    
    R_cristalino = [[] for k in range(10)]
    for k in range(10):
        R_cristalino[k] = exp[k+1]
    
    plt.figure()
    for k in range(9):
        plt.plot(landa_exp, R_cristalino[k], label=f'R measurement {k+1}', alpha=0.5)
       
    plt.xlabel(r'Wavelength (nm)', fontsize=18)
    plt.ylabel(r'Reflectance ($\theta = 0 \degree$)', fontsize=18)
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'Crystalline $Sb_2S_3$', fontsize=18)
    plt.grid(None)
    plt.xlim([420,800])
    plt.legend(fontsize=8)


def main_exp_amorfo():
    
    #FIGURA EXPERIMENTALES BRUTOS
    exp = np.loadtxt('R_amorfo_exp.txt', skiprows=1, delimiter='\t').transpose() 
    landa_exp = exp[0]
    
    R_amorfo = [[] for k in range(10)]
    for k in range(10):
        R_amorfo[k] = exp[k+1]
    
    plt.figure()
    for k in range(9):
        plt.plot(landa_exp, R_amorfo[k], label=f'R measurement {k+1}', alpha=0.5)
       
    plt.xlabel(r'Wavelength (nm)', fontsize=18)
    plt.ylabel(r'Reflectance ($\theta = 0 \degree$)', fontsize=18)
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'Amorphous $Sb_2S_3$', fontsize=18)
    plt.grid(None)
    plt.xlim([420,800])
    plt.legend(fontsize=8)
    
    
main_exp_amorfo()
main_exp_cristalino()


#%% EXPERIMENTAL ESPAÑOL

def main_exp_cristalino():
    
    #FIGURA EXPERIMENTALES BRUTOS
    exp = np.loadtxt('R_cristalino_exp.txt', skiprows=1, delimiter='\t').transpose() 
    landa_exp = exp[0]
    
    R_cristalino = [[] for k in range(10)]
    for k in range(10):
        R_cristalino[k] = exp[k+1]
    
    plt.figure()
    for k in range(9):
        plt.plot(landa_exp, R_cristalino[k], label=f'Medida {k+1}', alpha=0.5)
       
    plt.xlabel(r'$\lambda$ / nm', fontsize=18)
    plt.ylabel(r'Reflectancia ($\theta = 0 \degree$)', fontsize=18)
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'$Sb_2S_3$ Cristalino', fontsize=18)
    plt.grid(None)
    plt.xlim([420,800])
    plt.legend(fontsize=8)


def main_exp_amorfo():
    
    #FIGURA EXPERIMENTALES BRUTOS
    exp = np.loadtxt('R_amorfo_exp.txt', skiprows=1, delimiter='\t').transpose() 
    landa_exp = exp[0]
    
    R_amorfo = [[] for k in range(10)]
    for k in range(10):
        R_amorfo[k] = exp[k+1]
    
    plt.figure()
    for k in range(9):
        plt.plot(landa_exp, R_amorfo[k], label=f'Medida {k+1}', alpha=0.5)
       
    plt.xlabel(r'$\lambda$ / nm', fontsize=18)
    plt.ylabel(r'Reflectancia ($\theta = 0 \degree$)', fontsize=18)
    plt.ylim([0,1])
    plt.ylim([-0.05,1.05])
    plt.title(f'$Sb_2S_3$ Amorfo', fontsize=18)
    plt.grid(None)
    plt.xlim([420,800])
    plt.legend(fontsize=8)
    
    
main_exp_amorfo()
main_exp_cristalino()


#%% GRAFICA 55 GRADOS CON ERROR
    

def main_55_ambos():
    
    n = 300 #numero pixeles colormap
    landa = np.linspace(210,1500,n)
    theta = 55.0 #grados
    pol = 'p'


    d_amorfo = np.array([100,150,200,255,290,350,500])*1e-9

    d_cristalino = d_amorfo
        

    R_pack_cristalino = [[] for j in range(np.size(landa))]

    R_pack_amorfo = [[] for j in range(np.size(landa))]


    for k in range(0,np.size(d_amorfo)):
        
        R = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(amorfo.get_n(landa[i]), amorfo.get_k(landa[i]), d_amorfo[k], amorfo.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d_amorfo[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))

        # guardo el array de R's en el pack
        R_pack_amorfo[k] = R


    for k in range(0,np.size(d_cristalino)):
        
        R = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(cristalino.get_n(landa[i]), cristalino.get_k(landa[i]), d_cristalino[k], cristalino.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d_cristalino[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))

        # guardo el array de R's en el pack
        R_pack_cristalino[k] = R
                   
    
    R_log = [[] for j in range(np.size(d_cristalino))]
    for k in range(np.size(d_cristalino)):
        R_log[k] = 10*np.log10(np.array(R_pack_cristalino[k])/np.array(R_pack_amorfo[k]))
        

    '''
    Grafica del contraste a 55 grados
    '''
    
    # calculado
    
    plt.figure()
    
    for k in range(np.size(d_amorfo)):
        plt.plot(landa, R_log[k], label=f'd = {d_amorfo[k]*1e9:.0f} nm')
        

    plt.grid(None)
    plt.title(r'Contraste $Sb_2S_3$ ($\theta$ = 55$\degree$)', fontsize=18)
    plt.xlim(min(landa), max(landa))
    plt.ylim(-41,41)
    
    plt.ylabel(r'$C$ / dB', fontsize=18)
    plt.xlabel('$\lambda$ / nm', fontsize=18)
    plt.legend()
    plt.plot(landa, np.zeros_like(landa), '--k', linewidth=0.5)


main_55_ambos()


#%% COLORMAP CONTRASTE ESPESOR

def grafica2d_bwr(espesores, landa, R, magnitud, vmin, vmax):
    fig, ax = plt.subplots()

    # Escalar los datos de R en el rango de 0 a 1
    norm = colors.Normalize(vmin, vmax)
    im = ax.imshow(R, cmap='bwr', norm=norm, extent=[min(espesores), max(espesores), min(landa), max(landa)], aspect='auto', origin='lower')

    # Configurar etiquetas de los ejes
    ax.set_ylabel(r'$d$ / nm', fontsize=18)
    ax.set_xlabel('$\lambda$ / nm', fontsize=18)

    # Crear una barra de color para el mapa de escala de grises
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'{magnitud}', fontsize=18)

    # Mostrar la gráfica
    plt.grid(None)

def colormap_contraste():
    
    n = 400 #numero pixeles colormap
    landa = np.linspace(210,1500,n)
    theta = 55.0 #grados
    pol = 'p'


    d_amorfo = np.linspace(100,500,n)*1e-9

    d_cristalino = d_amorfo
        

    R_pack_cristalino = [[] for j in range(np.size(landa))]

    R_pack_amorfo = [[] for j in range(np.size(landa))]


    for k in range(0,np.size(d_amorfo)):
        
        R = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(amorfo.get_n(landa[i]), amorfo.get_k(landa[i]), d_amorfo[k], amorfo.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d_amorfo[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))

        # guardo el array de R's en el pack
        R_pack_amorfo[k] = R


    for k in range(0,np.size(d_cristalino)):
        
        R = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(cristalino.get_n(landa[i]), cristalino.get_k(landa[i]), d_cristalino[k], cristalino.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d_cristalino[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))

        # guardo el array de R's en el pack
        R_pack_cristalino[k] = R
                   
    
    R_log = [[] for j in range(np.size(d_cristalino))]
    for k in range(np.size(d_cristalino)):
        R_log[k] = 10*np.log10(np.array(R_pack_cristalino[k])/np.array(R_pack_amorfo[k]))
        
        
    plt.figure()
    grafica2d_bwr(landa, d_cristalino*1e9, R_log, r'$C$ / dB', vmin=-40, vmax=40)
    plt.title(r'Contraste $Sb_2S_3$ ($\theta$ = 55$\degree$)', fontsize=18)
        
    

    # '''
    # Grafica del contraste a 55 grados
    # '''
    
    # # calculado
    
    # plt.figure()
    
    # for k in range(np.size(d_amorfo)):
    #     plt.plot(landa, R_log[k], label=f'd = {d_amorfo[k]*1e9:.0f} nm')
        

    # plt.grid(None)
    # plt.title(r'Contraste $Sb_2S_3$ ($\theta$ = 55$\degree$)', fontsize=18)
    # plt.xlim(min(landa), max(landa))
    # plt.ylim(-41,41)
    
    # plt.ylabel(r'$C$ / dB', fontsize=18)
    # plt.xlabel('$\lambda$ / nm', fontsize=18)
    # plt.legend()
    # plt.plot(landa, np.zeros_like(landa), '--k', linewidth=0.5)


colormap_contraste()




#%% CONSTATES SB2S3


def main_sb2s3():
    
    data = np.loadtxt('Crystalline_Sb2S3_Bari.csv', skiprows=1, delimiter=';').transpose()
    landa = data[1]
    n = data[4]
    k = data[5]
    
    cristalino = tf.Material(n, k, landa, 'Cristalino')
    
    data = np.loadtxt('Amorphous_Sb2S3_Bari.csv', skiprows=1, delimiter=';').transpose()
    landa = data[1]
    n = data[4]
    k = data[5]
    
    amorfo = tf.Material(n, k, landa, 'Amorfo')
    
    color_cris = 'brown'

    color_amor = 'darkblue'
   
    landas = np.linspace(400,1500,1000)
    
    n_cris = cristalino.get_n(landas)
    n_amor = amorfo.get_n(landas)
    
    k_cris = cristalino.get_k(landas)
    k_amor = amorfo.get_k(landas)
    
    plt.figure()
    plt.plot(landas, n_amor, color=color_amor, label=r'$n_{am}$')
    plt.plot(landas, n_cris, color=color_cris, label=r'$n_{cr}$')
    plt.plot(landas, k_amor, '--', color=color_amor, label=r'$k_{am}$', alpha=0.5)
    plt.plot(landas, k_cris, '--', color=color_cris, label=r'$k_{cr}$', alpha=0.5)

    plt.ylim([-0.1,3.1])
    plt.xlim(min(landas),max(landas))
    plt.title(r"$\rm Sb_2S_3$", fontsize=18)
    plt.grid(None)    
    plt.legend()
    plt.xlabel('$\lambda$ / nm', fontsize=18)
    plt.ylabel('Índice de refracción complejo', fontsize=18)
    
    n_contrast = n_cris - n_amor
    k_contrast = k_cris - k_amor
    
    plt.figure()
    plt.plot(landas, n_contrast, color='green', label=r'$\Delta n$')
    plt.plot(landas, k_contrast, color='purple', label=r'$\Delta k$')

    plt.ylim([-0.01,0.7])
    plt.xlim(min(landas),max(landas))
    plt.title(r"$\rm Sb_2S_3$", fontsize=18)
    plt.grid(None)    
    plt.legend()
    plt.xlabel('$\lambda$ / nm', fontsize=18)
    plt.ylabel('Contraste de índice', fontsize=18) 


main_sb2s3()


#%% COLORES

data = np.loadtxt('Iluminante_A.txt', skiprows=0, delimiter='\t').transpose()
landas_a = data[0]
iluminante_a = data[1] 

data = np.loadtxt('Iluminante_D65.txt', skiprows=0, delimiter='\t').transpose()
landas_d65 = data[0]
iluminante_d65 = data[1] 

data = np.loadtxt('Observador_2_grados_CIE_1931.txt', skiprows=0, delimiter='\t').transpose()
landas_obs = data[0]
obs_1 = data[1] 
obs_2 = data[2] 
obs_3 = data[3] 


def XYZ_to_xyz(color_XYZ):
    
    X = color_XYZ[0]
    Y = color_XYZ[1]
    Z = color_XYZ[2]
    
    x=X/(X+Y+Z)
    y=Y/(X+Y+Z)
    z=Z/(X+Y+Z)
    
    return [x,y,z]



def calcula_xyzlab(reflectancia, landas, iluminante):

    if iluminante == "A":
        landas_ilu = landas_a
        espectro_ilum = iluminante_a
    elif iluminante == "D65":
        landas_ilu = landas_d65
        espectro_ilum = iluminante_d65
    else:
        print('No existe el iluminante')
        
    espectro_iluminante_interp = interp1d(landas_ilu, espectro_ilum)
    espectro_iluminante = espectro_iluminante_interp(landas_obs)
    
    espectro_R_interp = interp1d(landas, reflectancia)
    espectro_R = espectro_R_interp(landas_obs)

    K = 100/sum(espectro_iluminante*obs_2)
    X = K*sum(espectro_iluminante * espectro_R * obs_1)
    Y = K*sum(espectro_iluminante * espectro_R * obs_2)
    Z = K*sum(espectro_iluminante * espectro_R * obs_3)
    
    x=X/(X+Y+Z)
    y=Y/(X+Y+Z)
    z=Z/(X+Y+Z)

    # Coordenadas Lab
    ref_X = K * np.sum(espectro_iluminante * np.ones(espectro_iluminante.shape) * obs_1)
    ref_Y = K * np.sum(espectro_iluminante * np.ones(espectro_iluminante.shape) * obs_2)
    ref_Z = K * np.sum(espectro_iluminante * np.ones(espectro_iluminante.shape) * obs_3)
    
    var_X = X / ref_X
    var_Y = Y / ref_Y
    var_Z = Z / ref_Z
    
    if var_X > 216 / 24389:
        var_X = var_X ** (1 / 3)
    else:
        var_X = (841 / 108 * var_X) + (16 / 116)
    
    if var_Y > 216 / 24389:
        var_Y = var_Y ** (1 / 3)
    else:
        var_Y = (841 / 108 * var_Y) + (16 / 116)
    
    if var_Z > 216 / 24389:
        var_Z = var_Z ** (1 / 3)
    else:
        var_Z = (841 / 108 * var_Z) + (16 / 116)
    
    L = (116 * var_Y) - 16
    a = 500 * (var_X - var_Y)
    b = 200 * (var_Y - var_Z)
    C = math.sqrt(a**2 + b**2)
    h = math.atan2(b, a)
    
    return ([x,y,z],[L,a,b,C,h])



def distancia_color(lab_teo,lab_exp):
    L_t = lab_teo[0]
    a_t = lab_teo[1]
    b_t = lab_teo[2]
    
    L_e = lab_exp[0]
    a_e = lab_exp[1]
    b_e = lab_exp[2]
    
    return np.sqrt( (L_t - L_e)**2 + (a_t - a_e)**2 + (b_t - b_e)**2 )

    
def main_colores_cristalino():

    # intervalo de longitudes de onda que nos interesa:
        # más tarde esto tendrá que meterlo el usuario
    landa = np.arange(380,800) #nm
    theta = 0 #grados
    pol = 'p'
    
    # d = np.linspace(250,260)*1e-9
    d = np.linspace(280,300,200)*1e-9
    
    # creamos una lista vacia donde iremos guardando R(lambda, d) par cada espesor d
    R_pack = [[] for j in range(np.size(d))]
    
    #iteramos calculando las reflectancias para cada espesor y cada lambda
    for k in range(0,np.size(d)):
        
        R = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(cristalino.get_n(landa[i]), cristalino.get_k(landa[i]), d[k], cristalino.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))
    
        # guardo el array de R's en el pack
        R_pack[k] = R
   
    plt.figure()
    plt.style.use(['science', 'notebook', 'grid'])
    
    colores_xyz = [[] for i in range(np.size(d))]
    
    for i in range(np.size(d)):
        
        colores_xyz[i] = tf.calcula_color(R_pack[i], landa, "D65")
    
    plt.figure()
    tf.grafica_color(d*1e9, colores_xyz)
    # plot_thickness_color(d*1e9, colores_xyz)
    plt.title('Color $Sb_2S_3$ Cristalino (Iluminante D65)', fontsize=18)
    print(f'CRISTALINO, Para D65, d= {d[25]*1e9:.0f} nm el color xyz = {colores_xyz[25]}')
        
    
    # datos colorimetro XYZ mayusculas para cristalino DIA 1 FONDO BLANCO (los mejores)
    a = [21.19,21.36,21.62]
    b = [21.05,21.21,21.38]
    c = [20.63,20.78,20.58]
    d = [20.46,20.56,20.42]
    e = [21.76,21.99,22.08]
    color_exp_cristalino= [a,b,c,d,e]
    

    # # datos colorimetro XYZ mayusculas para cristalino DIA 2 FONDO NEGRO
    # a = [18.46,19.44,20.14]
    # b = [18.42,19.45,20.09]
    # c = [18.22,19.15,20.02]
    # d = [18.25,19.16,19.94]
    # e = [18.33,19.25,20.04]
    # color_exp_cristalino= [a,b,c,d,e]
    

    for color in color_exp_cristalino:
        print(f'CRIS EXP {XYZ_to_xyz(color)}') 
    
    color_exp_cris_xyz = [[] for j in range(5)]  
    for i in range(5):
        color_exp_cris_xyz[i] = XYZ_to_xyz(color_exp_cristalino[i])    
    
    plt.figure()
    tf.grafica_color([1], [color_exp_cris_xyz[0]])
    
    lab_teo = calcula_xyzlab(R_pack[25], landa, "D65")[0]
    lab_teo = XYZ_to_Lab(lab_teo)
    
    lab_exp = XYZ_to_Lab(color_exp_cris_xyz[0])
    
    print(f'CIELAB teorico 1: {lab_teo}')
    print(f'CIELAB colorimetro 1: {lab_exp}')
    
    print(f'La distancia es {distancia_color(lab_teo,lab_exp)}')
    

main_colores_cristalino()
    
def main_colores_amorfo():

    # intervalo de longitudes de onda que nos interesa:
        # más tarde esto tendrá que meterlo el usuario
    landa = np.arange(380,800) #nm
    theta = 0 #grados
    pol = 'p'
    
    d = np.linspace(280,300,200)*1e-9
    
    # creamos una lista vacia donde iremos guardando R(lambda, d) par cada espesor d
    R_pack = [[] for j in range(np.size(d))]
    
    #iteramos calculando las reflectancias para cada espesor y cada lambda
    for k in range(0,np.size(d)):
        
        R = [] 
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(amorfo.get_n(landa[i]), amorfo.get_k(landa[i]), d[k], amorfo.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))
    
        R_pack[k] = R


   
    plt.figure()
    plt.style.use(['science', 'notebook', 'grid'])
    
    colores_xyz = [[] for i in range(np.size(d))]
    
    for i in range(np.size(d)):
        
        colores_xyz[i] = tf.calcula_color(R_pack[i], landa, "D65")
    
    plt.figure()
    tf.grafica_color(d*1e9, colores_xyz)
    plt.title('Color $Sb_2S_3$ Amorfo (Iluminante D65)', fontsize=18)
    print(f'Para D65, d= {d[25]*1e9:.0f} nm el color xyz = {colores_xyz[25]}')


    # # datos colorimetro XYZ mayusculas para amorfo DIA 1  FONDO BLANCO
    # a = [40.87, 42.78, 13.07]
    # b = [39.61, 41.68, 12.94]
    # c = [41.00, 42.81, 13.00]
    # d = [41.53, 43.20, 13.25]
    # e = [41.24, 43.19, 13.22]
    
    # color_exp_amorfo= [a,b,c,d,e]

    # datos colorimetro XYZ mayusculas para amorfo DIA 2 FONDO NEGRO
    a = [14.62,16.48,8.95]
    b = [14.41,16.27,9.43]
    c = [14.46,16.52,9.96]
    d = [14.28,15.69,9.62]
    e = [14.22,15.90,9.61]
    color_exp_amorfo= [a,b,c,d,e]
    
    for color in color_exp_amorfo:
        print(f'AMORFO EXP {XYZ_to_xyz(color)}') 
    
    color_exp_amorfo_xyz = [[] for j in range(5)]  
    for i in range(5):
        color_exp_amorfo_xyz[i] = XYZ_to_xyz(color_exp_amorfo[i])    
    
    plt.figure()
    tf.grafica_color([1], [color_exp_amorfo_xyz[0]])
    
    lab_teo = calcula_xyzlab(R_pack[25], landa, "D65")[0]
    lab_teo = XYZ_to_Lab(lab_teo)
    
    lab_exp = XYZ_to_Lab(color_exp_amorfo_xyz[0])
    
    print(f'CIELAB teorico 1: {lab_teo}')
    print(f'CIELAB colorimetro 1: {lab_exp}')
    
    print(f'La distancia es {distancia_color(lab_teo,lab_exp)}')
    
    
main_colores_amorfo()

#%% GRAFICA R 55 GRADOS ADICIONAL

# experimental 55 grados

exp = np.loadtxt('amorfo_55.csv', skiprows=0, delimiter=';').transpose() 
landa_exp_amorfo = exp[0]
r55_exp_amorfo = exp[1]

exp = np.loadtxt('cristalino_55.csv', skiprows=0, delimiter=';').transpose() 
landa_exp_crist = exp[0]
r55_exp_crist = exp[1]


color_cris = 'brown'

color_amor = 'darkblue'


def r_55_amor():
    
    landa = np.arange(400,1500) #nm
    theta = 55 #grados
    pol = 'p'
    
    d = np.linspace(280,300)*1e-9

    R_pack = [[] for j in range(np.size(d))]
    T_pack = [[] for j in range(np.size(d))]
    A_pack = [[] for j in range(np.size(d))]
    
    for k in range(0,np.size(d)):
        
        R = []
        
        T = []
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(amorfo.get_n(landa[i]), amorfo.get_k(landa[i]), d[k], amorfo.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))
            
            T.append(stack.transmitancia(np.radians(theta), landa[i]*1e-9, pol))
            
        R_pack[k] = R
        T_pack[k] = T
        A_pack[k] = 1 - np.array(R) - np.array(T)
    
   
    plt.figure()
    for R in R_pack:
        plt.plot(landa, R, color=color_amor, alpha=0.05)
        
    plt.plot(landa, R_pack[25], color=color_amor, label='Simulación')
    plt.plot(landa_exp_amorfo, r55_exp_amorfo, '--k', label='Experimental')
        
    
    plt.xlabel(r'$\lambda$ / nm', fontsize=18)
    plt.ylabel(r'Reflectancia', fontsize=18)
    plt.ylim([0,1])
    plt.ylim([-0.01,0.5])
    plt.title(r'$Sb_2S_3$ Amorfo ($\theta$ = 55 $\degree$)', fontsize=18)
    plt.grid(None)
    plt.xlim([400,1500])
    plt.legend(fontsize=18)
    
    

def r_55_cris():
    
    landa = np.arange(400,1500) #nm
    theta = 55 #grados
    pol = 'p'
    
    d = np.linspace(245,265)*1e-9

    R_pack = [[] for j in range(np.size(d))]
    T_pack = [[] for j in range(np.size(d))]
    A_pack = [[] for j in range(np.size(d))]
    
    for k in range(0,np.size(d)):
        
        R = []
        
        T = []
        
        for i in range(0,np.size(landa)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(cristalino.get_n(landa[i]), cristalino.get_k(landa[i]), d[k], cristalino.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[i]), vidrio.get_k(landa[i]), d[k], vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta), landa[i]*1e-9, pol))
            
            T.append(stack.transmitancia(np.radians(theta), landa[i]*1e-9, pol))
            
        R_pack[k] = R
        T_pack[k] = T
        A_pack[k] = 1 - np.array(R) - np.array(T)
    
   
    plt.figure()
    for R in R_pack:
        plt.plot(landa, R, color=color_cris, alpha=0.05)
        
    plt.plot(landa, R_pack[25], color=color_cris, label='Simulación')
    plt.plot(landa_exp_crist, r55_exp_crist, '--k', label='Experimental')
        
    
    plt.xlabel(r'$\lambda$ / nm', fontsize=18)
    plt.ylabel(r'Reflectancia', fontsize=18)
    plt.ylim([0,1])
    plt.ylim([-0.01,0.5])
    plt.title(r'$Sb_2S_3$ Cristalino ($\theta$ = 55 $\degree$)', fontsize=18)
    plt.grid(None)
    plt.xlim([400,1500])
    plt.legend(fontsize=18)


r_55_cris()
r_55_amor()




#%% CONTRASTE 55 DOS SIMULACIONES

def contraste_55_simulado(d_amorfo,d_cristalino):

    n = 300 #numero pixeles colormap
    landa = np.arange(600,1500)
    theta = np.arange(54,56) #grados
    pol = 'p'

    d_amorfo = d_amorfo*1e-9

    d_cristalino = d_cristalino*1e-9

    R_pack_cristalino = [[] for j in range(np.size(landa))]

    R_pack_amorfo = [[] for j in range(np.size(landa))]


    for k in range(0,np.size(landa)):
        
        R = [] 
        
        for i in range(0,np.size(theta)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(amorfo.get_n(landa[k]), amorfo.get_k(landa[k]), d_amorfo, amorfo.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[k]), vidrio.get_k(landa[k]), d_amorfo, vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta[i]), landa[k]*1e-9, pol))

        # guardo el array de R's en el pack
        R_pack_amorfo[k] = R


    for k in range(0,np.size(landa)):
        
        R = [] 
        
        for i in range(0,np.size(theta)):
            
            stack = tf.Multicapa()
            
            stack.añade_capa(cristalino.get_n(landa[k]), cristalino.get_k(landa[k]), d_cristalino, cristalino.get_nombre)

            stack.añade_capa(vidrio.get_n(landa[k]), vidrio.get_k(landa[k]), d_cristalino, vidrio.get_nombre)
            
            R.append(stack.reflectancia(np.radians(theta[i]), landa[k]*1e-9, pol))

        # guardo el array de R's en el pack
        R_pack_cristalino[k] = R

         
    
    R_log = [[] for j in range(np.size(landa))]
    for k in range(np.size(landa)):
        R_log[k] = 10*np.log10(np.array(R_pack_cristalino[k])/np.array(R_pack_amorfo[k]))
        
    R_log = np.array(R_log)
    
    indice_55 = np.where((theta >= 55.0) & (theta <= 56))[0][0]
    
    R_log_55 = []
    
    for R in R_log:
        R_log_55.append(R[indice_55])
        
    return landa, R_log_55


def main_dos_55():
    landa, contraste_UNO = contraste_55_simulado(290,255)
    landa, contraste_DOS = contraste_55_simulado(298,243)
    # landa, contraste_DOS = contraste_55_simulado(300,245)

    # experimental 55 grados

    exp = np.loadtxt('amorfo_55.csv', skiprows=0, delimiter=';').transpose() 
    landa_exp_amorfo = exp[0]
    r55_exp_amorfo = exp[1]

    exp = np.loadtxt('cristalino_55.csv', skiprows=0, delimiter=';').transpose() 
    landa_exp_crist = exp[0]
    r55_exp_crist = exp[1]

    R_log_55_exp = 10*np.log10(np.array(r55_exp_crist)/np.array(r55_exp_amorfo))

    plt.figure()
    plt.plot(landa, contraste_UNO, '-g' , label = 'Simulación', alpha=0.7)
    plt.plot(landa, contraste_DOS, '-b' , label = 'Simulación (min. RMSE)', alpha=0.5)    
    plt.grid(None)
    plt.title(r'Contraste $Sb_2S_3$ ($\theta$ = 55$\degree$)', fontsize=18)
    plt.xlim(min(landa), max(landa))
    plt.ylim(-41,41)
    plt.plot(landa_exp_amorfo, R_log_55_exp, '-k', label = 'Experimental')
    plt.plot(landa, np.zeros_like(landa), '--k', linewidth=0.5)

    plt.ylabel(r'$C$ / dB', fontsize=18)
    plt.xlabel('$\lambda$ / nm', fontsize=18)
    plt.legend()
    
main_dos_55()



