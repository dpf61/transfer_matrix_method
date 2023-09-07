# -*- coding: utf-8 -*-
"""
Created on Tue May 23 12:34:08 2023

@author: dpf61

COLORIMETRIA
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import colour
import math

# import scienceplots

# colour.plotting.colour_style()

# plt.style.use({'figure.figsize': (10.24, 5.76)})

# plt.style.use('science')
plt.style.use(['science', 'notebook', 'grid'])



data = np.loadtxt('Iluminante_A.txt', skiprows=0, delimiter='\t').transpose()
landas_a = np.round(data[0])
iluminante_a = data[1] 


plt.figure()
plt.plot(landas_a, iluminante_a)
plt.title('Iluminante A')

data = np.loadtxt('Iluminante_D65.txt', skiprows=0, delimiter='\t').transpose()
landas_d65 = np.round(data[0])
iluminante_d65 = data[1] 


plt.figure()
plt.plot(landas_d65, iluminante_d65)
plt.title('Iluminante D65')

data = np.loadtxt('Observador_2_grados_CIE_1931.txt', skiprows=0, delimiter='\t').transpose()
landas_obs = np.round(data[0])
obs_1 = data[1] 
obs_2 = data[2] 
obs_3 = data[3] 


plt.figure()
plt.plot(landas_obs, obs_1, 'r', label=r'$\bar{x}(\lambda)$') 
plt.plot(landas_obs, obs_2, 'g', label=r'$\bar{y}(\lambda)$')
plt.plot(landas_obs, obs_3, 'b', label=r'$\bar{z}(\lambda)$')
plt.title('Funciones de igualación del color', fontsize=20)
plt.ylabel('Respuesta relativa', fontsize=20)
plt.xlabel(r'$\lambda$ / nm', fontsize=20)
plt.xlim([min(landas_obs), max(landas_obs)])
plt.ylim([-0.05, 1.85])
plt.legend(fontsize=20)
plt.savefig('funciones_dos_grados.png', dpi=600)

def calcula_color(reflectancia, landas, iluminante):
    
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
    X = K*sum(espectro_iluminante*espectro_R*obs_1)
    Y = K*sum(espectro_iluminante*espectro_R*obs_2)
    Z = K*sum(espectro_iluminante*espectro_R*obs_3)
    
    x=X/(X+Y+Z)
    y=Y/(X+Y+Z)
    z=Z/(X+Y+Z)
    
    return [x,y,z]


R = np.array([0.2, 0.5, 0.8])
l = np.array([300, 550, 900])
color_xyz=calcula_color(R, l, "A")

print(color_xyz)


def xyz_to_rgb(color_xyz):
    
    matriz = np.array([[0.4174,-0.1576,-0.0828],[-0.0909,0.2522,0.0157],[0.0009,-0.0025,0.1786]])
    
    return matriz@color_xyz
    
print(xyz_to_rgb(color_xyz)) 


def dibujar_espectro(lambdas, espectro, title):
    # Crear un diccionario con los valores de lambda y espectro
    espectro_dict = dict(zip(lambdas, espectro))
    
    # Crear un objeto de espectro a partir del diccionario
    espectro_obj = colour.SpectralDistribution(espectro_dict)
    
    # Dibujar el espectro
    fig, ax = colour.plotting.plot_single_sd(espectro_obj, out_of_gamut_clipping =True, title=title)
    
    # # Mostrar el gráfico
    # plt.show()
    
    return fig, ax

# Ejemplo de uso
lambdas = landas_d65
espectro = iluminante_d65

fig, ax = dibujar_espectro(lambdas, espectro, 'Iluminante')
ax.set_ylabel('Emisión espectral relativa', fontsize=20)
ax.set_xlabel(r'$\lambda$ / nm', fontsize=20)
fig.savefig('D65_ilu_big.png', dpi=600)

# lambdas = landas_a
# espectro = iluminante_a

# dibujar_espectro(lambdas, espectro, 'Iluminante A')




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


def grafica_color(eje_y, xyz_colors):
    # Convierte los colores xyz a coordenadas RGB para poder graficarlos
    rgb_colors = [colour.XYZ_to_sRGB(xyz_color) for xyz_color in xyz_colors]
    # Crea la gráfica
    for i in range(len(eje_y)):
        plt.hlines(eje_y[i], 0, 1, colors=[rgb_colors[i]], linewidth=10)
    plt.xlabel('Color')
    plt.ylabel('Espesor')
    plt.xlim(0,1)
    plt.ylim(np.min(eje_y), np.max(eje_y))
    plt.show()