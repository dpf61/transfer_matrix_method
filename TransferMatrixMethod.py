# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 21:55:23 2023

@author: dpf61

TRANSFER MATRIX METHOD
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from typing import Final
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D

plt.style.use(['science', 'notebook', 'grid'])

# Importamos la librería colour necesaria para transformar colores a RGB de
# forma que Matplotlib los reconozca en el estándar adecuado
import colour
from colour import XYZ_to_sRGB

'''
    CLASES
'''


class Material:
    '''
    Clase que implementa un material caracterizado por sus constantes ópticas:
        índice de refracción n(lambda) y coeficiente de extinción k(lambda).
    '''
    
    def __init__(self, n: list[float], k: list[float], landas: list[float], nombre: str):
        '''
        Método constructor de la clase Material

        Parameters
        ----------
        n : list[float]
            Valores n(lambda) de la parte real del índice de refracción.
            Unidades: adimensional.
        k : list[float]
            Valores k(lambda) de la parte imaginaria del índice de refracción.
            Unidades: adimensional.
        landas : list[float]
            Longitudes de onda (lambdas) de las funciones n(lambda) y
            k(lambda) introducidas.
            Unidades: NANÓMETROS.
        nombre : str
            Nombre del compuesto químico. Introducir fórmula química entre $$.

        Returns
        -------
        Objeto de la clase Material.

        '''
        self.n = n
        self.k = k
        self.landas = landas
        self.complexN = self.n - 1j*self.k
        self.nombre = nombre
        
    def get_n(self, landa: list[float]):
        
        # interpolacion que devuelve n en función de la landa
        funcion = interp1d(self.landas, self.n)
        
        return funcion(landa) #list
    
    def get_k(self, landa: list[float]):
        
        # interpolacion que devuelve n en función de la landa
        funcion = interp1d(self.landas, self.k)
        
        return funcion(landa) #list
    
    def get_nombre(self):
        return self.nombre
        
    def grafica_nk(self, intervalo_landas: list[float]):
        
        plt.plot(intervalo_landas, self.get_n(intervalo_landas), label=f'n($\lambda$) {self.nombre}')
        plt.plot(intervalo_landas, self.get_k(intervalo_landas), label=f'k($\lambda$) {self.nombre}')
        plt.xlabel(r'$\lambda$/nm', fontsize=18)
        plt.xlim([min(intervalo_landas), max(intervalo_landas)])
        plt.grid(None)
        
    def grafica_n(self, intervalo_landas: list[float]):
        
        plt.plot(intervalo_landas, self.get_n(intervalo_landas), label=f'{self.nombre}')
        plt.xlabel(r'$\lambda$/nm', fontsize=18)
        plt.xlim([min(intervalo_landas), max(intervalo_landas)])
        plt.ylabel(r'Índice de refracción, $n$', fontsize=18)
        plt.grid(None)

    def grafica_k(self, intervalo_landas: list[float]):
        
        plt.plot(intervalo_landas, self.get_k(intervalo_landas), label=f'{self.nombre}')
        plt.xlabel(r'$\lambda$/nm', fontsize=18)
        plt.xlim([min(intervalo_landas), max(intervalo_landas)])
        plt.ylabel(r'Coeficiente de extinción, $k$', fontsize=18)
        plt.grid(None)


class ThinFilm:
    '''
    Clase que implementa una thin film caracterizada por su índice de refracción,
    su coeficiente de extinción y su espesor.
    '''
    
    def __init__(self, n: float, k: float, d: float, nombre: str):
        self.n = n
        self.k = k
        self.complexN = self.n - 1j*self.k # + o - ik REVISARR
        self.d = d
        self.nombre = nombre
    
    def matriz_transfer(self, theta: float, landa: float, pol: str):
        """
        Parameters
        ----------
        theta : float (RADIANES)
            Ángulo de incidencia sobre esa capa. Si es una capa distinta a la
            inicial, ese ángulo vendrá dado por la ley de Snell.
        landa : float (METROS)
            Longitud de onda de la luz incidente.
        pol : str 
            Polarización de la luz incidente. "p" si paralela al PI y "s" si
            perpendicular
            
        Returns
        -------
        Matriz 2x2 de esa thin film en particular.
        """
        # Constantes:
        # Y_0: Admitancia del vacio, en Siemens o 1/Ohm        
        Y_0: Final[float] = 2.6544e-3
        
        
        # theta = np.radians(theta) # LO QUITO PORQUE DA ERROR AL METER UN INPUT COMPLEJO
        
        # Desfase experimentado por la onda al atravesar la capa
        delta = 2*np.pi * self.complexN * self.d * np.cos(theta)/landa
        
        if pol == "s":
            eta = Y_0*self.complexN*np.cos(theta)  
        else: # pol = "p"
            eta = Y_0*self.complexN/np.cos(theta)
            
        M = np.array([[np.cos(delta), 1j*np.sin(delta)/eta], [1j*eta*np.sin(delta), np.cos(delta)]])
        
        return M

class Multicapa:
    
    def __init__(self):
        self.multicapa = [] # Lista vacía para guardar los objetos ThinFilm
        
    def añade_capa(self, n: float, k: float, d: float, nombre: str):
        # Este método añade una ThinFilm al array de forma "manual"
        self.multicapa.append(ThinFilm(n, k, d, nombre))

    def crear_capas(self, N):
        # for i in range(N):
        #     # Crear un objeto ThinFilm con un nombre y un grosor aleatorios
        #     nombre = f"Layer {i+1}"
        #     thickness = np.random.randint(10, 100)
        #     capa = ThinFilm(1.5, 0.1, thickness, nombre)
        #     # Añadir el objeto a la lista del atributo layers
        #     self.capas.append(capa)
        
        # añado este código que pide las caracteristicas de las TF al usuario
        for i in range(0,N):
            print('Introduzca las características de la ThinFilm ', i)
            n = float(input("Índice de refracción: "))
            k = float(input("Coeficiente de extinción: "))
            d = float(input("Espesor (en metros): "))
            nombre = str(input("Nombre del material: "))
            self.multicapa.append(ThinFilm(n, k, d, nombre))
            
    def eta_m(self, theta_0: float, pol: str):
        # es util un metodo que devuelva el valor de eta_m
        Y_0: Final[float] = 2.6544e-3
        
        i = np.size(self.multicapa)-1 #ultimo indice, el del sustrato
        # Sacamos theta_m por snell
        theta_m = np.arcsin(1*np.sin(theta_0)/self.multicapa[i].complexN)
        
        if pol == "s":
            eta_m = Y_0*self.multicapa[i].complexN*np.cos(theta_m)
        else: 
            eta_m = Y_0*self.multicapa[i].complexN/np.cos(theta_m)
            
        return eta_m

    def tmm(self, theta: float, landa: float, pol: str):
        '''
        
        Parameters
        ----------
        theta : float (RADIANES)
            DESCRIPTION. 
        landa : float
            DESCRIPTION.
        pol : str
            DESCRIPTION.

        Returns
        -------
        [B, C]
            Array de dos componentes

        '''
        
        Y_0: Final[float] = 2.6544e-3
        
        # La admitancia inclinada viene dada por el angulo de incidencia y la polarización.
        # En el caso del vacío/aire:
        if pol == "s":
            eta_0 = Y_0*np.cos(theta)
        else:
            eta_0 = Y_0/np.cos(theta)
        
        producto = [[1,0],[0,1]] #matriz identidad como si fuera la del aire/vacío
        
        N = np.size(self.multicapa) #número de capas para el siguiente bucle
        
        for i in range(0, N):
            
            # Actualizamos el ángulo de incidencia en la capa mediante ley de Snell
            # DEBERIA METODOLOGIZAR ESTE PASO, QUEDA MUY SUCIO, ADEMAS EN LA EC (2.97)
            # PARECE MÁS SENCILLO CALCULAR LOS THETA_R DESDE THETA_0
            if i == 0:
                theta = np.arcsin( 1*np.sin(theta) / self.multicapa[i].complexN)
            else:
                theta = np.arcsin(self.multicapa[i-1].complexN * np.sin(theta) / self.multicapa[i].complexN)
            
            # Realizamos la multiplicación
            if i is not (N-1):
                producto = producto @ self.multicapa[i].matriz_transfer(theta, landa, pol)
            else: #es la capa de sustrato
                if pol == "s":
                    eta_m = Y_0*self.multicapa[i].complexN*np.cos(theta)
                else: 
                    eta_m = Y_0*self.multicapa[i].complexN/np.cos(theta)
                producto = producto @ [1, eta_m]
        
        return producto # devuelve la matriz de transferencia [B, C]
        
    def reflectancia(self, theta: float, landa: float, pol: str):
        
        Y_0: Final[float] = 2.6544e-3 #CAMBIAR ESTO Y PONERLO COMO UN SELF DE LA CLASE
        
        # La admitancia inclinada viene dada por el angulo de incidencia y la polarización.
        # En el caso del vacío/aire:
        if pol == "s":
            eta_0 = Y_0*np.cos(theta)
        else:
            eta_0 = Y_0/np.cos(theta)
            
        #Cálculo de R:
        
        tmm = self.tmm(theta, landa, pol)
        B = tmm[0]
        C = tmm[1]
        R = np.abs( (eta_0*B-C)/(eta_0*B+C) )**2
        # R = (eta_0*B-C)/(eta_0*B+C) * np.conjugate((eta_0*B-C)/(eta_0*B+C))
        
        return R
    
    def transmitancia(self, theta: float, landa: float, pol: str):
        
        Y_0: Final[float] = 2.6544e-3 #CAMBIAR ESTO Y PONERLO COMO UN SELF DE LA CLASE
        
        # La admitancia inclinada viene dada por el angulo de incidencia y la polarización.
        # En el caso del vacío/aire:
        if pol == "s":
            eta_0 = Y_0*np.cos(theta)
        else:
            eta_0 = Y_0/np.cos(theta)
            
        #Cálculo de R:
        
        tmm = self.tmm(theta, landa, pol)
        B = tmm[0]
        C = tmm[1]
        eta_m = self.eta_m(theta, pol)
        T = 4*eta_0*np.real(eta_m)/np.abs(eta_0*B+C)**2
        
        return T



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

def calcula_color_xyz(reflectancia, landas, iluminante):
    '''
    Método basado en la teoría de la página 126 de Fundamentos de Colorimetría.
    
    '''  
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
    
    # print(f'El valor de triestimulo es {[X,Y,Z]}')
    
    x=X/(X+Y+Z)
    y=Y/(X+Y+Z)
    z=Z/(X+Y+Z)
    
    # print(f'El valor de coordenada es {[x,y,z]}')
    
    return [x,y,z]


def grafica_color(espesores, colores_xyz):
    # Crear la figura y los ejes
    fig, ax = plt.subplots()

    # Generar líneas horizontales de colores en la gráfica
    for espesor, color_xyz in zip(espesores, colores_xyz):
        # Convertir de XYZ a RGB
        color_rgb = XYZ_to_sRGB(color_xyz, colour.CCS_ILLUMINANTS['CIE 1931 2 Degree Standard Observer']['D65'])
        # illuminant = np.array([0.34570, 0.35850])
        # color_rgb = XYZ_to_RGB(color_xyz, colour.RGB_Colourspace('Adobe RGB (1998)'), illuminant, None)
        # color_rgb = xyz2rgb(color_xyz)


        # Generar línea horizontal de color en la gráfica
        ax.axhline(y=espesor, color=color_rgb, linewidth = 2)

    # Configurar los límites y etiquetas de los ejes
    ax.set_xlim(0, 1)  # Límites del eje x
    ax.set_ylim(np.min(espesores), np.max(espesores))  # Límites del eje y
    # ax.set_xlabel('Color')  # Etiqueta del eje x
    ax.set_xticks([])
    ax.set_ylabel(r'Espesor, $d$ / nm', fontsize=18)  # Etiqueta del eje y
    
    
    
    
def calcula_color(reflectancia, landas, iluminante):
    '''
    Método basado en la teoría de la página 126 de Fundamentos de Colorimetría.
    
    '''  
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
    
    # print(f'El valor de triestimulo es {[X,Y,Z]}')
    
    x=X/(X+Y+Z)
    y=Y/(X+Y+Z)
    z=Z/(X+Y+Z)
    
    # print(f'El valor de coordenada es {[x,y,z]}')
    
    return [x,y,z]