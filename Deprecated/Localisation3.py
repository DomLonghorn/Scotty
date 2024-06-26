# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 08:37:48 2020

@author: VH Chen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy import interpolate as interpolate
import math
from scipy import constants as constants






def find_H_Cardano(K_magnitude,launch_angular_frequency,epsilon_para,epsilon_perp,epsilon_g,theta_m):  
    # This function is designed to be evaluated in post-procesing, hence it uses different inputs from the usual H

    wavenumber_K0 = launch_angular_frequency / constants.c
    n_ref_index = K_magnitude/wavenumber_K0
    sin_theta_m = np.sin(theta_m)
    cos_theta_m = np.cos(theta_m)
    
    
    D_11_component = epsilon_perp - n_ref_index**2*sin_theta_m**2
    D_22_component = epsilon_perp - n_ref_index**2
    D_bb_component = epsilon_para - n_ref_index**2*cos_theta_m**2
    D_12_component = epsilon_g
    D_1b_component = n_ref_index**2*sin_theta_m*cos_theta_m
    
    h_2_coefficient = - D_11_component - D_22_component - D_bb_component
#    h_2_coefficient = 0
    h_1_coefficient = D_11_component*D_bb_component + D_11_component*D_22_component + D_22_component*D_bb_component - D_12_component**2 - D_1b_component**2
#    h_0_coefficient = 0
    h_0_coefficient = D_22_component*D_1b_component**2 + D_bb_component*D_12_component**2 - D_11_component*D_22_component*D_bb_component
    
    h_t_coefficient = (
                            -2*h_2_coefficient**3
                            +9*h_2_coefficient*h_1_coefficient
                            -27*h_0_coefficient
                            +3*np.sqrt(3)*np.sqrt(
                                4*h_2_coefficient**3 * h_0_coefficient
                                -h_2_coefficient**2 * h_1_coefficient**2
                                -18*h_2_coefficient * h_1_coefficient * h_0_coefficient
                                +4*h_1_coefficient**3
                                +27*h_0_coefficient**2
                                +0j #to make the argument of the np.sqrt complex, so that the sqrt evaluates negative functions
                            )
                        )**(1/3)
    
    H_1_Cardano = h_t_coefficient/(3*2**(1/3)) - 2**(1/3) *(3*h_1_coefficient - h_2_coefficient**2)/(3*h_t_coefficient) - h_2_coefficient/3
    H_2_Cardano = - (1 - 1j*np.sqrt(3))/(6*2**(1/3))*h_t_coefficient + (1 + 1j*np.sqrt(3))*(3*h_1_coefficient - h_2_coefficient**2)/(3*2**(2/3)*h_t_coefficient) - h_2_coefficient/3
    H_3_Cardano = - (1 + 1j*np.sqrt(3))/(6*2**(1/3))*h_t_coefficient + (1 - 1j*np.sqrt(3))*(3*h_1_coefficient - h_2_coefficient**2)/(3*2**(2/3)*h_t_coefficient) - h_2_coefficient/3
    return H_1_Cardano, H_2_Cardano, H_3_Cardano

def find_coefficients(K_magnitude,launch_angular_frequency,epsilon_para,epsilon_perp,epsilon_g,theta_m):  
    # This function is designed to be evaluated in post-procesing, hence it uses different inputs from the usual H

    wavenumber_K0 = launch_angular_frequency / constants.c
    n_ref_index = K_magnitude/wavenumber_K0
    sin_theta_m = np.sin(theta_m)
    cos_theta_m = np.cos(theta_m)
    
    
    D_11_component = epsilon_perp - n_ref_index**2*sin_theta_m**2
    D_22_component = epsilon_perp - n_ref_index**2
    D_bb_component = epsilon_para - n_ref_index**2*cos_theta_m**2
    D_12_component = epsilon_g
    D_1b_component = n_ref_index**2*sin_theta_m*cos_theta_m
    
    h_2_coefficient = - D_11_component - D_22_component - D_bb_component
    h_1_coefficient = D_11_component*D_bb_component + D_11_component*D_22_component + D_22_component*D_bb_component - D_12_component**2 - D_1b_component**2
    h_0_coefficient = D_22_component*D_1b_component**2 + D_bb_component*D_12_component**2 - D_11_component*D_22_component*D_bb_component
    
    return h_0_coefficient, h_1_coefficient, h_2_coefficient

def find_D_matrix(K_magnitude,launch_angular_frequency,epsilon_para,epsilon_perp,epsilon_g,theta_m):
    wavenumber_K0 = launch_angular_frequency / constants.c
    n_ref_index = K_magnitude/wavenumber_K0
    sin_theta_m = np.sin(theta_m)
    cos_theta_m = np.cos(theta_m)
    
    
    D_11_component = epsilon_perp - n_ref_index**2*sin_theta_m**2
    D_22_component = epsilon_perp - n_ref_index**2
    D_bb_component = epsilon_para - n_ref_index**2*cos_theta_m**2
    D_12_component = epsilon_g
    D_1b_component = n_ref_index**2*sin_theta_m*cos_theta_m
    
    D_matrix = np.zeros([3,3],dtype='complex128')
    
    D_matrix[0][0] = D_11_component
    D_matrix[1][1] = D_22_component
    D_matrix[2][2] = D_bb_component
    D_matrix[0][2] = D_1b_component
    D_matrix[2][0] = D_1b_component
    D_matrix[0][1] = -1j*D_12_component
    D_matrix[1][0] = 1j*D_12_component

    return D_matrix

def find_dH_dKR_Cardano(K_magnitude, launch_angular_frequency, epsilon_para,epsilon_perp, epsilon_g, theta_m, delta_K_R):
    
    H_plus  = find_H_Cardano(K_magnitude,launch_angular_frequency,epsilon_para,epsilon_perp,epsilon_g,theta_m)
    H_minus = find_H(q_R, q_Z, K_R-delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    dH_dKR = (H_plus - H_minus) / (2 * delta_K_R)
    
    return dH_dKR


def find_dH_dKZ_Cardano(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_Z, 
                interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus  = find_H(q_R, q_Z, K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus = find_H(q_R, q_Z, K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    dH_dKZ = (H_plus - H_minus) / (2 * delta_K_Z)
    
    return dH_dKZ


def find_dH_dKzeta_Cardano(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_zeta, 
                   interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus  = find_H(q_R, q_Z, K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus = find_H(q_R, q_Z, K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    dH_dKzeta = (H_plus - H_minus) / (2 * delta_K_zeta)
    
    return dH_dKzeta













loadfile = np.load('data_input.npz')
launch_freq_GHz =loadfile['launch_freq_GHz']
loadfile.close()

loadfile = np.load('analysis_output.npz')
localisation_piece = loadfile['localisation_piece']
cutoff_index = loadfile['cutoff_index']
RZ_distance_along_line = loadfile['RZ_distance_along_line']
distance_along_line = loadfile['distance_along_line']
k_perp_1_backscattered = loadfile['k_perp_1_backscattered']
Psi_xx_output = loadfile['Psi_xx_output']
Psi_xy_output = loadfile['Psi_xy_output']
Psi_yy_output = loadfile['Psi_yy_output']
M_xx_output = loadfile['M_xx_output']
M_xy_output = loadfile['M_xy_output']
M_yy_output = loadfile['M_yy_output']
in_index = loadfile['in_index']
out_index = loadfile['out_index']
theta_m_output = loadfile['theta_m_output']
loadfile.close()

loadfile = np.load('data_output.npz')
g_magnitude_output = loadfile['g_magnitude_output']
q_R_array = loadfile['q_R_array']
K_R_array = loadfile['K_R_array']
K_zeta_initial = loadfile['K_zeta_initial']
K_Z_array = loadfile['K_Z_array']
epsilon_para_output = loadfile['epsilon_para_output']
epsilon_perp_output = loadfile['epsilon_perp_output']
epsilon_g_output = loadfile['epsilon_g_output']
loadfile.close()

delta_K_R = 0.1 #in the same units as K_R
delta_K_zeta = 0.1 #in the same units as K_zeta
delta_K_Z = 0.1 #in the same units as K_z



launch_angular_frequency = 2*math.pi*10.0**9 * launch_freq_GHz
K_magnitude_array = np.sqrt(K_R_array**2 + K_Z_array**2 + K_zeta_initial**2/q_R_array**2)

K_magnitude_array_plus_KR     = np.sqrt((K_R_array+delta_K_R)**2 + K_Z_array**2 + K_zeta_initial**2/q_R_array**2)
K_magnitude_array_minus_KR    = np.sqrt((K_R_array-delta_K_R)**2 + K_Z_array**2 + K_zeta_initial**2/q_R_array**2)
K_magnitude_array_plus_Kzeta  = np.sqrt(K_R_array**2 + K_Z_array**2 + (K_zeta_initial+delta_K_zeta)**2/q_R_array**2)
K_magnitude_array_minus_Kzeta = np.sqrt(K_R_array**2 + K_Z_array**2 + (K_zeta_initial-delta_K_zeta)**2/q_R_array**2)
K_magnitude_array_plus_KZ     = np.sqrt(K_R_array**2 + (K_Z_array+delta_K_Z)**2 + K_zeta_initial**2/q_R_array**2)
K_magnitude_array_minus_KZ    = np.sqrt(K_R_array**2 + (K_Z_array-delta_K_Z)**2 + K_zeta_initial**2/q_R_array**2)

numberOfDataPoints = len(K_magnitude_array)

H_1_Cardano_array = np.zeros(numberOfDataPoints,dtype='complex128')
H_2_Cardano_array = np.zeros(numberOfDataPoints,dtype='complex128')
H_3_Cardano_array = np.zeros(numberOfDataPoints,dtype='complex128')

H_1_Cardano_array_numerical = np.zeros(numberOfDataPoints,dtype='complex128')
H_2_Cardano_array_numerical = np.zeros(numberOfDataPoints,dtype='complex128')
H_3_Cardano_array_numerical = np.zeros(numberOfDataPoints,dtype='complex128')

H_1_Cardano_array_numerical2 = np.zeros(numberOfDataPoints,dtype='complex128')
H_2_Cardano_array_numerical2 = np.zeros(numberOfDataPoints,dtype='complex128')
H_3_Cardano_array_numerical2 = np.zeros(numberOfDataPoints,dtype='complex128')

h_0_coefficient_array = np.zeros(numberOfDataPoints)
h_1_coefficient_array = np.zeros(numberOfDataPoints)
h_2_coefficient_array = np.zeros(numberOfDataPoints)

for ii in range(0,numberOfDataPoints):
    K_magnitude = K_magnitude_array[ii]
    epsilon_para = epsilon_para_output[ii]
    epsilon_perp = epsilon_perp_output[ii]
    epsilon_g = epsilon_g_output[ii]
    theta_m = theta_m_output[ii]
    
    H_1_Cardano_array[ii],H_2_Cardano_array[ii],H_3_Cardano_array[ii] = find_H_Cardano(K_magnitude,launch_angular_frequency,epsilon_para,epsilon_perp,epsilon_g,theta_m)

    h_0_coefficient_array[ii],h_1_coefficient_array[ii],h_2_coefficient_array[ii] = find_coefficients(K_magnitude,launch_angular_frequency,epsilon_para,epsilon_perp,epsilon_g,theta_m)

    [H_1_Cardano_array_numerical2[ii],H_2_Cardano_array_numerical2[ii],H_3_Cardano_array_numerical2[ii]] = np.roots([1,h_2_coefficient_array[ii],h_1_coefficient_array[ii],h_0_coefficient_array[ii]])
    
    D_matrix = find_D_matrix(K_magnitude,launch_angular_frequency,epsilon_para,epsilon_perp,epsilon_g,theta_m) 
    
    [H_1_Cardano_array_numerical[ii],H_2_Cardano_array_numerical[ii],H_3_Cardano_array_numerical[ii]] = np.linalg.eigvals(D_matrix)
    
plt.figure()
plt.plot(distance_along_line,abs(H_1_Cardano_array),'r')    
plt.plot(distance_along_line,abs(H_2_Cardano_array),'g')    
plt.plot(distance_along_line,abs(H_3_Cardano_array),'b')    
plt.ylim(0,1)

plt.figure()
plt.plot(distance_along_line,abs(H_1_Cardano_array_numerical),'r')    
plt.plot(distance_along_line,abs(H_2_Cardano_array_numerical),'g')    
plt.plot(distance_along_line,abs(H_3_Cardano_array_numerical),'b')    
plt.ylim(0,1)

plt.figure()
plt.plot(distance_along_line,abs(H_1_Cardano_array_numerical2),'r')    
plt.plot(distance_along_line,abs(H_2_Cardano_array_numerical2),'g')    
plt.plot(distance_along_line,abs(H_3_Cardano_array_numerical2),'b')    
plt.ylim(0,1)


plt.figure()
plt.plot(distance_along_line,h_0_coefficient_array,'r')    
plt.plot(distance_along_line,h_1_coefficient_array,'g')    
plt.plot(distance_along_line,h_2_coefficient_array,'b')    

plt.figure()
plt.plot(distance_along_line,abs(H_1_Cardano_array))    
plt.plot(distance_along_line,abs(H_1_Cardano_array_numerical))    
plt.plot(distance_along_line,abs(H_1_Cardano_array_numerical2))    
