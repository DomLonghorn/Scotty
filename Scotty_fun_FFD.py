# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 10:44:34 2018

Functions for finding derivatives of H using forward finite difference.

@author: chenv
Valerian Hongjie Hall-Chen
valerian@hall-chen.com

Run in Python 3,  does not work in Python 2
"""

#import numpy as np
#from scipy import constants as constants
from Scotty_fun_general import find_H


def find_dH_dR(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_R, 
               interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_0 = find_H(q_R, 
                 q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_1 = find_H(q_R+delta_R, 
                 q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_2 = find_H(q_R+2*delta_R, 
                 q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    dH_dR = ( (-3/2)*H_0 + (2)*H_1 + (-1/2)*H_2 ) / (delta_R)
    
    return dH_dR


def find_dH_dZ(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_Z, 
               interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_0 = find_H(q_R, q_Z, 
                 K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_1 = find_H(q_R, q_Z+delta_Z, 
                 K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_2 = find_H(q_R, q_Z+2*delta_Z, 
                 K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    dH_dZ = ( (-3/2)*H_0 + (2)*H_1 + (-1/2)*H_2 ) / (delta_Z)
    
    return dH_dZ


def find_dH_dKR(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_R, 
                interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus  = find_H(q_R, q_Z, K_R+delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus = find_H(q_R, q_Z, K_R-delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    dH_dKR = (H_plus - H_minus) / (2 * delta_K_R)
    
    return dH_dKR


def find_dH_dKZ(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_Z, 
                interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus  = find_H(q_R, q_Z, K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus = find_H(q_R, q_Z, K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    dH_dKZ = (H_plus - H_minus) / (2 * delta_K_Z)
    
    return dH_dKZ


def find_dH_dKzeta(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_zeta, 
                   interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus  = find_H(q_R, q_Z, K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus = find_H(q_R, q_Z, K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    dH_dKzeta = (H_plus - H_minus) / (2 * delta_K_zeta)
    
    return dH_dKzeta
      

def find_d2H_dR2(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_R, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_0 = find_H(q_R, 
                 q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_1 = find_H(q_R+delta_R, 
                 q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_2 = find_H(q_R+2*delta_R, 
                 q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_3 = find_H(q_R+3*delta_R, 
                 q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dR2 = ( (2)*H_0 + (-5)*H_1 + (4)*H_2 + (-1)*H_3 ) / (delta_R**2)
    
    return d2H_dR2


def find_d2H_dZ2(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_Z, 
                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus    = find_H(q_R, q_Z+delta_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_current = find_H(q_R, q_Z,         K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus   = find_H(q_R, q_Z-delta_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dZ2 = (H_plus - 2*H_current + H_minus) / (delta_Z**2)
    
    return d2H_dZ2


def find_d2H_dR_dZ(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_R, delta_Z, 
                   interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_R_plus_Z   = find_H(q_R+delta_R, q_Z+delta_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                               interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_R_minus_Z  = find_H(q_R+delta_R, q_Z-delta_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                               interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z) 
    H_minus_R_plus_Z  = find_H(q_R-delta_R, q_Z+delta_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                               interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_R_minus_Z = find_H(q_R-delta_R, q_Z-delta_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                               interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dR_dZ = (H_plus_R_plus_Z - H_plus_R_minus_Z - H_minus_R_plus_Z + H_minus_R_minus_Z) / (4 * delta_R * delta_Z)
    
    return d2H_dR_dZ


def find_d2H_dKR2(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_R, 
                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus    = find_H(q_R, q_Z, K_R+delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_current = find_H(q_R, q_Z, K_R,           K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus   = find_H(q_R, q_Z, K_R-delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKR2 = (H_plus - 2*H_current + H_minus) / (delta_K_R**2)
    
    return d2H_dKR2


def find_d2H_dKR_dKzeta(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_R, delta_K_zeta, 
                        interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_R_plus_K_zeta   = find_H(q_R, q_Z, K_R+delta_K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_R_minus_K_zeta  = find_H(q_R, q_Z, K_R+delta_K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_R_plus_K_zeta  = find_H(q_R, q_Z, K_R-delta_K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_R_minus_K_zeta = find_H(q_R, q_Z, K_R-delta_K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKR_dKzeta = (H_plus_K_R_plus_K_zeta - H_plus_K_R_minus_K_zeta - H_minus_K_R_plus_K_zeta + H_minus_K_R_minus_K_zeta) / (4*delta_K_R*delta_K_zeta)
    
    return d2H_dKR_dKzeta


def find_d2H_dKR_dKz(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_R, delta_K_Z, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_R_plus_K_z   = find_H(q_R, q_Z, K_R+delta_K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                                   interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_R_minus_K_z  = find_H(q_R, q_Z, K_R+delta_K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                                   interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_R_plus_K_z  = find_H(q_R, q_Z, K_R-delta_K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                                   interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_R_minus_K_z = find_H(q_R, q_Z, K_R-delta_K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                                   interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKR_dKz = (H_plus_K_R_plus_K_z - H_plus_K_R_minus_K_z - H_minus_K_R_plus_K_z + H_minus_K_R_minus_K_z) / (4*delta_K_R*delta_K_Z)
    
    return d2H_dKR_dKz


def find_d2H_dKzeta2(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_zeta, 
                     interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus    = find_H(q_R, q_Z, K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_current = find_H(q_R, q_Z, K_R, K_zeta,            K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus   = find_H(q_R, q_Z, K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKzeta2 = (H_plus - 2*H_current + H_minus) / (delta_K_zeta**2)    
    
    return d2H_dKzeta2

def find_d2H_dKzeta_dKz(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_zeta, delta_K_Z, 
                        interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_zeta_plus_K_z   = find_H(q_R, q_Z, K_R, K_zeta+delta_K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_zeta_minus_K_z  = find_H(q_R, q_Z, K_R, K_zeta+delta_K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_zeta_plus_K_z  = find_H(q_R, q_Z, K_R, K_zeta-delta_K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_zeta_minus_K_z = find_H(q_R, q_Z, K_R, K_zeta-delta_K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKzeta_dKz = (H_plus_K_zeta_plus_K_z - H_plus_K_zeta_minus_K_z - H_minus_K_zeta_plus_K_z + H_minus_K_zeta_minus_K_z) / (4*delta_K_zeta*delta_K_Z)    
    
    return d2H_dKzeta_dKz


def find_d2H_dKZ2(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_Z, 
                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus    = find_H(q_R, q_Z, K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_current = find_H(q_R, q_Z, K_R, K_zeta, K_Z,           launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus   = find_H(q_R, q_Z, K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKZ2 = (H_plus - 2*H_current + H_minus) / (delta_K_Z**2)   
    
    return d2H_dKZ2


def find_d2H_dKR_dR(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_R, delta_R, 
                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_R_plus_R   = find_H(q_R+delta_R, q_Z, K_R+delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_R_minus_R  = find_H(q_R-delta_R, q_Z, K_R+delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_R_plus_R  = find_H(q_R+delta_R, q_Z, K_R-delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_R_minus_R = find_H(q_R-delta_R, q_Z, K_R-delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKR_dR = (H_plus_K_R_plus_R - H_plus_K_R_minus_R - H_minus_K_R_plus_R + H_minus_K_R_minus_R) / (4*delta_K_R*delta_R)
    
    return d2H_dKR_dR
 
    
def find_d2H_dKR_dZ(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_R, delta_Z, 
                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_R_plus_Z   = find_H(q_R, q_Z+delta_Z, K_R+delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_R_minus_Z  = find_H(q_R, q_Z-delta_Z, K_R+delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_R_plus_Z  = find_H(q_R, q_Z+delta_Z, K_R-delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_R_minus_Z = find_H(q_R, q_Z-delta_Z, K_R-delta_K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKR_dZ = (H_plus_K_R_plus_Z - H_plus_K_R_minus_Z - H_minus_K_R_plus_Z + H_minus_K_R_minus_Z) / (4*delta_K_R*delta_Z)
    
    return d2H_dKR_dZ
  
    
def find_d2H_dKzeta_dR(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_zeta, delta_R, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_zeta_plus_R   = find_H(q_R+delta_R, q_Z, K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_zeta_minus_R  = find_H(q_R-delta_R, q_Z, K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_zeta_plus_R  = find_H(q_R+delta_R, q_Z, K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_zeta_minus_R = find_H(q_R-delta_R, q_Z, K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKzeta_dR = (H_plus_K_zeta_plus_R - H_plus_K_zeta_minus_R - H_minus_K_zeta_plus_R + H_minus_K_zeta_minus_R) / (4*delta_K_zeta*delta_R)  
    
    return d2H_dKzeta_dR

 
def find_d2H_dKzeta_dZ(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_zeta, delta_Z, 
                       interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_zeta_plus_Z   = find_H(q_R, q_Z+delta_Z, K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_zeta_minus_Z  = find_H(q_R, q_Z-delta_Z, K_R, K_zeta+delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_zeta_plus_Z  = find_H(q_R, q_Z+delta_Z, K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_zeta_minus_Z = find_H(q_R, q_Z-delta_Z, K_R, K_zeta-delta_K_zeta, K_Z, launch_angular_frequency, mode_flag, 
                                  interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)    
    d2H_dKzeta_dZ = (H_plus_K_zeta_plus_Z - H_plus_K_zeta_minus_Z - H_minus_K_zeta_plus_Z + H_minus_K_zeta_minus_Z) / (4*delta_K_zeta*delta_Z)  
    
    return d2H_dKzeta_dZ


def find_d2H_dKZ_dR(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_Z, delta_R, 
                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_Z_plus_R   = find_H(q_R+delta_R, q_Z, K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_Z_minus_R  = find_H(q_R-delta_R, q_Z, K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_Z_plus_R  = find_H(q_R+delta_R, q_Z, K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_Z_minus_R = find_H(q_R-delta_R, q_Z, K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKZ_dR = (H_plus_K_Z_plus_R - H_plus_K_Z_minus_R - H_minus_K_Z_plus_R + H_minus_K_Z_minus_R) / (4*delta_K_Z*delta_R)  
    
    return d2H_dKZ_dR


def find_d2H_dKZ_dZ(q_R, q_Z, K_R, K_zeta, K_Z, launch_angular_frequency, mode_flag, delta_K_Z, delta_Z, 
                    interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z):
    
    H_plus_K_Z_plus_Z   = find_H(q_R, q_Z+delta_Z, K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_plus_K_Z_minus_Z  = find_H(q_R, q_Z-delta_Z, K_R, K_zeta, K_Z+delta_K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_Z_plus_Z  = find_H(q_R, q_Z+delta_Z, K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    H_minus_K_Z_minus_Z = find_H(q_R, q_Z-delta_Z, K_R, K_zeta, K_Z-delta_K_Z, launch_angular_frequency, mode_flag, 
                                 interp_poloidal_flux, interp_density_1D, interp_B_R, interp_B_T, interp_B_Z)
    d2H_dKZ_dZ = (H_plus_K_Z_plus_Z - H_plus_K_Z_minus_Z - H_minus_K_Z_plus_Z + H_minus_K_Z_minus_Z) / (4*delta_K_Z*delta_Z)  
    
    return d2H_dKZ_dZ


