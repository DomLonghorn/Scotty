# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 10:44:34 2018

@author: VH Hall-Chen
Valerian Hongjie Hall-Chen
valerian.chen@gmail.com
valerian.hall-chen.com

Run in Python 3, does not work in Python 2

Version history
v1 - Renaming of variables, cleaning up

Plan
Output everything to a file, and then do the analysis on that file.

1) Check that gradK_xi and such are done correctly, check that K_mag is calculated correctly when K_zeta is nonzero
2) Check that the calculation of Psi makes sense (and the rotation angle)
3) Check that K initial's calculation makes sense


Notes
- The loading of the input files was taken from integral_5 and modified.
- I should launch the beam inside the last closed flux surface
- K**2 = K_R**2 + K_z**2 + (K_zeta/r_R)**2, and K_zeta is constant (mode number). See 14 Sep 2018 notes.

Coordinates
X,Y,Z - Lab cartesian coordinates
R,zeta,Z - Lab cylindrical coordinates
x,y,g - Beam coordinates

Units
- SI units
- Distance in m
- Angles in rad
- electron cyclotron frequency positive
- K normalised such that K = 1 in vacuum. (Not implemented yet)
- Distance not normalised yet, should give it thought
- Start in vacuum, otherwise Psi_3D_beam_initial_cartersian does not get done properly
-

"""

import numpy as np
import math
from scipy import interpolate as interpolate
from scipy import integrate as integrate
from scipy import constants as constants
from scipy import linalg as linalg
import matplotlib.pyplot as plt
import os
from netCDF4 import Dataset
import bisect
import time

from Scotty_fun_general import read_floats_into_list_until, find_nearest, contract_special, find_x0, find_H, find_waist
from Scotty_fun_general import find_inverse_2D, find_Psi_3D_lab, find_q_lab_Cartesian, find_K_lab_Cartesian, find_K_lab, find_Psi_3D_lab_Cartesian
from Scotty_fun_general import find_normalised_plasma_freq, find_normalised_gyro_freq
from Scotty_fun_general import find_epsilon_para, find_epsilon_perp,find_epsilon_g
from Scotty_fun_general import find_dbhat_dR, find_dbhat_dZ
from Scotty_fun_general import find_d_poloidal_flux_dR, find_d_poloidal_flux_dZ,find_Psi_3D_plasma
from Scotty_fun_general import find_dB_dR_FFD, find_dB_dZ_FFD, find_d2B_dR2_FFD, find_d2B_dZ2_FFD, find_d2B_dR_dZ_FFD          
from Scotty_fun_general import find_dB_dR_CFD, find_dB_dZ_CFD, find_d2B_dR2_CFD, find_d2B_dZ2_CFD#, find_d2B_dR_dZ_FFD          
from Scotty_fun_general import find_d2_poloidal_flux_dR2, find_d2_poloidal_flux_dZ2
from Scotty_fun_general import find_H_Cardano, find_D

from Scotty_fun_evolution import ray_evolution_2D_fun, beam_evolution_fun


from Scotty_fun_FFD import find_dH_dR, find_dH_dZ # \nabla H
from Scotty_fun_CFD import find_dH_dKR, find_dH_dKZ, find_dH_dKzeta # \nabla_K H
from Scotty_fun_FFD import find_d2H_dR2, find_d2H_dZ2, find_d2H_dR_dZ # \nabla \nabla H
from Scotty_fun_CFD import find_d2H_dKR2, find_d2H_dKR_dKzeta, find_d2H_dKR_dKZ, find_d2H_dKzeta2, find_d2H_dKzeta_dKZ, find_d2H_dKZ2 # \nabla_K \nabla_K H
from Scotty_fun_mix import find_d2H_dKR_dR, find_d2H_dKR_dZ, find_d2H_dKzeta_dR, find_d2H_dKzeta_dZ, find_d2H_dKZ_dR, find_d2H_dKZ_dZ # \nabla_K \nabla H
from Scotty_fun_FFD import find_dpolflux_dR, find_dpolflux_dZ # For find_B if using efit files directly

def beam_me_up(poloidal_launch_angle_Torbeam,
               toroidal_launch_angle_Torbeam,
               launch_freq_GHz,
               mode_flag,
               vacuumLaunch_flag,
               launch_beam_width,
               launch_beam_radius_of_curvature,
               launch_position,
               find_B_method,
               efit_time_index=None,
               vacuum_propagation_flag=False,
               Psi_BC_flag = False,
               poloidal_flux_enter=None,
               input_filename_suffix='',
               output_filename_suffix='',
               figure_flag=True,
               plasmaLaunch_K=np.zeros(3),
               plasmaLaunch_Psi_3D_lab_Cartesian=np.zeros([3,3]),
               density_fit_parameters=None,
               acceleration_flag=False
               ):

    """
    find_B_method: 1) 'efit' finds B from efit files directly 2) 'torbeam' finds B from topfile
    """

    delta_R = -0.0001 #in the same units as data_R_coord
    delta_Z = 0.0001 #in the same units as data_Z_coord
    delta_K_R = 0.1 #in the same units as K_R
    delta_K_zeta = 0.1 #in the same units as K_zeta
    delta_K_Z = 0.1 #in the same units as K_z

    #major_radius = 0.9


    print('Beam trace me up, Scotty!')
    # ------------------------------




    # ------------------------------
     # Input data #
    # ------------------------------

    # # For analytical profile
    # # Specify the parameters
    #B_toroidal_max = 0.67 # in Tesla (?)
    #B_poloidal_max = 0.01 # in Tesla
    #
    #core_ne = 4.0 # units of 10^19 / m-3
    ##core_ne = 0
    #
    #aspect_ratio = 3.0 # major_radius/minor_radius
    #minor_radius = 0.5 # in meters
    #
    ## Calculates other parameters
    #major_radius = aspect_ratio * minor_radius
    #
    ## Generate ne
    #ne_data_length = 101
    #ne_data_fludata_R_coord = np.linspace(0,1,ne_data_length)
    #ne_data_density_array = np.zeros(ne_data_length)
    #ne_data_density_array = n_e_fun(ne_data_fludata_R_coord,core_ne)
    #
    #buffer_factor = 1.5
    #data_R_coord_length = 130
    #data_R_coord_start = major_radius - buffer_factor*minor_radius # in meters
    #data_R_coord_end = major_radius + buffer_factor*minor_radius
    #data_Z_coord_length = 65
    #data_Z_coord_start = -buffer_factor*minor_radius
    #data_Z_coord_end = buffer_factor*minor_radius
    #
    #data_R_coord = np.linspace(data_R_coord_start,data_R_coord_end,data_R_coord_length)
    #data_Z_coord = np.linspace(data_Z_coord_start,data_Z_coord_end,data_Z_coord_length)
    #
    #B_r = np.zeros([data_R_coord_length,data_Z_coord_length])
    #B_z = np.zeros([data_R_coord_length,data_Z_coord_length])
    #B_t = np.zeros([data_R_coord_length,data_Z_coord_length])
    #psi = np.zeros([data_R_coord_length,data_Z_coord_length])
    #
    #x_grid, z_grid = np.meshgrid(data_R_coord, data_Z_coord,indexing='ij')
    #data_B_T_grid = B_toroidal_fun(B_toroidal_max, data_R_coord, data_Z_coord, major_radius)
    #data_B_R_grid = B_r_fun(B_poloidal_max, data_R_coord, data_Z_coord, major_radius, minor_radius)
    #data_B_Z_grid = B_z_fun(B_poloidal_max, data_R_coord, data_Z_coord, major_radius, minor_radius)
    #data_poloidal_flux_grid = psi_fun(x_grid,z_grid, major_radius, minor_radius)
    ## ------------------------------

    # Tidying up the input data
    launch_angular_frequency = 2*math.pi*10.0**9 * launch_freq_GHz
    wavenumber_K0 = launch_angular_frequency / constants.c


    # Experimental Profile----------
    input_files_path ='D:\\Dropbox\\VHChen2020\\Data\\Input_Files_29Apr2019\\'
#    input_files_path ='D:\\Dropbox\\VHChen2018\\Data\\Input_Files_29Apr2019\\'
#    input_files_path ='D:\\Dropbox\\VHChen2019\\Code - Scotty\\Benchmark_9\\Torbeam\\'
#    input_files_path = os.path.dirname(os.path.abspath(__file__)) + '\\'
    
    if density_fit_parameters is None:
        print('ne(psi): loading from input file')
        print('Warning: Scale factor of 1.1 used')
        # Importing data from ne.dat
    #    ne_filename = input_files_path + 'ne' +input_filename_suffix+ '_smoothed.dat'
        ne_filename = input_files_path + 'ne' +input_filename_suffix+ '_fitted.dat'
    #    ne_filename = input_files_path + 'ne' +input_filename_suffix+ '.dat'
        
        ne_data = np.fromfile(ne_filename,dtype=float, sep='   ') # electron density as a function of poloidal flux label
        
    #    ne_data_length = int(ne_data[0])
        ne_data_density_array = 1.1*ne_data[2::2] # in units of 10.0**19 m-3
        ne_data_radialcoord_array = ne_data[1::2]
        ne_data_poloidal_flux_array = ne_data_radialcoord_array**2 # Loading radial coord for now, makes it easier to benchmark with Torbeam. Hence, have to convert to poloidal flux
    #    ne_data_normalised_poloidal_flux_array = ne_data[1::2] # My new IDL file outputs the flux density directly, instead of radialcoord        
        
        interp_density_1D = interpolate.interp1d(ne_data_poloidal_flux_array, ne_data_density_array,
                                                 kind='cubic', axis=-1, copy=True, bounds_error=False,
                                                 fill_value=0, assume_sorted=False) # density is 0 outside the LCFS, hence the fill_value. Use 'linear' instead of 'cubic' if the density data has a discontinuity in the first derivative    
        def find_density_1D(poloidal_flux, interp_density_1D=interp_density_1D):
            density = interp_density_1D(poloidal_flux)
            return density
    elif len(density_fit_parameters) == 4:
        print('ne(psi): Using order_1_polynomial*tanh')
        print('Warning: Scale factor of 1.1 used')
        ne_data_density_array=None # So that saving the input data later does not complain
        ne_data_radialcoord_array=None # So that saving the input data later does not complain
        
        def find_density_1D(poloidal_flux, poloidal_flux_enter=poloidal_flux_enter,density_fit_parameters=density_fit_parameters):
            density_fit = 1.1*(density_fit_parameters[0]*poloidal_flux + density_fit_parameters[1])*np.tanh(density_fit_parameters[2] * poloidal_flux + density_fit_parameters[3])
            is_inside = poloidal_flux <= poloidal_flux_enter # Boolean array
            density = is_inside * density_fit # The Boolean array sets stuff outside poloidal_flux_enter to zero
            return density     
    elif len(density_fit_parameters) == 3:
        print('ne(psi): using constant*tanh')
        ne_data_density_array=None # So that saving the input data later does not complain
        ne_data_radialcoord_array=None # So that saving the input data later does not complain
        
        def find_density_1D(poloidal_flux, poloidal_flux_enter=poloidal_flux_enter,density_fit_parameters=density_fit_parameters):
            density_fit = density_fit_parameters[0]*np.tanh( density_fit_parameters[1] * (poloidal_flux - density_fit_parameters[2]) )
            is_inside = poloidal_flux <= poloidal_flux_enter # Boolean array
            density = is_inside * density_fit # The Boolean array sets stuff outside poloidal_flux_enter to zero
            return density         
    else:
        print('density_fit_parameters has an invalid length')
    
    # This part of the code defines find_B_R, find_B_T, find_B_zeta
    interp_order = 5 # For the 2D interpolation functions
    interp_smoothing = 2 # For the 2D interpolation functions. For no smoothing, set to 0
    
    if find_B_method == 'torbeam':    
        print('Using Torbeam input files for B and poloidal flux')
        # topfile
        #Others: inbeam.dat, Te.dat (not currently used in this code)
    
        topfile_filename = input_files_path + 'topfile' +input_filename_suffix
    
        with open(topfile_filename) as f:
            while not 'X-coordinates' in f.readline(): pass # Start reading only from X-coords onwards
            data_R_coord = read_floats_into_list_until('Z-coordinates', f)
            data_Z_coord = read_floats_into_list_until('B_R', f)
            data_B_R_grid = read_floats_into_list_until('B_t', f)
            data_B_T_grid = read_floats_into_list_until('B_Z', f)
            data_B_Z_grid = read_floats_into_list_until('psi', f)
            data_poloidal_flux_grid = read_floats_into_list_until('you fall asleep', f)
        # ------------------------------


        data_B_R_grid = np.transpose((np.asarray(data_B_R_grid)).reshape(len(data_Z_coord),len(data_R_coord), order='C'))
        data_B_T_grid = np.transpose((np.asarray(data_B_T_grid)).reshape(len(data_Z_coord),len(data_R_coord), order='C'))
        data_B_Z_grid = np.transpose((np.asarray(data_B_Z_grid)).reshape(len(data_Z_coord),len(data_R_coord), order='C'))
        data_poloidal_flux_grid = np.transpose((np.asarray(data_poloidal_flux_grid)).reshape(len(data_Z_coord),len(data_R_coord), order='C'))
        # -------------------
    
        # Interpolation functions declared
        interp_B_R = interpolate.RectBivariateSpline(data_R_coord,data_Z_coord,data_B_R_grid, bbox=[None, None, None, None], kx=interp_order, ky=interp_order, s=interp_smoothing)
        interp_B_T = interpolate.RectBivariateSpline(data_R_coord,data_Z_coord,data_B_T_grid, bbox=[None, None, None, None], kx=interp_order, ky=interp_order, s=interp_smoothing)
        interp_B_Z = interpolate.RectBivariateSpline(data_R_coord,data_Z_coord,data_B_Z_grid, bbox=[None, None, None, None], kx=interp_order, ky=interp_order, s=interp_smoothing)
        
        def find_B_R(q_R,q_Z):
            B_R = interp_B_R(q_R,q_Z)
            return B_R
    
        def find_B_T(q_R,q_Z):
            B_T = interp_B_T(q_R,q_Z)
            return B_T
        
        def find_B_Z(q_R,q_Z):
            B_Z = interp_B_Z(q_R,q_Z)
            return B_Z

        interp_poloidal_flux = interpolate.RectBivariateSpline(data_R_coord,data_Z_coord,data_poloidal_flux_grid, bbox=[None, None, None, None], kx=interp_order, ky=interp_order, s=interp_smoothing)
        #    interp_poloidal_flux = interpolate.interp2d(data_R_coord, data_Z_coord, np.transpose(data_poloidal_flux_grid), kind='cubic',
        #                                           copy=True, bounds_error=False, fill_value=None) # Flux extrapolated outside region
        
        efit_time = None # To prevent the data-saving routines from complaining later on
        
    elif find_B_method == 'efit':
        print('Using efit output files directly for B and poloidal flux')

        dataset = Dataset('efitOut_29908.nc')
        
        time_index = efit_time_index
        time_array = dataset.variables['time'][:]
#        equilibriumStatus_array = dataset.variables['equilibriumStatus'][:]
        
        efit_time = time_array[time_index]
        print(efit_time)
        
        output_group = dataset.groups['output']
#        input_group = dataset.groups['input']
        
        profiles2D = output_group.groups['profiles2D']
        data_unnormalised_poloidal_flux_grid = profiles2D.variables['poloidalFlux'][time_index][:][:] #unnormalised, as a function of R and Z
        data_R_coord = profiles2D.variables['r'][time_index][:]
        data_Z_coord = profiles2D.variables['z'][time_index][:]
        
        # radialProfiles = output_group.groups['radialProfiles']
        # Bt_array = radialProfiles.variables['Bt'][time_index][:]
        # r_array_B = radialProfiles.variables['r'][time_index][:]
        
        # separatrixGeometry = output_group.groups['separatrixGeometry']
        # geometricAxis = separatrixGeometry.variables['geometricAxis'][time_index] # R,Z location of the geometric axis
        
        # globalParameters = output_group.groups['globalParameters']
        # bvacRgeom = globalParameters.variables['bvacRgeom'][time_index] # Vacuum B field (= B_zeta, in vacuum) at the geometric axis
        
        fluxFunctionProfiles = output_group.groups['fluxFunctionProfiles']
        poloidalFlux = fluxFunctionProfiles.variables['normalizedPoloidalFlux'][:]
        unnormalizedPoloidalFlux = fluxFunctionProfiles.variables['poloidalFlux'][time_index][:] # poloidalFlux as a function of normalised poloidal flux
        rBphi = fluxFunctionProfiles.variables['rBphi'][time_index][:] 

        dataset.close()

        interp_rBphi = interpolate.interp1d(poloidalFlux, rBphi, 
                                                 kind='cubic', axis=-1, copy=True, bounds_error=False, 
                                                 fill_value=rBphi[-1], assume_sorted=False)
 
        # normalised_polflux = polflux_const_m * poloidalFlux + polflux_const_c (think y = mx + c)        
        [polflux_const_m, polflux_const_c] = np.polyfit(unnormalizedPoloidalFlux,poloidalFlux,1) #linear fit
        data_poloidal_flux_grid = data_unnormalised_poloidal_flux_grid*polflux_const_m + polflux_const_c
        
        interp_poloidal_flux = interpolate.RectBivariateSpline(data_R_coord,data_Z_coord,data_poloidal_flux_grid, bbox=[None, None, None, None], kx=interp_order, ky=interp_order, s=interp_smoothing)

            # Extrapolation assumes the last element of rBphi corresponds to poloidalflux = 1
        rBphi_gradient_for_extrapolation = (rBphi[-1] - rBphi[-2])/(poloidalFlux[-1] - poloidalFlux[-2])
        last_poloidal_flux = poloidalFlux[-1]
        last_rBphi = rBphi[-1]
        
        def find_B_R(q_R,q_Z,delta_R=delta_R,interp_poloidal_flux=interp_poloidal_flux,polflux_const_m=polflux_const_m):
            dpolflux_dZ = find_dpolflux_dZ(q_R,q_Z,delta_R,interp_poloidal_flux)
            B_R = -  dpolflux_dZ / (polflux_const_m * q_R)
            return B_R
        
        def find_B_T(q_R,q_Z,
                     last_poloidal_flux=last_poloidal_flux,
                     last_rBphi=last_rBphi,
                     rBphi_gradient_for_extrapolation=rBphi_gradient_for_extrapolation,
                     interp_poloidal_flux=interp_poloidal_flux,
                     interp_rBphi=interp_rBphi):
            # B_T from EFIT

            polflux = interp_poloidal_flux(q_R,q_Z, grid=False)

            is_inside = polflux <= last_poloidal_flux # Boolean array. True if it's inside the range where the data is available
            is_outside = ~is_inside
            
            rBphi = ( is_inside *   interp_rBphi(polflux)
                    + is_outside * (rBphi_gradient_for_extrapolation*(polflux-last_poloidal_flux) + last_rBphi)
                    )
                
            B_T = rBphi/q_R

            return B_T
        
        def find_B_Z(q_R,q_Z,delta_Z=delta_Z,interp_poloidal_flux=interp_poloidal_flux,polflux_const_m=polflux_const_m):
            dpolflux_dR = find_dpolflux_dR(q_R,q_Z,delta_Z,interp_poloidal_flux)
            B_Z = dpolflux_dR / (polflux_const_m * q_R)
            return B_Z


    else:
        print('Invalid find_B_method')
        














    ## -------------------
    ## Launch parameters    
    ## -------------------
    if vacuumLaunch_flag:
        print('Beam launched from outside the plasma')
        
        K_R_launch    = -wavenumber_K0 * np.cos( toroidal_launch_angle_Torbeam/180.0*math.pi ) * np.cos( poloidal_launch_angle_Torbeam/180.0*math.pi ) # K_R
        K_zeta_launch = -wavenumber_K0 * np.sin( toroidal_launch_angle_Torbeam/180.0*math.pi ) * np.cos( poloidal_launch_angle_Torbeam/180.0*math.pi ) * launch_position[0]# K_zeta
        K_Z_launch    = -wavenumber_K0 * np.sin( poloidal_launch_angle_Torbeam/180.0*math.pi ) # K_z    
        launch_K = np.array([K_R_launch,K_zeta_launch,K_Z_launch])
        
        poloidal_launch_angle_Rz = (180.0+poloidal_launch_angle_Torbeam)/180.0*np.pi
        poloidal_rotation_angle = (90.0+poloidal_launch_angle_Torbeam)/180.0*np.pi
        toroidal_launch_angle_Rz = (180.0+toroidal_launch_angle_Torbeam)/180.0*np.pi 
    
        Psi_w_beam_launch_cartersian = np.array(
                [
                [ wavenumber_K0/launch_beam_radius_of_curvature+2j*launch_beam_width**(-2), 0],
                [ 0, wavenumber_K0/launch_beam_radius_of_curvature+2j*launch_beam_width**(-2)]
                ]
                )    
        identity_matrix_2D = np.array(
                [
                [ 1, 0],
                [ 0, 1]
                ]
                )    
        
        rotation_matrix_pol = np.array( [
                [ np.cos(poloidal_rotation_angle), 0, np.sin(poloidal_rotation_angle) ],
                [ 0, 1, 0 ],
                [ -np.sin(poloidal_rotation_angle), 0, np.cos(poloidal_rotation_angle) ]
                ] )
        
        rotation_matrix_tor = np.array( [
                [ np.cos(toroidal_launch_angle_Torbeam/180.0*math.pi), np.sin(toroidal_launch_angle_Torbeam/180.0*math.pi), 0 ],
                [ -np.sin(toroidal_launch_angle_Torbeam/180.0*math.pi), np.cos(toroidal_launch_angle_Torbeam/180.0*math.pi), 0 ],
                [ 0,0,1 ]
                ] )
        
        rotation_matrix         = np.matmul(rotation_matrix_pol,rotation_matrix_tor)
        rotation_matrix_inverse = np.transpose(rotation_matrix)
        
        Psi_3D_beam_launch_cartersian = np.array([
                [ Psi_w_beam_launch_cartersian[0][0], Psi_w_beam_launch_cartersian[0][1], 0 ],
                [ Psi_w_beam_launch_cartersian[1][0], Psi_w_beam_launch_cartersian[1][1], 0 ],
                [ 0, 0, 0 ]
                ])
        Psi_3D_lab_launch_cartersian = np.matmul( rotation_matrix_inverse, np.matmul(Psi_3D_beam_launch_cartersian, rotation_matrix) )
        Psi_3D_lab_launch = find_Psi_3D_lab(Psi_3D_lab_launch_cartersian,launch_position[0],launch_position[1],K_R_launch,K_zeta_launch)
            
        if vacuum_propagation_flag:
            Psi_w_beam_inverse_launch_cartersian = find_inverse_2D(Psi_w_beam_launch_cartersian)
           
            # Finds entry point
            search_Z_end = launch_position[2] - launch_position[0]*np.tan(np.radians(poloidal_launch_angle_Torbeam))
            numberOfCoarseSearchPoints = 50
            R_coarse_search_array = np.linspace(launch_position[0],0,numberOfCoarseSearchPoints)
            Z_coarse_search_array = np.linspace(launch_position[2],search_Z_end,numberOfCoarseSearchPoints)
            poloidal_flux_coarse_search_array = np.zeros(numberOfCoarseSearchPoints)
            poloidal_flux_coarse_search_array = interp_poloidal_flux(R_coarse_search_array,Z_coarse_search_array,grid=False)
            meets_flux_condition_array = poloidal_flux_coarse_search_array < 0.9*poloidal_flux_enter
            dummy_array = np.array(range(numberOfCoarseSearchPoints))
            indices_inside_for_sure_array = dummy_array[meets_flux_condition_array]
            first_inside_index = indices_inside_for_sure_array[0]
            numberOfFineSearchPoints = 7500
            R_fine_search_array = np.linspace(launch_position[0],R_coarse_search_array[first_inside_index],numberOfFineSearchPoints)
            Z_fine_search_array = np.linspace(launch_position[2],Z_coarse_search_array[first_inside_index],numberOfFineSearchPoints)
            poloidal_fine_search_array = np.zeros(numberOfFineSearchPoints)
            poloidal_fine_search_array = interp_poloidal_flux(R_fine_search_array,Z_fine_search_array,grid=False)
            entry_index = find_nearest(poloidal_fine_search_array,poloidal_flux_enter)
            if poloidal_fine_search_array[entry_index] > poloidal_flux_enter:
                # The first point needs to be in the plasma
                # If the first point is outside, then there will be errors when the gradients are calculated
                entry_index = entry_index + 1
            entry_position = np.zeros(3) # R,Z
            entry_position[0] = R_fine_search_array[entry_index]
            entry_position[1] = K_zeta_launch/K_R_launch * ( 1/launch_position[0] - 1/entry_position[0] )
            entry_position[2] = Z_fine_search_array[entry_index]

            distance_from_launch_to_entry = np.sqrt(
                                                    launch_position[0]**2 
                                                    + entry_position[0]**2 
                                                    - 2 * launch_position[0] * entry_position[0] * np.cos(entry_position[1] - launch_position[1])
                                                    + (launch_position[2] - entry_position[2])**2
                                                    )
            # Entry point found

            
            # Calculate entry parameters from launch parameters
            # That is, find beam at start of plasma given its parameters at the antenna
            K_lab_launch = np.array([K_R_launch,K_zeta_launch,K_Z_launch])
            K_lab_Cartesian_launch = find_K_lab_Cartesian(K_lab_launch,launch_position)
            K_lab_Cartesian_entry = K_lab_Cartesian_launch
            entry_position_Cartesian = find_q_lab_Cartesian(entry_position)
            K_lab_entry = find_K_lab(K_lab_Cartesian_entry,entry_position_Cartesian)
            
            K_R_entry    = K_lab_entry[0] # K_R
            K_zeta_entry = K_lab_entry[1]
            K_Z_entry    = K_lab_entry[2] # K_z
        
            Psi_w_beam_inverse_entry_cartersian = distance_from_launch_to_entry/(wavenumber_K0)*identity_matrix_2D + Psi_w_beam_inverse_launch_cartersian
            Psi_w_beam_entry_cartersian = find_inverse_2D(Psi_w_beam_inverse_entry_cartersian)
        
            Psi_3D_beam_entry_cartersian = np.array([
                    [ Psi_w_beam_entry_cartersian[0][0], Psi_w_beam_entry_cartersian[0][1], 0 ],
                    [ Psi_w_beam_entry_cartersian[1][0], Psi_w_beam_entry_cartersian[1][1], 0 ],
                    [ 0, 0, 0 ]
                    ]) # 'entry' is still in vacuum, so the components of Psi along g are all 0 (since \nabla H = 0)
            
        
            Psi_3D_lab_entry_cartersian = np.matmul( rotation_matrix_inverse, np.matmul(Psi_3D_beam_entry_cartersian, rotation_matrix) )
        
        
            # Convert to cylindrical coordinates
            Psi_3D_lab_entry = find_Psi_3D_lab(Psi_3D_lab_entry_cartersian,entry_position[0],entry_position[1],K_R_entry,K_zeta_entry)

            # -------------------
            # Find initial parameters in plasma
            # -------------------
            K_R_initial        = K_R_entry
            K_zeta_initial     = K_zeta_entry
            K_Z_initial        = K_Z_entry
            initial_position   = entry_position
            if Psi_BC_flag: # Use BCs
                dH_dKR_initial    = find_dH_dKR(initial_position[0], initial_position[2], K_R_initial, K_zeta_initial, K_Z_initial,
                                             launch_angular_frequency, mode_flag, delta_K_R, 
                                             interp_poloidal_flux, find_density_1D, find_B_R, find_B_T, find_B_Z)
                dH_dKzeta_initial = find_dH_dKzeta(initial_position[0], initial_position[2], K_R_initial, K_zeta_initial, K_Z_initial,
                                             launch_angular_frequency, mode_flag, delta_K_zeta, 
                                             interp_poloidal_flux, find_density_1D, find_B_R, find_B_T, find_B_Z) 
                dH_dKZ_initial    = find_dH_dKZ(initial_position[0], initial_position[2], K_R_initial, K_zeta_initial, K_Z_initial,
                                             launch_angular_frequency, mode_flag, delta_K_Z, 
                                             interp_poloidal_flux, find_density_1D, find_B_R, find_B_T, find_B_Z)
                dH_dR_initial     = find_dH_dR(initial_position[0], initial_position[2], K_R_initial, K_zeta_initial, K_Z_initial, 
                                            launch_angular_frequency, mode_flag, delta_R, 
                                            interp_poloidal_flux, find_density_1D, find_B_R, find_B_T, find_B_Z)
                dH_dZ_initial     = find_dH_dZ(initial_position[0], initial_position[2], K_R_initial, K_zeta_initial, K_Z_initial, 
                                            launch_angular_frequency, mode_flag, delta_Z, 
                                            interp_poloidal_flux, find_density_1D, find_B_R, find_B_T, find_B_Z)
                d_poloidal_flux_d_R_boundary = find_d_poloidal_flux_dR(initial_position[0], initial_position[2], delta_R, interp_poloidal_flux)
                d_poloidal_flux_d_Z_boundary = find_d_poloidal_flux_dZ(initial_position[0], initial_position[2], delta_R, interp_poloidal_flux)

                Psi_3D_lab_initial = find_Psi_3D_plasma(Psi_3D_lab_entry,
                                                        dH_dKR_initial, dH_dKzeta_initial, dH_dKZ_initial,
                                                        dH_dR_initial, dH_dZ_initial,
                                                        d_poloidal_flux_d_R_boundary.item(), d_poloidal_flux_d_Z_boundary.item()) #. item() to convert variable from type ndarray to float, such that the array elements all have the same type
            else: # Do not use BCs
                Psi_3D_lab_initial = Psi_3D_lab_entry
            
        else: #Run solver from the launch position, no analytical vacuum propagation
            Psi_3D_lab_initial = Psi_3D_lab_launch
            K_R_initial        = K_R_launch
            K_zeta_initial     = K_zeta_launch
            K_Z_initial        = K_Z_launch
            initial_position   = launch_position            
            
            Psi_3D_lab_entry = None
            distance_from_launch_to_entry=None
            Psi_3D_lab_entry_cartersian = np.full_like(Psi_3D_lab_launch,fill_value=np.nan)
            

    else:
        print('Beam launched from inside the plasma')
        K_R_launch = plasmaLaunch_K[0]
        K_zeta_launch = plasmaLaunch_K[1]
        K_Z_launch = plasmaLaunch_K[2]
        
        Psi_3D_lab_initial = find_Psi_3D_lab(plasmaLaunch_Psi_3D_lab_Cartesian,launch_position[0],plasmaLaunch_K[0],plasmaLaunch_K[1])
        K_R_initial        = K_R_launch
        K_zeta_initial     = K_zeta_launch
        K_Z_initial        = K_Z_launch
        initial_position   = launch_position
        
#        print(K_R_initial)
#        print(K_zeta_launch)
#        print(K_Z_launch)
    # -------------------



    # -------------------

    # Initial conditions for the solver
    beam_parameters_initial = np.zeros(17)
    # This used to be complex, with a length of 11, but the solver throws a warning saying that something is casted to real
    # It seems to be fine, bu
    
    beam_parameters_initial[0]  = initial_position[0] # q_R
    beam_parameters_initial[1]  = initial_position[1] # q_zeta (This should be = 0)
    beam_parameters_initial[2]  = initial_position[2] # q_Z
    
    beam_parameters_initial[3]  = K_R_initial
    beam_parameters_initial[4]  = K_Z_initial
    
    beam_parameters_initial[5]  = np.real(Psi_3D_lab_initial[0,0]) # Psi_RR
    beam_parameters_initial[6]  = np.real(Psi_3D_lab_initial[1,1]) # Psi_zetazeta 
    beam_parameters_initial[7]  = np.real(Psi_3D_lab_initial[2,2]) # Psi_ZZ
    beam_parameters_initial[8]  = np.real(Psi_3D_lab_initial[0,1]) # Psi_Rzeta 
    beam_parameters_initial[9]  = np.real(Psi_3D_lab_initial[0,2]) # Psi_RZ 
    beam_parameters_initial[10] = np.real(Psi_3D_lab_initial[1,2]) # Psi_zetaZ 
    
    beam_parameters_initial[11] = np.imag(Psi_3D_lab_initial[0,0]) # Psi_RR
    beam_parameters_initial[12] = np.imag(Psi_3D_lab_initial[1,1]) # Psi_zetazeta 
    beam_parameters_initial[13] = np.imag(Psi_3D_lab_initial[2,2]) # Psi_ZZ
    beam_parameters_initial[14] = np.imag(Psi_3D_lab_initial[0,1]) # Psi_Rzeta 
    beam_parameters_initial[15] = np.imag(Psi_3D_lab_initial[0,2]) # Psi_RZ 
    beam_parameters_initial[16] = np.imag(Psi_3D_lab_initial[1,2]) # Psi_zetaZ  
    
    ray_parameters_initial    = np.real(beam_parameters_initial[0:5])
    ray_parameters_2D_initial = np.delete(ray_parameters_initial,1) # Remove q_zeta
    # -------------------

    # Define events for the solver
    def event_leave_plasma(tau, beam_parameters, K_zeta, 
                           launch_angular_frequency, mode_flag,
                           delta_R, delta_Z, delta_K_R, delta_K_zeta, delta_K_Z,
                           interp_poloidal_flux, find_density_1D, 
                           find_B_R, find_B_T, find_B_Z,
                           poloidal_flux_enter=poloidal_flux_enter):
        
        q_R          = beam_parameters[0]
        q_Z          = beam_parameters[1]
        
        poloidal_flux = interp_poloidal_flux(q_R,q_Z)
        
        poloidal_flux_difference = poloidal_flux - poloidal_flux_enter
        # goes from negative to positive when leaving the plasma
        return poloidal_flux_difference
    event_leave_plasma.terminal = True # Stop the solver when the beam leaves the plasma
    event_leave_plasma.direction = 1.0 # positive value, when function result goes from negative to positive
    
    event_leave_plasma.terminal = True # Stop the solver when the beam leaves the plasma
    event_leave_plasma.direction = 1.0 # positive value, when function result goes from negative to positive
    # -------------------


    # Propagate the beam
        # Calls scipy's initial value problem solver   
        
    tau_max = 10**5 # If the ray hasn't left the plasma by the time this tau is reached, the solver gives up
    solver_arguments = (K_zeta_initial, 
                        launch_angular_frequency, mode_flag,
                        delta_R, delta_Z, delta_K_R, delta_K_zeta, delta_K_Z,
                        interp_poloidal_flux, find_density_1D, 
                        find_B_R, find_B_T, find_B_Z) # Stuff the solver needs to evolve beam_parameters

    """
    Propagate a ray. Quickly finds tau at which the ray leaves the plasma, as well as estimates location of cut-off
    """
    solver_start_time = time.time()
    
    solver_ray_output = integrate.solve_ivp(
                        ray_evolution_2D_fun, [0,tau_max], ray_parameters_2D_initial, 
                        method='RK45', t_eval=None, dense_output=(acceleration_flag), 
                        events=event_leave_plasma, vectorized=False, args=solver_arguments
                    ) # Compute dense output if acceleration_flag, do not compute if not


    solver_end_time = time.time()
    print('Time taken (ray solver)', solver_end_time-solver_start_time,'s')


    solver_ray_status = solver_ray_output.status
    if solver_ray_status != 1:
        print('Warning: Ray has not left plasma. Increase tau_max')
    
    ray_parameters_2D = solver_ray_output.y    
    K_magnitude_ray = (  (np.real(ray_parameters_2D[2,:]))**2 
                        + K_zeta_initial**2/(np.real(ray_parameters_2D[0,:]))**2 
                        + (np.real(ray_parameters_2D[3,:]))**2)**(0.5)
    
    tau_ray = solver_ray_output.t
    index_cutoff_estimate = K_magnitude_ray.argmin()

    tau_leave = np.squeeze(solver_ray_output.t_events)

        ##
    
    



        ##    
    
    if acceleration_flag:
        """
        - Propagates a ray with the appropriate points
        - Calculate second derivatives along the ray
        - Get beam properties from those derivatives
        """   
        numberOfTauPointsFine = 1000 + 1
        tau_points_fine = np.linspace(0,0.995*tau_leave,numberOfTauPointsFine)
        # # factor of 0.995 so that the last point would definitely still be in the plasma
        
        ray_parameters_2D_fine = solver_ray_output.sol(tau_points_fine)
        
        [q_R_fine, q_Z_fine, K_R_fine, K_Z_fine] = ray_parameters_2D_fine
        
   
        # \nabla \nabla H
        d2H_dR2_fine   = find_d2H_dR2(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                 launch_angular_frequency, mode_flag,
                                 delta_R,
                                 interp_poloidal_flux, find_density_1D,
                                 find_B_R, find_B_T, find_B_Z
                                 )
        d2H_dR_dZ_fine = find_d2H_dR_dZ(q_R_fine ,q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                   launch_angular_frequency, mode_flag,
                                   delta_R, delta_Z,
                                   interp_poloidal_flux, find_density_1D,
                                   find_B_R, find_B_T, find_B_Z
                                   )
        d2H_dZ2_fine   = find_d2H_dZ2(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                 launch_angular_frequency, mode_flag,
                                 delta_Z,
                                 interp_poloidal_flux, find_density_1D,
                                 find_B_R, find_B_T, find_B_Z
                                 )
        # grad_grad_H = np.squeeze(np.array([
        #     [d2H_dR2.item()  , 0.0, d2H_dR_dZ.item()],
        #     [0.0             , 0.0, 0.0             ],
        #     [d2H_dR_dZ.item(), 0.0, d2H_dZ2.item()  ] 
        #     ]))
        grad_grad_H_fine        = np.zeros([numberOfTauPointsFine,3,3])
        grad_grad_H_fine[:,0,0] = d2H_dR2_fine
        grad_grad_H_fine[:,2,2] = d2H_dZ2_fine
        grad_grad_H_fine[:,2,0] = d2H_dR_dZ_fine
        grad_grad_H_fine[:,0,2] = grad_grad_H_fine[:,2,0]
    
        # \nabla_K \nabla H
        d2H_dKR_dR_fine    = find_d2H_dKR_dR(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                        launch_angular_frequency, mode_flag,
                                        delta_K_R, delta_R,
                                        interp_poloidal_flux, find_density_1D,
                                        find_B_R, find_B_T, find_B_Z
                                        )
        d2H_dKZ_dZ_fine    = find_d2H_dKZ_dZ(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                        launch_angular_frequency, mode_flag,
                                        delta_K_Z, delta_Z,
                                        interp_poloidal_flux, find_density_1D,
                                        find_B_R, find_B_T, find_B_Z
                                        )
        d2H_dKR_dZ_fine    = find_d2H_dKR_dZ(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                        launch_angular_frequency, mode_flag,
                                        delta_K_R, delta_Z,
                                        interp_poloidal_flux, find_density_1D,
                                        find_B_R, find_B_T, find_B_Z
                                        )
        d2H_dKzeta_dZ_fine = find_d2H_dKzeta_dZ(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine, 
                                           launch_angular_frequency, mode_flag,
                                           delta_K_zeta, delta_Z,
                                           interp_poloidal_flux, find_density_1D,
                                           find_B_R, find_B_T, find_B_Z
                                           )
        d2H_dKzeta_dR_fine = find_d2H_dKzeta_dR(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                           launch_angular_frequency, mode_flag,
                                           delta_K_zeta, delta_R,
                                           interp_poloidal_flux, find_density_1D,
                                           find_B_R, find_B_T, find_B_Z
                                           )
        d2H_dKZ_dR_fine    = find_d2H_dKZ_dR(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                        launch_angular_frequency, mode_flag,
                                        delta_K_Z, delta_R,
                                        interp_poloidal_flux, find_density_1D,
                                        find_B_R, find_B_T, find_B_Z
                                        )
        # gradK_grad_H = np.squeeze(np.array([
        #     [d2H_dKR_dR.item(),    0.0, d2H_dKR_dZ.item()   ],
        #     [d2H_dKzeta_dR.item(), 0.0, d2H_dKzeta_dZ.item()],
        #     [d2H_dKZ_dR.item(),    0.0, d2H_dKZ_dZ.item()   ]
        #     ]))
        gradK_grad_H_fine        = np.zeros([numberOfTauPointsFine,3,3])
        gradK_grad_H_fine[:,0,0] = d2H_dKR_dR_fine
        gradK_grad_H_fine[:,2,2] = d2H_dKZ_dZ_fine
        gradK_grad_H_fine[:,0,2] = d2H_dKR_dZ_fine
        gradK_grad_H_fine[:,1,0] = d2H_dKzeta_dR_fine
        gradK_grad_H_fine[:,1,2] = d2H_dKzeta_dZ_fine
        gradK_grad_H_fine[:,2,0] = d2H_dKZ_dR_fine
        grad_gradK_H_fine        = np.transpose(gradK_grad_H_fine,axes=[0,2,1])

        # \nabla_K \nabla_K H
        d2H_dKR2_fine       = find_d2H_dKR2(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                       launch_angular_frequency, mode_flag,
                                       delta_K_R,
                                       interp_poloidal_flux, find_density_1D,
                                       find_B_R, find_B_T, find_B_Z
                                       )
        d2H_dKzeta2_fine    = find_d2H_dKzeta2(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                          launch_angular_frequency, mode_flag,
                                          delta_K_zeta,
                                          interp_poloidal_flux, find_density_1D,
                                          find_B_R, find_B_T, find_B_Z
                                          )
        d2H_dKZ2_fine       = find_d2H_dKZ2(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                       launch_angular_frequency, mode_flag,
                                       delta_K_Z,
                                       interp_poloidal_flux, find_density_1D,
                                       find_B_R, find_B_T, find_B_Z
                                       )
        d2H_dKR_dKzeta_fine = find_d2H_dKR_dKzeta(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                             launch_angular_frequency, mode_flag,
                                             delta_K_R, delta_K_zeta,
                                             interp_poloidal_flux, find_density_1D,
                                             find_B_R, find_B_T, find_B_Z
                                             )
        d2H_dKR_dKZ_fine    = find_d2H_dKR_dKZ(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                          launch_angular_frequency, mode_flag,
                                          delta_K_R, delta_K_Z,
                                          interp_poloidal_flux, find_density_1D,
                                          find_B_R, find_B_T, find_B_Z
                                          )
        d2H_dKzeta_dKZ_fine = find_d2H_dKzeta_dKZ(q_R_fine, q_Z_fine, K_R_fine, K_zeta_initial, K_Z_fine,
                                             launch_angular_frequency, mode_flag,
                                             delta_K_zeta, delta_K_Z,
                                             interp_poloidal_flux, find_density_1D,
                                             find_B_R, find_B_T, find_B_Z
                                             )
        # gradK_gradK_H = np.squeeze(np.array([
        #     [d2H_dKR2.item()      , d2H_dKR_dKzeta.item(), d2H_dKR_dKZ.item()   ],
        #     [d2H_dKR_dKzeta.item(), d2H_dKzeta2.item()   , d2H_dKzeta_dKZ.item()],
        #     [d2H_dKR_dKZ.item()   , d2H_dKzeta_dKZ.item(), d2H_dKZ2.item()      ]
        #     ]))
        gradK_gradK_H_fine        = np.zeros([numberOfTauPointsFine,3,3])
        gradK_gradK_H_fine[:,0,0] = d2H_dKR2_fine
        gradK_gradK_H_fine[:,1,1] = d2H_dKzeta2_fine
        gradK_gradK_H_fine[:,2,2] = d2H_dKZ2_fine
        gradK_gradK_H_fine[:,0,1] = d2H_dKR_dKzeta_fine
        gradK_gradK_H_fine[:,0,2] = d2H_dKR_dKZ_fine
        gradK_gradK_H_fine[:,1,2] = d2H_dKzeta_dKZ_fine
        gradK_gradK_H_fine[:,1,0] = gradK_grad_H_fine[:,0,1]
        gradK_gradK_H_fine[:,2,0] = gradK_grad_H_fine[:,0,2]
        gradK_gradK_H_fine[:,2,1] = gradK_grad_H_fine[:,1,2]

        Psi_3D_fine        = np.zeros([numberOfTauPointsFine,3,3], dtype='complex128')

        Psi_3D_fine[0,0,0] = beam_parameters_initial[5]  + 1j*beam_parameters_initial[11] # d (Psi_RR) / d tau
        Psi_3D_fine[0,1,1] = beam_parameters_initial[6]  + 1j*beam_parameters_initial[12] # d (Psi_zetazeta) / d tau
        Psi_3D_fine[0,2,2] = beam_parameters_initial[7]  + 1j*beam_parameters_initial[13] # d (Psi_ZZ) / d tau
        Psi_3D_fine[0,0,1] = beam_parameters_initial[8]  + 1j*beam_parameters_initial[14] # d (Psi_Rzeta) / d tau
        Psi_3D_fine[0,0,2] = beam_parameters_initial[9]  + 1j*beam_parameters_initial[15] # d (Psi_RZ) / d tau
        Psi_3D_fine[0,1,2] = beam_parameters_initial[10] + 1j*beam_parameters_initial[16] # d (Psi_zetaZ) / d tau
        Psi_3D_fine[0,1,0] = Psi_3D_fine[0,0,1]
        Psi_3D_fine[0,2,0] = Psi_3D_fine[0,0,2]
        Psi_3D_fine[0,2,1] = Psi_3D_fine[0,1,2]

        for ii in range(0,numberOfTauPointsFine-1):
            print(ii)
            Psi_3D_fine[ii+1,:,:]= (( - grad_grad_H_fine[ii,:,:].squeeze()
                    - np.matmul(
                                Psi_3D_fine[ii,:,:].squeeze(), gradK_grad_H_fine[ii,:,:].squeeze()
                                )
                    - np.matmul(
                                grad_gradK_H_fine[ii,:,:].squeeze(), Psi_3D_fine[ii,:,:].squeeze()
                                )
                    - np.matmul(np.matmul(
                                Psi_3D_fine[ii,:,:].squeeze(), gradK_gradK_H_fine[ii,:,:].squeeze()
                                ),
                                Psi_3D_fine[ii,:,:].squeeze()
                                )
                    ) * (tau_points_fine[ii+1] - tau_points_fine[ii])
                    +
                    Psi_3D_fine[ii,:,:])
            
        print(Psi_3D_fine[0,:,:])
        print(Psi_3D_fine[1,:,:])
        print(Psi_3D_fine[2,:,:])
            
        plt.figure()
        plt.subplot(3,3,1)
        plt.plot(tau_points_fine,np.real(Psi_3D_fine[:,0,0]))
        plt.subplot(3,3,2)
        plt.plot(tau_points_fine,np.real(Psi_3D_fine[:,0,1]))        
        plt.subplot(3,3,3)
        plt.plot(tau_points_fine,np.real(Psi_3D_fine[:,0,2]))        
        plt.subplot(3,3,5)
        plt.plot(tau_points_fine,np.real(Psi_3D_fine[:,1,1]))        
        plt.subplot(3,3,6)
        plt.plot(tau_points_fine,np.real(Psi_3D_fine[:,1,2]))        
        plt.subplot(3,3,9)
        plt.plot(tau_points_fine,np.real(Psi_3D_fine[:,2,2]))

        plt.figure()
        plt.subplot(3,3,1)
        plt.plot(tau_points_fine,np.imag(Psi_3D_fine[:,0,0]))
        plt.subplot(3,3,2)
        plt.plot(tau_points_fine,np.imag(Psi_3D_fine[:,0,1]))        
        plt.subplot(3,3,3)
        plt.plot(tau_points_fine,np.imag(Psi_3D_fine[:,0,2]))        
        plt.subplot(3,3,5)
        plt.plot(tau_points_fine,np.imag(Psi_3D_fine[:,1,1]))        
        plt.subplot(3,3,6)
        plt.plot(tau_points_fine,np.imag(Psi_3D_fine[:,1,2]))        
        plt.subplot(3,3,9)
        plt.plot(tau_points_fine,np.imag(Psi_3D_fine[:,2,2]))
                
        plt.figure()
        plt.subplot(3,3,1)
        plt.plot(tau_points_fine,gradK_grad_H_fine[:,0,0],'r')
        plt.subplot(3,3,2)
        plt.subplot(3,3,3)
        plt.plot(tau_points_fine,gradK_grad_H_fine[:,0,2],'r')
        plt.subplot(3,3,4)
        plt.plot(tau_points_fine,gradK_grad_H_fine[:,1,0],'r')
        plt.subplot(3,3,5)
        plt.subplot(3,3,6)
        plt.plot(tau_points_fine,gradK_grad_H_fine[:,1,2],'r')
        plt.subplot(3,3,7)
        plt.plot(tau_points_fine,gradK_grad_H_fine[:,2,0],'r')
        plt.subplot(3,3,8)
        plt.subplot(3,3,9)
        plt.title('d2H_dKZ_dZ')
        plt.plot(tau_points_fine,gradK_grad_H_fine[:,2,2],'r')
        
        plt.figure()
        plt.subplot(3,3,1)
        plt.plot(tau_points_fine,grad_grad_H_fine[:,0,0],'or')
        plt.title('d2H_dR2')
        plt.subplot(3,3,2)
        plt.subplot(3,3,3)
        plt.plot(tau_points_fine,grad_grad_H_fine[:,0,2],'r')
        plt.title('d2H_dR_dZ')
        plt.subplot(3,3,4)
        plt.subplot(3,3,5)
        plt.subplot(3,3,6)
        plt.subplot(3,3,7)
        #plt.plot(tau_points_fine,grad_grad_H_fine[:,2,0],'r')
        plt.subplot(3,3,8)
        plt.subplot(3,3,9)
        plt.plot(tau_points_fine,grad_grad_H_fine[:,2,2],'r')
        plt.title('d2H_dZ2')
        
        
        plt.figure()
        plt.subplot(3,3,1)
        plt.plot(tau_points_fine,gradK_gradK_H_fine[:,0,0],'r')
        plt.subplot(3,3,2)
        plt.plot(tau_points_fine,gradK_gradK_H_fine[:,0,1],'r')
        plt.subplot(3,3,3)
        plt.plot(tau_points_fine,gradK_gradK_H_fine[:,0,2],'r')
        plt.subplot(3,3,4)
        plt.subplot(3,3,5)
        plt.plot(tau_points_fine,gradK_gradK_H_fine[:,1,1],'r')
        plt.subplot(3,3,6)
        plt.plot(tau_points_fine,gradK_gradK_H_fine[:,1,2],'r')
        plt.subplot(3,3,7)
        plt.subplot(3,3,8)
        plt.subplot(3,3,9)
        plt.title('d2H_dKZ2')
        plt.plot(tau_points_fine,gradK_gradK_H_fine[:,2,2],'r')  
        
        # numberOfDataPoints = len(tau_array)
        
        # q_R_array = np.real(beam_parameters[0,:])
        # q_zeta_array = np.real(beam_parameters[1,:])
        # q_Z_array = np.real(beam_parameters[2,:])
        
        # K_R_array = np.real(beam_parameters[3,:])
        # K_Z_array = np.real(beam_parameters[4,:])
        
        # Psi_3D_output = np.zeros([numberOfDataPoints,3,3],dtype='complex128')
        # Psi_3D_output[:,0,0] = beam_parameters[5,:]  + 1j*beam_parameters[11,:] # d (Psi_RR) / d tau
        # Psi_3D_output[:,1,1] = beam_parameters[6,:]  + 1j*beam_parameters[12,:] # d (Psi_zetazeta) / d tau
        # Psi_3D_output[:,2,2] = beam_parameters[7,:]  + 1j*beam_parameters[13,:] # d (Psi_ZZ) / d tau
        # Psi_3D_output[:,0,1] = beam_parameters[8,:]  + 1j*beam_parameters[14,:] # d (Psi_Rzeta) / d tau
        # Psi_3D_output[:,0,2] = beam_parameters[9,:]  + 1j*beam_parameters[15,:] # d (Psi_RZ) / d tau
        # Psi_3D_output[:,1,2] = beam_parameters[10,:] + 1j*beam_parameters[16,:] # d (Psi_zetaZ) / d tau
        # Psi_3D_output[:,1,0] = Psi_3D_output[:,0,1]
        # Psi_3D_output[:,2,0] = Psi_3D_output[:,0,2]
        # Psi_3D_output[:,2,1] = Psi_3D_output[:,1,2]
        
    else:
        """
        - Propagates another ray to find the cut-off location
        - suffix _fine for variables in this section
        - I guess I could've used 'dense_output' for the section above and got the information from there, and one day I'll go and see which method works better/faster
        """
        ray_parameters_2D_initial_fine = ray_parameters_2D[:,index_cutoff_estimate-1]
        tau_start_fine                 = tau_ray[index_cutoff_estimate-1]
        tau_end_fine                   = tau_ray[index_cutoff_estimate+1]
        tau_points_fine                = np.linspace(tau_start_fine,tau_end_fine,1001)
        
        solver_ray_output_fine = integrate.solve_ivp(
                            ray_evolution_2D_fun, [tau_start_fine,tau_end_fine], ray_parameters_2D_initial_fine, 
                            method='RK45', t_eval=tau_points_fine, dense_output=False, 
                            events=event_leave_plasma, vectorized=False, args=solver_arguments
                        )
        ray_parameters_2D_fine = solver_ray_output_fine.y  
        K_magnitude_ray_fine = (  (np.real(ray_parameters_2D_fine[2,:]))**2 
                                + K_zeta_initial**2/(np.real(ray_parameters_2D_fine[0,:]))**2 
                                + (np.real(ray_parameters_2D_fine[3,:]))**2)**(0.5)
        index_cutoff_fine = K_magnitude_ray_fine.argmin()
        tau_ray_fine    = solver_ray_output_fine.t
        tau_cutoff_fine = tau_ray_fine[index_cutoff_fine]        
        
        
        """
        - Propagates the beam
        - But first, tell the solver which values of tau we want for the output
        """        
        tau_points = np.linspace(0,0.995*tau_leave,101).tolist() 
        # I'm using a list because it makes the insertion of tau_cutoff_fine easier
        # factor of 0.995 so that the last point would definitely still be in the plasma
        bisect.insort(tau_points,tau_cutoff_fine) # This function modifies tau_points via a side-effect
        tau_points = np.array(tau_points)

        solver_start_time = time.time()
    
        solver_beam_output = integrate.solve_ivp(
                            beam_evolution_fun, [0,tau_leave], beam_parameters_initial, 
                            method='RK45', t_eval=tau_points, dense_output=False, 
                            events=None, vectorized=False, args=solver_arguments
                        )
        
        solver_end_time = time.time()
        print('Time taken (ray solver)', solver_end_time-solver_start_time,'s')    
        
        beam_parameters = solver_beam_output.y
        tau_array = solver_beam_output.t
        solver_status = solver_beam_output.status
        
        numberOfDataPoints = len(tau_array)
        
        q_R_array = np.real(beam_parameters[0,:])
        q_zeta_array = np.real(beam_parameters[1,:])
        q_Z_array = np.real(beam_parameters[2,:])
        
        K_R_array = np.real(beam_parameters[3,:])
        K_Z_array = np.real(beam_parameters[4,:])
        
        Psi_3D_output = np.zeros([numberOfDataPoints,3,3],dtype='complex128')
        Psi_3D_output[:,0,0] = beam_parameters[5,:]  + 1j*beam_parameters[11,:] # d (Psi_RR) / d tau
        Psi_3D_output[:,1,1] = beam_parameters[6,:]  + 1j*beam_parameters[12,:] # d (Psi_zetazeta) / d tau
        Psi_3D_output[:,2,2] = beam_parameters[7,:]  + 1j*beam_parameters[13,:] # d (Psi_ZZ) / d tau
        Psi_3D_output[:,0,1] = beam_parameters[8,:]  + 1j*beam_parameters[14,:] # d (Psi_Rzeta) / d tau
        Psi_3D_output[:,0,2] = beam_parameters[9,:]  + 1j*beam_parameters[15,:] # d (Psi_RZ) / d tau
        Psi_3D_output[:,1,2] = beam_parameters[10,:] + 1j*beam_parameters[16,:] # d (Psi_zetaZ) / d tau
        Psi_3D_output[:,1,0] = Psi_3D_output[:,0,1]
        Psi_3D_output[:,2,0] = Psi_3D_output[:,0,2]
        Psi_3D_output[:,2,1] = Psi_3D_output[:,1,2]

    print('Main loop complete')
    # -------------------


    ## -------------------
    ## This saves the data generated by the main loop and the input data
    ## -------------------
    print('Saving data')    
    np.savez('data_input' + output_filename_suffix, 
              data_poloidal_flux_grid=data_poloidal_flux_grid,
              data_R_coord=data_R_coord, data_Z_coord=data_Z_coord,
              poloidal_launch_angle_Torbeam=poloidal_launch_angle_Torbeam,
              toroidal_launch_angle_Torbeam=toroidal_launch_angle_Torbeam,
              launch_freq_GHz=launch_freq_GHz,
              mode_flag=mode_flag,
              launch_beam_width=launch_beam_width,
              launch_beam_radius_of_curvature=launch_beam_radius_of_curvature,
              launch_position=launch_position,
              launch_K=launch_K,
              ne_data_density_array=ne_data_density_array,ne_data_radialcoord_array=ne_data_radialcoord_array,
              efit_time=efit_time
             )    
    np.savez('solver_output' + output_filename_suffix, 
             solver_status=solver_status,
             tau_array=tau_array,
             q_R_array=q_R_array,
             q_zeta_array = q_zeta_array,
             q_Z_array = q_Z_array,
             K_R_array = K_R_array,
             K_Z_array = K_Z_array,
             Psi_3D_output = Psi_3D_output
             )   
    
    if solver_status == -1:
        # If the solver doesn't finish, end the function here
        print('Solver did not reach completion')
        return
    
    ## -------------------
    ## Generates additional data along the path of the beam
    ## -------------------

        # Calculate various properties along the ray
    poloidal_flux_output = interp_poloidal_flux(q_R_array,q_Z_array,grid=False)
    electron_density_output = find_density_1D(poloidal_flux_output)

        # Calculates nabla_K H
    dH_dKR_output    = find_dH_dKR(q_R_array, q_Z_array, K_R_array, K_zeta_initial, K_Z_array,
                                   launch_angular_frequency, mode_flag,
                                   delta_K_R,
                                   interp_poloidal_flux, find_density_1D,
                                   find_B_R, find_B_T, find_B_Z
                                   )
    dH_dKzeta_output = find_dH_dKzeta(q_R_array, q_Z_array, K_R_array, K_zeta_initial, K_Z_array,
                                      launch_angular_frequency, mode_flag, 
                                      delta_K_zeta,
                                      interp_poloidal_flux, find_density_1D,
                                      find_B_R, find_B_T, find_B_Z
                                      )
    dH_dKZ_output    = find_dH_dKZ(q_R_array, q_Z_array, K_R_array, K_zeta_initial, K_Z_array,
                                   launch_angular_frequency, mode_flag,
                                   delta_K_Z,
                                   interp_poloidal_flux, find_density_1D,
                                   find_B_R, find_B_T, find_B_Z
                                   )

    
        # Calculates g_hat
    g_hat_output = np.zeros([numberOfDataPoints,3])
    g_magnitude_output = (q_R_array**2 * dH_dKzeta_output**2 + dH_dKR_output**2 + dH_dKZ_output**2)**0.5
    g_hat_output[:,0] = dH_dKR_output                / g_magnitude_output # g_hat_R
    g_hat_output[:,1] = q_R_array * dH_dKzeta_output / g_magnitude_output # g_hat_zeta
    g_hat_output[:,2] = dH_dKZ_output                / g_magnitude_output # g_hat_Z

        # Calculates b_hat and grad_b_hat
    b_hat_output = np.zeros([numberOfDataPoints,3])
    B_R_output = find_B_R(q_R_array,q_Z_array)
    B_T_output = find_B_T(q_R_array,q_Z_array)
    B_Z_output = find_B_Z(q_R_array,q_Z_array)
    B_magnitude = np.sqrt(B_R_output**2 + B_T_output**2 + B_Z_output**2)
    b_hat_output[:,0] = B_R_output / B_magnitude
    b_hat_output[:,1] = B_T_output / B_magnitude
    b_hat_output[:,2] = B_Z_output / B_magnitude

    grad_bhat_output = np.zeros([numberOfDataPoints,3,3])
    dbhat_dR = find_dbhat_dR(q_R_array, q_Z_array, delta_R, find_B_R, find_B_T, find_B_Z)
    dbhat_dZ = find_dbhat_dZ(q_R_array, q_Z_array, delta_Z, find_B_R, find_B_T, find_B_Z)
    grad_bhat_output[:,0,:] = dbhat_dR.T # Transpose dbhat_dR so that it has the right shape
    grad_bhat_output[:,2,:] = dbhat_dZ.T
    grad_bhat_output[:,1,0] = - B_T_output / (B_magnitude * q_R_array)
    grad_bhat_output[:,1,1] = B_R_output / (B_magnitude * q_R_array)


    # x_hat and y_hat
    y_hat_output = np.zeros([numberOfDataPoints,3])
    x_hat_output = np.zeros([numberOfDataPoints,3])
    y_output = np.cross(b_hat_output,g_hat_output) 
    x_output = np.cross(g_hat_output,y_output) 
    y_output_magnitude = np.linalg.norm(y_output,axis=1)
    x_output_magnitude = np.linalg.norm(x_output,axis=1)
    y_hat_output[:,0] = y_output[:,0] / y_output_magnitude
    y_hat_output[:,1] = y_output[:,1] / y_output_magnitude
    y_hat_output[:,2] = y_output[:,2] / y_output_magnitude
    x_hat_output[:,0] = x_output[:,0] / x_output_magnitude
    x_hat_output[:,1] = x_output[:,1] / x_output_magnitude
    x_hat_output[:,2] = x_output[:,2] / x_output_magnitude

    ## -------------------
    ## Not useful for physics or data analysis
    ## But good for checking whether things are working properly
    ## -------------------
        
        # Gradients of poloidal flux along the ray
    # dpolflux_dR_FFD_debugging   =
    # dpolflux_dZ_FFD_debugging   =
    # d2polflux_dR2_FFD_debugging =
    # d2polflux_dZ2_FFD_debugging = 

        # Gradients of the total magnetic field along the ray
    # dB_dR_FFD_debugging     =
    # dB_dZ_FFD_debugging     =
    # d2B_dR2_FFD_debugging   =
    # d2B_dZ2_FFD_debugging   =
    # d2B_dR_dZ_FFD_debugging =


    ## -------------------
    ## This saves the data generated by the main loop and the input data
    ## Input data saved at this point in case something is changed between loading and the end of the main loop, this allows for comparison
    ## The rest of the data is save further down, after the analysis generates them.
    ## Just in case the analysis fails to run, at least one can get the data from the main loop
    ## -------------------
    if vacuumLaunch_flag:
        # np.savez('data_input' + output_filename_suffix, 
        #           # tau_step=tau_step, 
        #           data_poloidal_flux_grid=data_poloidal_flux_grid,
        #           data_R_coord=data_R_coord, data_Z_coord=data_Z_coord,
        #           poloidal_launch_angle_Torbeam=poloidal_launch_angle_Torbeam,
        #           toroidal_launch_angle_Torbeam=toroidal_launch_angle_Torbeam,
        #           launch_freq_GHz=launch_freq_GHz,
        #           mode_flag=mode_flag,
        #           launch_beam_width=launch_beam_width,
        #           launch_beam_curvature=launch_beam_curvature,
        #           launch_position=launch_position,
        #           launch_K=launch_K,
        #           ne_data_density_array=ne_data_density_array,ne_data_radialcoord_array=ne_data_radialcoord_array
        #           )    
        np.savez('data_output' + output_filename_suffix, 
                  tau_array=tau_array, q_R_array=q_R_array, q_zeta_array=q_zeta_array, q_Z_array=q_Z_array,
                  K_R_array=K_R_array, K_zeta_initial=K_zeta_initial, K_Z_array=K_Z_array,
                  Psi_3D_output=Psi_3D_output, Psi_3D_lab_launch=Psi_3D_lab_launch,
                  Psi_3D_lab_entry=Psi_3D_lab_entry,
                  distance_from_launch_to_entry=distance_from_launch_to_entry,
                  g_hat_output=g_hat_output,g_magnitude_output=g_magnitude_output,
                  B_magnitude=B_magnitude,
                  B_R_output = B_R_output,
                  B_T_output = B_T_output,
                  B_Z_output = B_Z_output,
                  x_hat_output=x_hat_output,y_hat_output=y_hat_output,
                  b_hat_output=b_hat_output,
                  grad_bhat_output=grad_bhat_output,
                  dH_dKR_output=dH_dKR_output,dH_dKzeta_output=dH_dKzeta_output,dH_dKZ_output=dH_dKZ_output,
                  # dH_dR_output=dH_dR_output,dH_dZ_output=dH_dZ_output,
                  # grad_grad_H_output=grad_grad_H_output,gradK_grad_H_output=gradK_grad_H_output,gradK_gradK_H_output=gradK_gradK_H_output,
                  # d_poloidal_flux_dR_output=d_poloidal_flux_dR_output,
                  # d_poloidal_flux_dZ_output=d_poloidal_flux_dZ_output,
                   # epsilon_para_output=epsilon_para_output,epsilon_perp_output=epsilon_perp_output,epsilon_g_output=epsilon_g_output,
                  # electron_density_output=electron_density_output,H_output=H_output,
                   poloidal_flux_output=poloidal_flux_output,
                  # dB_dR_FFD_debugging=dB_dR_FFD_debugging,dB_dZ_FFD_debugging=dB_dZ_FFD_debugging,
                  # d2B_dR2_FFD_debugging=d2B_dR2_FFD_debugging,d2B_dZ2_FFD_debugging=d2B_dZ2_FFD_debugging,d2B_dR_dZ_FFD_debugging=d2B_dR_dZ_FFD_debugging,
                  # poloidal_flux_debugging_1R=poloidal_flux_debugging_1R,
                  # poloidal_flux_debugging_2R=poloidal_flux_debugging_2R,
                  # poloidal_flux_debugging_3R=poloidal_flux_debugging_3R,
                  # poloidal_flux_debugging_1Z=poloidal_flux_debugging_1Z,
                  # poloidal_flux_debugging_2Z=poloidal_flux_debugging_2Z,
                  # poloidal_flux_debugging_3Z=poloidal_flux_debugging_3Z,
                  # poloidal_flux_debugging_2R_2Z=poloidal_flux_debugging_2R_2Z,
                  # electron_density_debugging_1R=electron_density_debugging_1R,
                  # electron_density_debugging_2R=electron_density_debugging_2R,
                  # electron_density_debugging_3R=electron_density_debugging_3R,
                  # electron_density_debugging_1Z=electron_density_debugging_1Z,
                  # electron_density_debugging_2Z=electron_density_debugging_2Z,
                  # electron_density_debugging_3Z=electron_density_debugging_3Z,
                  # electron_density_debugging_2R_2Z=electron_density_debugging_2R_2Z,
                  # dpolflux_dR_FFD_debugging=dpolflux_dR_FFD_debugging,    
                  # dpolflux_dZ_FFD_debugging=dpolflux_dZ_FFD_debugging,
                  # d2polflux_dR2_FFD_debugging=d2polflux_dR2_FFD_debugging,
                  # d2polflux_dZ2_FFD_debugging=d2polflux_dZ2_FFD_debugging, 
                  )
    else:
          np.savez('data_input' + output_filename_suffix, 
                  # tau_step=tau_step, 
                  data_poloidal_flux_grid=data_poloidal_flux_grid,
                  data_R_coord=data_R_coord, data_Z_coord=data_Z_coord,
                  launch_freq_GHz=launch_freq_GHz,
                  mode_flag=mode_flag,
                  launch_position=launch_position,
                  plasmaLaunch_K=plasmaLaunch_K,
                  plasmaLaunch_Psi_3D_lab_Cartesian=plasmaLaunch_Psi_3D_lab_Cartesian,
                  ne_data_density_array=ne_data_density_array,ne_data_radialcoord_array=ne_data_radialcoord_array
                  )  
          np.savez('data_output' + output_filename_suffix, 
                  tau_array=tau_array, q_R_array=q_R_array, q_zeta_array=q_zeta_array, q_Z_array=q_Z_array,
                  K_R_array=K_R_array, K_zeta_initial=K_zeta_initial, K_Z_array=K_Z_array,
                  Psi_3D_output=Psi_3D_output, Psi_3D_lab_launch=Psi_3D_lab_launch,
                  g_hat_output=g_hat_output,g_magnitude_output=g_magnitude_output,
                  B_magnitude=B_magnitude,
                  B_R_output = B_R_output,
                  B_T_output = B_T_output,
                  B_Z_output = B_Z_output,
                  x_hat_output=x_hat_output,y_hat_output=y_hat_output,
                  b_hat_output=b_hat_output,
                  grad_bhat_output=grad_bhat_output,
                  dH_dKR_output=dH_dKR_output,dH_dKzeta_output=dH_dKzeta_output,dH_dKZ_output=dH_dKZ_output,
                  # dH_dR_output=dH_dR_output,dH_dZ_output=dH_dZ_output,
                  # grad_grad_H_output=grad_grad_H_output,gradK_grad_H_output=gradK_grad_H_output,gradK_gradK_H_output=gradK_gradK_H_output,
                  # d_poloidal_flux_dR_output=d_poloidal_flux_dR_output,
                  # d_poloidal_flux_dZ_output=d_poloidal_flux_dZ_output,
                  # epsilon_para_output=epsilon_para_output,epsilon_perp_output=epsilon_perp_output,epsilon_g_output=epsilon_g_output,
                  # electron_density_output=electron_density_output,H_output=H_output
                  )     
        
    print('Data saved')
    # -------------------



    ## -------------------
    ## Process the data from the main loop to give a bunch of useful stuff
    ## -------------------
    print('Analysing data')

        # Calculates various useful stuff
    [q_X_array,q_Y_array,_] = find_q_lab_Cartesian([q_R_array,q_zeta_array,q_Z_array])
    point_spacing =  ( (np.diff(q_X_array))**2 + (np.diff(q_Y_array))**2 + (np.diff(q_Z_array))**2  )**0.5
    distance_along_line =  np.cumsum(point_spacing)
    distance_along_line = np.append(0,distance_along_line)
    RZ_point_spacing = ( (np.diff(q_Z_array))**2 + (np.diff(q_R_array))**2  )**0.5
    RZ_distance_along_line =  np.cumsum(RZ_point_spacing)
    RZ_distance_along_line = np.append(0,RZ_distance_along_line)
    
        # Calculates the index of the minimum magnitude of K
        # That is, finds when the beam hits the cut-off        
    K_magnitude_array = (K_R_array**2 + K_zeta_initial**2/q_R_array**2 + K_Z_array**2)**(0.5)        
        
    cutoff_index = find_nearest(abs(K_magnitude_array),  0) # Index of the cutoff, at the minimum value of K, use this with other arrays

    cyclotron_freq_output = launch_angular_frequency * find_normalised_gyro_freq(B_magnitude, launch_angular_frequency)

    
#         # Calculates when the beam 'enters' and 'leaves' the plasma
#         # Here, entry and exit refer to crossing poloidal_flux_enter, if specified
#         # Otherwise, entry and exit refer to crossing the LCFS, poloidal_flux = 1.0
#     if poloidal_flux_enter is None:
#         in_out_poloidal_flux = 1.0
#     else:
#         in_out_poloidal_flux = poloidal_flux_enter
#     poloidal_flux_a = poloidal_flux_output[0:cutoff_index]
#     poloidal_flux_b = poloidal_flux_output[cutoff_index::]
#     in_index = find_nearest(poloidal_flux_a,in_out_poloidal_flux)
#     out_index = cutoff_index + find_nearest(poloidal_flux_b,in_out_poloidal_flux)


    
        # Calcuating the corrections to make M from Psi
        # Corrections that are small in mismatch are ignored in the calculation of M
        # However, these small terms are still calculated in this section, so that their actual size can be checked, if need be
    ray_curvature_kappa_output = np.zeros([numberOfDataPoints,3])

    k_perp_1_backscattered = -2*K_magnitude_array

    y_hat_Cartesian = np.zeros([numberOfDataPoints,3])
    x_hat_Cartesian = np.zeros([numberOfDataPoints,3])
    y_hat_Cartesian[:,0] = y_hat_output[:,0]*np.cos(q_zeta_array ) - y_hat_output[:,1]*np.sin(q_zeta_array )
    y_hat_Cartesian[:,1] = y_hat_output[:,0]*np.sin(q_zeta_array ) + y_hat_output[:,1]*np.cos(q_zeta_array )
    y_hat_Cartesian[:,2] = y_hat_output[:,2]
    x_hat_Cartesian[:,0] = x_hat_output[:,0]*np.cos(q_zeta_array ) - x_hat_output[:,1]*np.sin(q_zeta_array )
    x_hat_Cartesian[:,1] = x_hat_output[:,0]*np.sin(q_zeta_array ) + x_hat_output[:,1]*np.cos(q_zeta_array )
    x_hat_Cartesian[:,2] = x_hat_output[:,2]
    
    Psi_3D_Cartesian = find_Psi_3D_lab_Cartesian(Psi_3D_output, q_R_array, q_zeta_array, K_R_array, K_zeta_initial)
    Psi_xx_output = contract_special(x_hat_Cartesian,contract_special(Psi_3D_Cartesian,x_hat_Cartesian))
    Psi_xy_output = contract_special(x_hat_Cartesian,contract_special(Psi_3D_Cartesian,y_hat_Cartesian))
    Psi_yy_output = contract_special(y_hat_Cartesian,contract_special(Psi_3D_Cartesian,y_hat_Cartesian))
    
    Psi_xx_entry = np.dot(x_hat_Cartesian[0,:],np.dot(Psi_3D_lab_entry_cartersian,x_hat_Cartesian[0,:]))
    Psi_xy_entry = np.dot(x_hat_Cartesian[0,:],np.dot(Psi_3D_lab_entry_cartersian,y_hat_Cartesian[0,:]))
    Psi_yy_entry = np.dot(y_hat_Cartesian[0,:],np.dot(Psi_3D_lab_entry_cartersian,y_hat_Cartesian[0,:]))

    xhat_dot_grad_bhat = contract_special(x_hat_output,grad_bhat_output)
    yhat_dot_grad_bhat = contract_special(y_hat_output,grad_bhat_output)
    ray_curvature_kappa_output[:,0] = (1/g_magnitude_output) * np.gradient(g_hat_output[:,0],tau_array)
    ray_curvature_kappa_output[:,1] = (1/g_magnitude_output) * np.gradient(g_hat_output[:,1],tau_array)
    ray_curvature_kappa_output[:,2] = (1/g_magnitude_output) * np.gradient(g_hat_output[:,2],tau_array)
        
    xhat_dot_grad_bhat_dot_xhat_output = contract_special(xhat_dot_grad_bhat,x_hat_output)
    xhat_dot_grad_bhat_dot_yhat_output = contract_special(xhat_dot_grad_bhat,y_hat_output)
    xhat_dot_grad_bhat_dot_ghat_output = contract_special(xhat_dot_grad_bhat,g_hat_output)
    yhat_dot_grad_bhat_dot_xhat_output = contract_special(yhat_dot_grad_bhat,x_hat_output)
    yhat_dot_grad_bhat_dot_yhat_output = contract_special(yhat_dot_grad_bhat,y_hat_output)
    yhat_dot_grad_bhat_dot_ghat_output = contract_special(yhat_dot_grad_bhat,g_hat_output)
    kappa_dot_xhat_output = contract_special(ray_curvature_kappa_output,xhat_dot_grad_bhat)
    kappa_dot_yhat_output = contract_special(ray_curvature_kappa_output,yhat_dot_grad_bhat)
    # TODO: Calculate the other small corrections
    
    M_xx_output = Psi_xx_output + (k_perp_1_backscattered/2) * xhat_dot_grad_bhat_dot_ghat_output
    M_xy_output = Psi_xy_output + (k_perp_1_backscattered/2) * yhat_dot_grad_bhat_dot_ghat_output
    M_yy_output = Psi_yy_output
    
        # Calculates the localisation, wavenumber resolution, and mismatch attenuation pieces
    det_M_w_analysis = M_xx_output*M_yy_output - M_xy_output**2
    M_w_inv_xx_output =   M_yy_output / det_M_w_analysis
    M_w_inv_xy_output = - M_xy_output / det_M_w_analysis
    M_w_inv_yy_output =   M_xx_output / det_M_w_analysis
    
    delta_k_perp_2 = np.sqrt( - 1 / np.imag(M_w_inv_yy_output) )
    delta_theta_m  = np.sqrt( 
                              np.imag(M_w_inv_yy_output) / ( (np.imag(M_w_inv_xy_output))**2 - np.imag(M_w_inv_xx_output)*np.imag(M_w_inv_yy_output) ) 
                            ) / (2 * K_magnitude_array)
    
    sin_theta_m_analysis = np.zeros(numberOfDataPoints)
    sin_theta_m_analysis[:] = (b_hat_output[:,0]*K_R_array[:] + b_hat_output[:,1]*K_zeta_initial/q_R_array[:] + b_hat_output[:,2]*K_Z_array[:]) / (K_magnitude_array[:]) # B \cdot K / (abs (B) abs(K))

    theta_m_output = np.sign(sin_theta_m_analysis)*np.arcsin(abs(sin_theta_m_analysis)) # Assumes the mismatch angle is never smaller than -90deg or bigger than 90deg
    
    print('theta_m', theta_m_output[cutoff_index])
    print('mismatch attenuation', np.exp(-theta_m_output[cutoff_index]**2/ (2*delta_theta_m[cutoff_index]**2)))
#    print(cutoff_index)


        # This part is used to make some nice plots when post-processing
    R_midplane_points = np.linspace(data_R_coord[0],data_R_coord[-1],1000)
    poloidal_flux_on_midplane = interp_poloidal_flux(R_midplane_points,0,grid=False) # poloidal flux at R and z=0
        
    
        # Calculates localisation (start)
    # Ray piece of localisation as a function of distance along ray
    K_magnitude_array_plus_KR     = np.sqrt((K_R_array+delta_K_R)**2 + K_Z_array**2 + K_zeta_initial**2/q_R_array**2)
    K_magnitude_array_minus_KR    = np.sqrt((K_R_array-delta_K_R)**2 + K_Z_array**2 + K_zeta_initial**2/q_R_array**2)
    K_magnitude_array_plus_Kzeta  = np.sqrt(K_R_array**2 + K_Z_array**2 + (K_zeta_initial+delta_K_zeta)**2/q_R_array**2)
    K_magnitude_array_minus_Kzeta = np.sqrt(K_R_array**2 + K_Z_array**2 + (K_zeta_initial-delta_K_zeta)**2/q_R_array**2)
    K_magnitude_array_plus_KZ     = np.sqrt(K_R_array**2 + (K_Z_array+delta_K_Z)**2 + K_zeta_initial**2/q_R_array**2)
    K_magnitude_array_minus_KZ    = np.sqrt(K_R_array**2 + (K_Z_array-delta_K_Z)**2 + K_zeta_initial**2/q_R_array**2)
    
    epsilon_para_output = find_epsilon_para(electron_density_output, launch_angular_frequency)
    epsilon_perp_output = find_epsilon_perp(electron_density_output, B_magnitude, launch_angular_frequency)
    epsilon_g_output = find_epsilon_g(electron_density_output, B_magnitude, launch_angular_frequency) 
    
    H_1_Cardano_array,H_2_Cardano_array,H_3_Cardano_array = find_H_Cardano(K_magnitude_array,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
    
    # In my experience, the H_3_Cardano expression corresponds to the O mode, and the H_2_Cardano expression corresponds to the X-mode
    # ALERT: This may not always be the case! Check the output figure to make sure that the appropriate solution is indeed 0 along the ray
    if mode_flag == 1:
        _,_,H_Cardano_plus_KR_array     = find_H_Cardano(K_magnitude_array_plus_KR    ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,_,H_Cardano_minus_KR_array    = find_H_Cardano(K_magnitude_array_minus_KR   ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,_,H_Cardano_plus_Kzeta_array  = find_H_Cardano(K_magnitude_array_plus_Kzeta ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,_,H_Cardano_minus_Kzeta_array = find_H_Cardano(K_magnitude_array_minus_Kzeta,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,_,H_Cardano_plus_KZ_array     = find_H_Cardano(K_magnitude_array_plus_KZ    ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,_,H_Cardano_minus_KZ_array    = find_H_Cardano(K_magnitude_array_minus_KZ   ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
    elif mode_flag == -1:
        _,H_Cardano_plus_KR_array,_     = find_H_Cardano(K_magnitude_array_plus_KR    ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,H_Cardano_minus_KR_array,_    = find_H_Cardano(K_magnitude_array_minus_KR   ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,H_Cardano_plus_Kzeta_array,_  = find_H_Cardano(K_magnitude_array_plus_Kzeta ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,H_Cardano_minus_Kzeta_array,_ = find_H_Cardano(K_magnitude_array_minus_Kzeta,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,H_Cardano_plus_KZ_array,_     = find_H_Cardano(K_magnitude_array_plus_KZ    ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)
        _,H_Cardano_minus_KZ_array,_    = find_H_Cardano(K_magnitude_array_minus_KZ   ,launch_angular_frequency,epsilon_para_output,epsilon_perp_output,epsilon_g_output,theta_m_output)

    g_R_Cardano    = np.real(H_Cardano_plus_KR_array - H_Cardano_minus_KR_array) / (2 * delta_K_R)
    g_zeta_Cardano = np.real(H_Cardano_plus_Kzeta_array - H_Cardano_minus_Kzeta_array) / (2 * delta_K_zeta)
    g_Z_Cardano    = np.real(H_Cardano_plus_KZ_array - H_Cardano_minus_KZ_array) / (2 * delta_K_Z)
    
    g_magnitude_Cardano = np.sqrt(g_R_Cardano**2 + g_zeta_Cardano**2 + g_Z_Cardano**2)  
    
    ##
    # From here on, we use the shorthand
        # loc: localisation
        # l_lc: distance from cutoff (l - l_c). Distance along the ray
        # cum: cumulative. As such, cum_loc is the cumulative integral of the localisation
        # p: polarisation    
        # r: ray
        # b: beam
        # s: spectrum
    # Otherwise, variable names get really unwieldly
    ##
    
    # localisation_ray = g_magnitude_Cardano[0]**2/g_magnitude_Cardano**2
        # The first point of the beam may be very slightly in the plasma, so I have used the vacuum expression for the group velocity instead
    loc_r = (2 * constants.c / launch_angular_frequency)**2 / g_magnitude_Cardano**2
    
    # Spectrum piece of localisation as a function of distance along ray      
    spectrum_power_law_coefficient = 13/3 # Turbulence cascade
    loc_s = ( k_perp_1_backscattered / (-2*wavenumber_K0) )**(-spectrum_power_law_coefficient)
    
    # Beam piece of localisation as a function of distance along ray   
    det_imag_Psi_w_analysis = np.imag(Psi_xx_output)*np.imag(Psi_yy_output) - np.imag(Psi_xy_output)**2 # Determinant of the imaginary part of Psi_w
    det_real_Psi_w_analysis = np.real(Psi_xx_output)*np.real(Psi_yy_output) - np.real(Psi_xy_output)**2 # Determinant of the real part of Psi_w. Not needed for the calculation, but gives useful insight

    beam_waist_y = find_waist(launch_beam_width, wavenumber_K0, 1/launch_beam_radius_of_curvature) # Assumes circular beam at launch
    
    loc_b = (beam_waist_y/np.sqrt(2))* det_imag_Psi_w_analysis / ( abs(det_M_w_analysis) * np.sqrt(-np.imag(M_w_inv_yy_output)) )
    # --
    
    # Polarisation piece of localisation as a function of distance along ray   
    # Polarisation e
        # eigenvector corresponding to eigenvalue = 0 (H=0)
        # First, find the components of the tensor D
        # Refer to 21st Dec 2020 notes for more
        # Note that e \cdot e* = 1
    [D_11_component, D_22_component, 
     D_bb_component, D_12_component, 
     D_1b_component] = find_D(K_magnitude_array,launch_angular_frequency,
                              epsilon_para_output,epsilon_perp_output,epsilon_g_output,
                              theta_m_output)
    
        # Dispersion tensor
    D_tensor = np.zeros([numberOfDataPoints,3,3],dtype='complex128')
    D_tensor[:,0,0] = D_11_component
    D_tensor[:,1,1] = D_22_component
    D_tensor[:,2,2] = D_bb_component
    D_tensor[:,0,1] = -1j * D_12_component
    D_tensor[:,1,0] =  1j * D_12_component                              
    D_tensor[:,0,2] = D_1b_component
    D_tensor[:,2,0] = D_1b_component


    H_eigvals, e_eigvecs = np.linalg.eigh(D_tensor)

    # In my experience, H_eigvals[:,2] corresponds to the O mode, and H_eigvals[:,1] corresponds to the X-mode
    # ALERT: This may not always be the case! Check the output figure to make sure that the appropriate solution is indeed 0 along the ray
    if mode_flag == 1:
        H_solver = H_eigvals[:,1]
        e_hat_output = e_eigvecs[:,:,1]
    elif mode_flag == -1:
        H_solver = H_eigvals[:,0]
        e_hat_output = e_eigvecs[:,:,0]

      # equilibrium dielectric tensor - identity matrix. \bm{\epsilon}_{eq} - \bm{1}
    epsilon_minus_identity = np.zeros([numberOfDataPoints,3,3],dtype='complex128')
    epsilon_minus_identity[:,0,0] =    epsilon_perp_output - np.ones(numberOfDataPoints)
    epsilon_minus_identity[:,1,1] =    epsilon_perp_output - np.ones(numberOfDataPoints)
    epsilon_minus_identity[:,2,2] =    epsilon_para_output - np.ones(numberOfDataPoints)
    epsilon_minus_identity[:,0,1] = -1j * epsilon_g_output
    epsilon_minus_identity[:,1,0] =  1j * epsilon_g_output

    loc_p = abs(contract_special(np.conjugate(e_hat_output), contract_special(epsilon_minus_identity,e_hat_output)))**2 / (electron_density_output*10**19)**2
    
    # Note that K_1 = K cos theta_m, K_2 = 0, K_b = K sin theta_m, as a result of cold plasma dispersion
    K_hat_dot_e_hat = (
          e_hat_output[:,0] * np.cos(theta_m_output)
        + e_hat_output[:,2] * np.sin(theta_m_output)
        )

    K_hat_dot_e_hat_sq = np.conjugate(K_hat_dot_e_hat) * K_hat_dot_e_hat    
    # --
    
    ## TODO: Come back and see if the naming of variables makes sense and is consistent
    
    l_lc = distance_along_line-distance_along_line[cutoff_index] # Distance from cutoff
    
    
    # plt.figure()
    # plt.plot(l_lc,theta_m_output)
    # plt.axhline(constants.e**4 / (launch_angular_frequency**2 *constants.epsilon_0*constants.m_e)**2,c='k')
    
        # Combining the various localisation pieces to get some overall localisation
    # loc_p_r_s   =                  loc_p * loc_r * loc_s
    loc_b_r_s = loc_b * loc_r * loc_s
    loc_b_r   = loc_b * loc_r
    
        # Finds the 1/e2 values (localisation)
    loc_b_r_s_max_over_e2 = loc_b_r_s.max() / (np.e)**2 # loc_b_r_s.max() / 2.71**2
    loc_b_r_max_over_e2   = loc_b_r.max() / (np.e)**2 # loc_b_r.max() / 2.71**2
    
        # Gives the inter-e2 range (analogous to interquartile range) in l-lc
    loc_b_r_s_delta_l_1 = find_x0(l_lc[0:cutoff_index], loc_b_r_s[0:cutoff_index], loc_b_r_s_max_over_e2)
    loc_b_r_s_delta_l_2 = find_x0(l_lc[cutoff_index::], loc_b_r_s[cutoff_index::], loc_b_r_s_max_over_e2)
    loc_b_r_s_delta_l = np.array([loc_b_r_s_delta_l_1, loc_b_r_s_delta_l_2])# The 1/e2 distances,  (l - l_c)
    loc_b_r_s_half_width_l = (loc_b_r_s_delta_l_2 - loc_b_r_s_delta_l_1)/2
    loc_b_r_delta_l_1 = find_x0(l_lc[0:cutoff_index], loc_b_r[0:cutoff_index], loc_b_r_max_over_e2)
    loc_b_r_delta_l_2 = find_x0(l_lc[cutoff_index::], loc_b_r[cutoff_index::], loc_b_r_max_over_e2)
    loc_b_r_delta_l = np.array([loc_b_r_delta_l_1, loc_b_r_delta_l_2])# The 1/e2 distances,  (l - l_c)   
    loc_b_r_half_width_l = (loc_b_r_delta_l_1 - loc_b_r_delta_l_2)/2

        # Estimates the inter-e2 range (analogous to interquartile range) in kperp1, from l-lc
        # Bear in mind that since abs(kperp1) is minimised at cutoff, one really has to use that in addition to these.
    loc_b_r_s_delta_kperp1_1 = find_x0(k_perp_1_backscattered[0:cutoff_index],l_lc[0:cutoff_index],loc_b_r_s_delta_l_1)
    loc_b_r_s_delta_kperp1_2 = find_x0(k_perp_1_backscattered[cutoff_index::], l_lc[cutoff_index::], loc_b_r_s_delta_l_2)
    loc_b_r_s_delta_kperp1 = np.array([loc_b_r_s_delta_kperp1_1,loc_b_r_s_delta_kperp1_2])
    loc_b_r_delta_kperp1_1 = find_x0(k_perp_1_backscattered[0:cutoff_index], l_lc[0:cutoff_index], loc_b_r_delta_l_1)
    loc_b_r_delta_kperp1_2 = find_x0(k_perp_1_backscattered[cutoff_index::], l_lc[cutoff_index::], loc_b_r_delta_l_2)
    loc_b_r_delta_kperp1 = np.array([loc_b_r_delta_kperp1_1,loc_b_r_delta_kperp1_2])


        # Calculate the cumulative integral of the localisation pieces
    cum_loc_b_r_s = integrate.cumtrapz(loc_b_r_s, distance_along_line, initial=0)    
    cum_loc_b_r_s = (cum_loc_b_r_s - max(cum_loc_b_r_s)/2)
    cum_loc_b_r = integrate.cumtrapz(loc_b_r, distance_along_line, initial=0)
    cum_loc_b_r = (cum_loc_b_r - max(cum_loc_b_r)/2)
    
        # Finds the 1/e2 values (cumulative integral of localisation)
    # cum_loc_b_r_s_max_over_e2_1 = cum_loc_b_r_s.min() * (1 - 1 / (np.e)**2)
    # cum_loc_b_r_s_max_over_e2_2 = cum_loc_b_r_s.max() * (1 - 1 / (np.e)**2)
    # cum_loc_b_r_max_over_e2_1 = cum_loc_b_r.min() * (1 - 1 / (np.e)**2)
    # cum_loc_b_r_max_over_e2_2 = cum_loc_b_r.max() * (1 - 1 / (np.e)**2)
    cum_loc_b_r_s_max_over_e2 = cum_loc_b_r_s.max() * (1 - 1 / (np.e)**2)
    cum_loc_b_r_max_over_e2 = cum_loc_b_r.max() * (1 - 1 / (np.e)**2)
    
        # Gives the inter-e range (analogous to interquartile range) in l-lc
    cum_loc_b_r_s_delta_l_1 = find_x0(l_lc, cum_loc_b_r_s, -cum_loc_b_r_s_max_over_e2)
    cum_loc_b_r_s_delta_l_2 = find_x0(l_lc, cum_loc_b_r_s,  cum_loc_b_r_s_max_over_e2)
    cum_loc_b_r_s_delta_l = np.array([cum_loc_b_r_s_delta_l_1,cum_loc_b_r_s_delta_l_2])
    cum_loc_b_r_s_half_width = (cum_loc_b_r_s_delta_l_2 - cum_loc_b_r_s_delta_l_1)/2
    cum_loc_b_r_delta_l_1 = find_x0(l_lc, cum_loc_b_r, -cum_loc_b_r_max_over_e2)
    cum_loc_b_r_delta_l_2 = find_x0(l_lc, cum_loc_b_r,  cum_loc_b_r_max_over_e2)
    cum_loc_b_r_delta_l = np.array([cum_loc_b_r_delta_l_1,cum_loc_b_r_delta_l_2])
    cum_loc_b_r_half_width = (cum_loc_b_r_delta_l_2 - cum_loc_b_r_delta_l_1)/2

        # Gives the inter-e2 range (analogous to interquartile range) in kperp1. 
        # Bear in mind that since abs(kperp1) is minimised at cutoff, one really has to use that in addition to these.
    cum_loc_b_r_s_delta_kperp1_1 = find_x0(k_perp_1_backscattered[0:cutoff_index], cum_loc_b_r_s[0:cutoff_index], -cum_loc_b_r_s_max_over_e2)
    cum_loc_b_r_s_delta_kperp1_2 = find_x0(k_perp_1_backscattered[cutoff_index::], cum_loc_b_r_s[cutoff_index::],  cum_loc_b_r_s_max_over_e2)
    cum_loc_b_r_s_delta_kperp1   = np.array([cum_loc_b_r_s_delta_kperp1_1,cum_loc_b_r_s_delta_kperp1_2])
    cum_loc_b_r_delta_kperp1_1 = find_x0(k_perp_1_backscattered[0:cutoff_index], cum_loc_b_r[0:cutoff_index], -cum_loc_b_r_max_over_e2)
    cum_loc_b_r_delta_kperp1_2 = find_x0(k_perp_1_backscattered[cutoff_index::], cum_loc_b_r[cutoff_index::],  cum_loc_b_r_max_over_e2)
    cum_loc_b_r_delta_kperp1   = np.array([cum_loc_b_r_delta_kperp1_1,cum_loc_b_r_delta_kperp1_2])
    
        # Gives the mode l-lc for backscattering        
    loc_b_r_s_max_index = find_nearest(loc_b_r_s, loc_b_r_s.max())    
    loc_b_r_s_max_l_lc = distance_along_line[loc_b_r_s_max_index] - distance_along_line[cutoff_index]
    loc_b_r_max_index = find_nearest(loc_b_r, loc_b_r.max())    
    loc_b_r_max_l_lc = distance_along_line[loc_b_r_max_index] - distance_along_line[cutoff_index]
    
        # Gives the mean l-lc for backscattering    
    cum_loc_b_r_s_mean_l_lc = np.trapz(loc_b_r_s*distance_along_line, distance_along_line) / np.trapz(loc_b_r_s, distance_along_line) - distance_along_line[cutoff_index]
    cum_loc_b_r_mean_l_lc = np.trapz(loc_b_r*distance_along_line, distance_along_line) / np.trapz(loc_b_r, distance_along_line) - distance_along_line[cutoff_index]

        # Gives the median l-lc for backscattering
    cum_loc_b_r_s_delta_l_0 = find_x0(l_lc, cum_loc_b_r_s, 0) 
    cum_loc_b_r_delta_l_0 = find_x0(l_lc, cum_loc_b_r, 0)

        # Due to the divergency of the ray piece, the mode kperp1 for backscattering is exactly that at the cut-off

        # Gives the mean kperp1 for backscattering    
    cum_loc_b_r_s_mean_kperp1 = np.trapz(loc_b_r_s*k_perp_1_backscattered, k_perp_1_backscattered) / np.trapz(loc_b_r_s, k_perp_1_backscattered)
    cum_loc_b_r_mean_kperp1   = np.trapz(loc_b_r  *k_perp_1_backscattered, k_perp_1_backscattered) / np.trapz(loc_b_r  , k_perp_1_backscattered)

        # Gives the median kperp1 for backscattering
    cum_loc_b_r_s_delta_kperp1_0 = find_x0(k_perp_1_backscattered, cum_loc_b_r_s, 0)    
    cum_loc_b_r_delta_kperp1_0 = find_x0(k_perp_1_backscattered[0:cutoff_index], cum_loc_b_r[0:cutoff_index], 0) # Only works if point is before cutoff. To fix.

        # To make the plots look nice
    k_perp_1_backscattered_plot = np.append(-2*wavenumber_K0, k_perp_1_backscattered)    
    k_perp_1_backscattered_plot = np.append(k_perp_1_backscattered_plot, -2*wavenumber_K0)
    cum_loc_b_r_s_plot = np.append(cum_loc_b_r_s[0], cum_loc_b_r_s)
    cum_loc_b_r_s_plot = np.append(cum_loc_b_r_s_plot, cum_loc_b_r_s[-1])
    cum_loc_b_r_plot = np.append(cum_loc_b_r[0], cum_loc_b_r)
    cum_loc_b_r_plot = np.append(cum_loc_b_r_plot, cum_loc_b_r[-1])
    
    


    


    # integrated_localisation_b_p_r_delta_kperp1_0 = find_x0(k_perp_1_backscattered[0:cutoff_index],integrated_localisation_b_p_r[0:cutoff_index],0)



#     # -------------------




    ## -------------------
    ## This saves the data generated by the analysis after the main loop
    ## -------------------
    print('Saving analysis data')
    np.savez('analysis_output' + output_filename_suffix, 
              Psi_xx_output = Psi_xx_output, Psi_xy_output = Psi_xy_output, Psi_yy_output = Psi_yy_output,
              Psi_xx_entry=Psi_xx_entry, Psi_xy_entry=Psi_xy_entry, Psi_yy_entry=Psi_yy_entry,
              Psi_3D_Cartesian=Psi_3D_Cartesian,x_hat_Cartesian=x_hat_Cartesian,y_hat_Cartesian=y_hat_Cartesian,
              M_xx_output = M_xx_output, M_xy_output = M_xy_output, M_yy_output = M_yy_output,
              M_w_inv_xx_output=M_w_inv_xx_output, M_w_inv_xy_output=M_w_inv_xy_output,M_w_inv_yy_output=M_w_inv_yy_output,
              xhat_dot_grad_bhat_dot_xhat_output=xhat_dot_grad_bhat_dot_xhat_output,
              xhat_dot_grad_bhat_dot_yhat_output=xhat_dot_grad_bhat_dot_yhat_output,
              xhat_dot_grad_bhat_dot_ghat_output=xhat_dot_grad_bhat_dot_ghat_output,
              yhat_dot_grad_bhat_dot_xhat_output=yhat_dot_grad_bhat_dot_xhat_output,
              yhat_dot_grad_bhat_dot_yhat_output=yhat_dot_grad_bhat_dot_yhat_output,
              yhat_dot_grad_bhat_dot_ghat_output=yhat_dot_grad_bhat_dot_ghat_output,
              kappa_dot_xhat_output=kappa_dot_xhat_output,
              kappa_dot_yhat_output=kappa_dot_yhat_output,
              delta_k_perp_2=delta_k_perp_2, delta_theta_m=delta_theta_m,
              theta_m_output=theta_m_output,
              RZ_distance_along_line=RZ_distance_along_line,
              distance_along_line=distance_along_line,
              k_perp_1_backscattered = k_perp_1_backscattered,K_magnitude_array=K_magnitude_array,
              cutoff_index=cutoff_index,
              x_hat_output=x_hat_output,
              y_hat_output=y_hat_output,
              b_hat_output=b_hat_output,
              g_hat_output=g_hat_output,
              # in_index=in_index,out_index=out_index,
              poloidal_flux_on_midplane=poloidal_flux_on_midplane,R_midplane_points=R_midplane_points,
              loc_b=loc_b,loc_p=loc_p,
              loc_r=loc_r,loc_s=loc_s,
              loc_b_r_s=loc_b_r_s, loc_b_r=loc_b_r,
              loc_b_r_s_max_over_e2=loc_b_r_s_max_over_e2,loc_b_r_max_over_e2=loc_b_r_max_over_e2,
              loc_b_r_s_delta_l=loc_b_r_s_delta_l,loc_b_r_delta_l=loc_b_r_delta_l, # The 1/e2 distances,  (l - l_c)
              loc_b_r_s_delta_kperp1=loc_b_r_s_delta_kperp1,loc_b_r_delta_kperp1=loc_b_r_delta_kperp1,  # The 1/e2 distances, kperp1, estimated from (l - l_c)
              cum_loc_b_r_s=cum_loc_b_r_s, cum_loc_b_r=cum_loc_b_r,
              k_perp_1_backscattered_plot=k_perp_1_backscattered_plot,
              cum_loc_b_r_s_plot=cum_loc_b_r_s_plot, cum_loc_b_r_plot=cum_loc_b_r_plot,
              cum_loc_b_r_s_max_over_e2=cum_loc_b_r_s_max_over_e2,cum_loc_b_r_max_over_e2=cum_loc_b_r_max_over_e2,
              cum_loc_b_r_s_delta_l=cum_loc_b_r_s_delta_l,cum_loc_b_r_delta_l=cum_loc_b_r_delta_l, # The cumloc 1/e2 distances, (l - l_c)
              cum_loc_b_r_s_delta_kperp1=cum_loc_b_r_s_delta_kperp1,cum_loc_b_r_delta_kperp1=cum_loc_b_r_delta_kperp1, # The cumloc 1/e2 distances, kperp1
              loc_b_r_s_max_l_lc=loc_b_r_s_max_l_lc, loc_b_r_max_l_lc=loc_b_r_max_l_lc, # mode l-lc
              cum_loc_b_r_s_mean_l_lc=cum_loc_b_r_s_mean_l_lc, cum_loc_b_r_mean_l_lc=cum_loc_b_r_mean_l_lc, # mean l-lc
              cum_loc_b_r_s_delta_l_0=cum_loc_b_r_s_delta_l_0, cum_loc_b_r_delta_l_0=cum_loc_b_r_delta_l_0, # median l-lc
              cum_loc_b_r_s_mean_kperp1=cum_loc_b_r_s_mean_kperp1, cum_loc_b_r_mean_kperp1=cum_loc_b_r_mean_kperp1, # mean kperp1
              cum_loc_b_r_s_delta_kperp1_0=cum_loc_b_r_s_delta_kperp1_0, cum_loc_b_r_delta_kperp1_0=cum_loc_b_r_delta_kperp1_0, # median kperp1
              # det_imag_Psi_w_analysis=det_imag_Psi_w_analysis,det_real_Psi_w_analysis=det_real_Psi_w_analysis,det_M_w_analysis=det_M_w_analysis
              )
    print('Analysis data saved')
    # -------------------

    
    ## -------------------
    ## This saves some simple figures
    ## Allows one to quickly gain an insight into what transpired in the simulation
    ## -------------------
    if figure_flag:
        print('Making figures')
        output_figurename_suffix = output_filename_suffix + '.png'
        
        """
        Plots the beam path on the R Z plane
        """
        plt.figure()
        plt.title('Rz')
        plt.xlabel('R / m') # x-direction
        plt.ylabel('z / m')
    
        contour_levels = np.linspace(0,1.0,11)
        CS = plt.contour(data_R_coord, data_Z_coord, np.transpose(data_poloidal_flux_grid), contour_levels,vmin=0,vmax=1,cmap='plasma_r')
        plt.clabel(CS, inline=1, fontsize=10) # Labels the flux surfaces
        plt.plot(
                np.concatenate([[launch_position[0],initial_position[0]],q_R_array ]),
                np.concatenate([[launch_position[2],initial_position[2]],q_Z_array ]),
                '--.k') # Central (reference) ray
        #cutoff_contour = plt.contour(x_grid, z_grid, normalised_plasma_freq_grid,
        #                             levels=1,vmin=1,vmax=1,linewidths=5,colors='grey')
        plt.xlim(data_R_coord[0],data_R_coord[-1])
        plt.ylim(data_Z_coord[0],data_Z_coord[-1])

        plt.savefig('Ray1_' + output_figurename_suffix)
        plt.close()
        
        """
        Plots Cardano's and np.linalg's solutions to the actual dispersion relation
        Useful to check whether the solution which = 0 along the path changes
        """
        plt.figure()
        plt.plot(l_lc,abs(H_eigvals[:,0]),'ro')
        plt.plot(l_lc,abs(H_eigvals[:,1]),'go')
        plt.plot(l_lc,abs(H_eigvals[:,2]),'bo')
        plt.plot(l_lc,abs(H_1_Cardano_array),'r')    
        plt.plot(l_lc,abs(H_2_Cardano_array),'g')    
        plt.plot(l_lc,abs(H_3_Cardano_array),'b')    
        plt.savefig('H_' + output_figurename_suffix)
        plt.close()        

        # Commented out because this does not work properly
        # """
        # Plots Psi before and after the BCs are applied 
        # """
        # K_magnitude_entry = np.sqrt(K_R_entry**2 + K_zeta_entry**2 * entry_position[0]**2 + K_Z_entry**2)
        
        # Psi_w_entry = np.array([
        # [Psi_xx_entry,Psi_xy_entry],
        # [Psi_xy_entry,Psi_yy_entry]
        # ])
        
        # Psi_w_initial = np.array([
        #         [Psi_xx_output[0],Psi_xy_output[0]],
        #         [Psi_xy_output[0],Psi_yy_output[0]]
        #         ])
        
        # [Psi_w_entry_real_eigval_a, Psi_w_entry_real_eigval_b], Psi_w_entry_real_eigvec = np.linalg.eig(np.real(Psi_w_entry))
        # [Psi_w_entry_imag_eigval_a, Psi_w_entry_imag_eigval_b], Psi_w_entry_imag_eigvec = np.linalg.eig(np.imag(Psi_w_entry))
        # Psi_w_entry_real_eigvec_a = Psi_w_entry_real_eigvec[:,0]
        # Psi_w_entry_real_eigvec_b = Psi_w_entry_real_eigvec[:,1]
        # Psi_w_entry_imag_eigvec_a = Psi_w_entry_imag_eigvec[:,0]
        # Psi_w_entry_imag_eigvec_b = Psi_w_entry_imag_eigvec[:,1]
        
        # [Psi_w_initial_real_eigval_a, Psi_w_initial_real_eigval_b], Psi_w_initial_real_eigvec = np.linalg.eig(np.real(Psi_w_initial))
        # [Psi_w_initial_imag_eigval_a, Psi_w_initial_imag_eigval_b], Psi_w_initial_imag_eigvec = np.linalg.eig(np.imag(Psi_w_initial))
        # Psi_w_initial_real_eigvec_a = Psi_w_initial_real_eigvec[:,0]
        # Psi_w_initial_real_eigvec_b = Psi_w_initial_real_eigvec[:,1]
        # Psi_w_initial_imag_eigvec_a = Psi_w_initial_imag_eigvec[:,0]
        # Psi_w_initial_imag_eigvec_b = Psi_w_initial_imag_eigvec[:,1]
        
        # numberOfPlotPoints = 50
        # sin_array = np.sin(np.linspace(0,2*np.pi,numberOfPlotPoints))
        # cos_array = np.cos(np.linspace(0,2*np.pi,numberOfPlotPoints))
        
        # width_ellipse_entry = np.zeros([numberOfPlotPoints,2])
        # width_ellipse_initial = np.zeros([numberOfPlotPoints,2])
        # rad_curv_ellipse_entry = np.zeros([numberOfPlotPoints,2])
        # rad_curv_ellipse_initial = np.zeros([numberOfPlotPoints,2])
        # for ii in range(0,numberOfPlotPoints):
        #     width_ellipse_entry[ii,:] = np.sqrt(2/Psi_w_entry_imag_eigval_a)*Psi_w_entry_imag_eigvec_a*sin_array[ii] + np.sqrt(2/Psi_w_entry_imag_eigval_b)*Psi_w_entry_imag_eigvec_b*cos_array[ii]
        #     width_ellipse_initial[ii,:] = np.sqrt(2/Psi_w_initial_imag_eigval_a)*Psi_w_initial_imag_eigvec_a*sin_array[ii] + np.sqrt(2/Psi_w_initial_imag_eigval_b)*Psi_w_initial_imag_eigvec_b*cos_array[ii]
        
        #     rad_curv_ellipse_entry[ii,:] = (K_magnitude_entry/Psi_w_entry_real_eigval_a)*Psi_w_entry_real_eigvec_a*sin_array[ii] + (K_magnitude_entry/Psi_w_entry_real_eigval_b)*Psi_w_entry_real_eigvec_b*cos_array[ii]
        #     rad_curv_ellipse_initial[ii,:] = (K_magnitude_array[0]/Psi_w_initial_real_eigval_a)*Psi_w_initial_real_eigvec_a*sin_array[ii] + (K_magnitude_array[0]/Psi_w_initial_real_eigval_b)*Psi_w_initial_real_eigvec_b*cos_array[ii]

        
        # plt.figure()
        # plt.subplot(1,2,1)
        # plt.plot(width_ellipse_entry[:,0],width_ellipse_entry[:,1])
        # plt.plot(width_ellipse_initial[:,0],width_ellipse_initial[:,1])
        # plt.gca().set_aspect('equal', adjustable='box')
        
        # plt.subplot(1,2,2)
        # plt.plot(rad_curv_ellipse_entry[:,0],rad_curv_ellipse_entry[:,1])
        # plt.plot(rad_curv_ellipse_initial[:,0],rad_curv_ellipse_initial[:,1])
        # plt.gca().set_aspect('equal', adjustable='box')
        # plt.savefig('BC_' + output_filename_suffix)
        # plt.close()

        
        print('Figures have been saved')
    ## -------------------


    return None
