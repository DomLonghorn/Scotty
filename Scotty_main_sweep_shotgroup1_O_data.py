# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 10:44:34 2018

@author: VH Hall-Chen
Valerian Hongjie Hall-Chen
valerian@hall-chen.com


For shot 29908, the EFIT++ times are efit_times = np.linspace(0.155,0.25,20)
I want efit_times[np.arange(0,10)*2 + 1]. 160ms, 170ms, ..., 250ms
"""
from Scotty_beam_me_up import beam_me_up
from Scotty_fun_general import find_q_lab_Cartesian, find_q_lab, find_K_lab_Cartesian, find_K_lab, find_waist, find_Rayleigh_length, genray_angles_from_mirror_angles
from Scotty_fun_general import propagate_beam

from scipy import constants
import math
import numpy as np
import sys
import os

from Scotty_init_bruv import get_parameters_for_Scotty


equil_times = np.array([0.16,0.19,0.22,0.25])
# equil_times = np.array([0.210])
mirror_rotations = np.linspace(3,-6,46)
# mirror_rotations = np.linspace(3,1.2,10)
mirror_tilt = -4
launch_freqs_GHz = np.array([55.0,57.5,60.0,62.5,67.5,70.0,72.5,75.0])


for ii, equil_time in enumerate(equil_times):
    for jj, mirror_rotation in enumerate(mirror_rotations):
        for kk, launch_freq_GHz in enumerate(launch_freqs_GHz):
            args_dict, kwargs_dict = get_parameters_for_Scotty(
                                          'DBS_NSTX_MAST',
                                          launch_freq_GHz = launch_freq_GHz,
                                          mirror_rotation = mirror_rotation, # angle, in deg
                                          mirror_tilt     = mirror_tilt, # angle, in deg
                                          find_B_method   = 'EFITpp', # EFITpp, UDA_saved, UDA, torbeam
                                          find_ne_method  = 'poly3',
                                          equil_time      = equil_time,
                                          shot            = 29908,
                                          user            = 'Valerian_desktop'
                                         )
            
            if args_dict['launch_freq_GHz'] > 52.5:
                args_dict['mode_flag'] = 1
            else:
                args_dict['mode_flag'] = -1
    
            if args_dict['mode_flag'] == 1:
                mode_string = 'O'
            elif args_dict['mode_flag'] == -1:
                mode_string = 'X' 

            kwargs_dict['output_filename_suffix'] = (
                                        '_r' + f'{mirror_rotation:.1f}'
                                        '_t' + f'{mirror_tilt:.1f}'
                                      + '_f' + f'{launch_freq_GHz:.1f}'
                                      + '_'  + mode_string
                                      + '_'  + f'{equil_time*1000:.3g}' + 'ms'
                                          )      

            kwargs_dict['figure_flag'] = False
            kwargs_dict['output_path'] = 'D:\\Dropbox\\VHChen2021\\Data - Scotty\\Run 20\\'
            kwargs_dict['density_fit_parameters'] = None
            kwargs_dict['ne_data_path'] = 'D:\\Dropbox\\VHChen2021\\Data - Equilibrium\MAST\\'
            kwargs_dict['input_filename_suffix'] = '_shotgroup1_avr_' + f'{equil_time*1000:.0f}' +'ms'

            if equil_time == 0.16:
                kwargs_dict['poloidal_flux_enter'] = 1.17280433**2
            if equil_time == 0.19:
                kwargs_dict['poloidal_flux_enter'] = 1.20481476**2
            if equil_time == 0.22:
                kwargs_dict['poloidal_flux_enter'] = 1.16129185**2
            if equil_time == 0.25:
                kwargs_dict['poloidal_flux_enter'] = 1.15082068**2

            kwargs_dict['delta_R'] = -0.00001
            kwargs_dict['delta_Z'] = -0.00001
            kwargs_dict['delta_K_R'] = 0.1
            kwargs_dict['delta_K_zeta'] = 0.1
            kwargs_dict['delta_K_Z'] = 0.1
            kwargs_dict['interp_smoothing'] = 2.0
            kwargs_dict['len_tau'] = 1002
            kwargs_dict['rtol'] = 1e-4
            kwargs_dict['atol'] = 1e-7

            if ii == 0 and jj == 0 and kk == 0:
                kwargs_dict['verbose_output_flag'] = True
            else:
                kwargs_dict['verbose_output_flag'] = False  
                
                
            data_output = kwargs_dict['output_path'] + 'data_output' + kwargs_dict['output_filename_suffix'] + '.npz'
            analysis_output = kwargs_dict['output_path'] + 'analysis_output' + kwargs_dict['output_filename_suffix'] + '.npz'
            if os.path.exists(data_output) and os.path.exists(analysis_output):
                continue
            else:   
                beam_me_up(**args_dict, **kwargs_dict)                 
                
            # beam_me_up(**args_dict, **kwargs_dict)
    
