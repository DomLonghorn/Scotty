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
from Scotty_fun_general import (
    find_q_lab_Cartesian,
    find_q_lab,
    find_K_lab_Cartesian,
    find_K_lab,
    find_waist,
    find_Rayleigh_length,
    genray_angles_from_mirror_angles,
)
from Scotty_fun_general import propagate_beam

from scipy import constants
import math
import numpy as np
import sys

from Scotty_init_bruv import get_parameters_for_Scotty




B_p_a = 0.20
args_dict, kwargs_dict = get_parameters_for_Scotty(
                                'DBS_synthetic'
                                )
args_dict['mode_flag'] = -1
args_dict["launch_position"] = np.array(
        [1.3, 0, -0.0]
    )
# args_dict["launch_position"] = np.array(
#         [2.1, 0, 0.0]
#     )
args_dict["poloidal_launch_angle_Torbeam"] = 1.0
args_dict["toroidal_launch_angle_Torbeam"] = 0.0
kwargs_dict['B_p_a'] = B_p_a
kwargs_dict["vacuum_propagation_flag"] = False
kwargs_dict["vacuumLaunch_flag"] = False
kwargs_dict["Psi_BC_flag"] = True
kwargs_dict['output_filename_suffix'] = (
                            '_X2'
                                )      

beam_me_up(**args_dict, **kwargs_dict)

# for tor_launch_angle in tor_launch_angles:
#     args_dict, kwargs_dict = get_parameters_for_Scotty(
#                                   'DBS_synthetic_Hot'
#                                   )
    
#     args_dict['toroidal_launch_angle_Torbeam'] = tor_launch_angle
#     kwargs_dict['B_p_a'] = 0.0
#     kwargs_dict['output_filename_suffix'] = (
#                                 '_t' + f'{tor_launch_angle:.2f}'
#                                   )      
  
#     beam_me_up(**args_dict, **kwargs_dict)
