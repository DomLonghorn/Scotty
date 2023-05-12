import numpy as np
import matplotlib.pyplot as plt
from Scotty_fun_general import (
    find_waist,
    find_distance_from_waist,
    find_q_lab_Cartesian,
    find_nearest,
    contract_special,
)
from Scotty_fun_general import (
    find_normalised_plasma_freq,
    find_normalised_gyro_freq,
    make_unit_vector_from_cross_product,
    find_vec_lab_Cartesian,
)

import math
from scipy import constants, integrate
import sys

suffix = "_ColdOMode10.20"

# Loads in the output variables needed
loadfile = np.load("data_output" + suffix + ".npz")
B_magnitude = loadfile["B_magnitude"]
B_R_output = loadfile["B_R_output"]
B_T_output = loadfile["B_T_output"]
B_Z_output = loadfile["B_Z_output"]
K_Z_array=loadfile["K_Z_array"]
K_R_array=loadfile["K_R_array"]

K_zeta_initial=loadfile["K_zeta_initial"]
b_hat_output = loadfile["b_hat_output"]
x_hat_output = loadfile["x_hat_output"]
q_R_array = loadfile["q_R_array"]
q_zeta_array = loadfile["q_zeta_array"]
q_Z_array = loadfile["q_Z_array"]
loadfile.close()

loadfile = np.load('analysis_output' + suffix + '.npz')
theta_m_array=loadfile["theta_m_output"]
loadfile.close()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    print(idx)
    return array[idx]
#Index I want is 1

print(find_nearest(q_R_array, 1.646))

B_R = B_R_output[1]
B_T = B_T_output[1]
B_Z = B_Z_output[1]
K_R = K_R_array[1]
K_Z = K_Z_array[1]
K_zeta = K_zeta_initial
q_R = q_R_array[1]

K_magnitude = np.sqrt(K_Z**2+K_R**2+K_zeta_initial**2)
B_Total = np.sqrt(B_R**2 + B_T**2 + B_Z**2)

b_hat = np.array([B_R, B_T, B_Z]) / B_Total #k_para as defined in the Booker fornulation
K_hat = np.array([K_R, K_zeta / q_R, K_Z]) / K_magnitude #k_perp

sin_theta_m_sq = (np.dot(b_hat, K_hat)) ** 2  
theta_m = np.arcsin(np.sqrt(sin_theta_m_sq))

K_perp = K_magnitude * np.cos(theta_m)
K_para = K_magnitude * np.sin(theta_m)

print("K_perp", K_perp)
print("K_para", K_para)
print("q_R", q_R_array[1])
print("q_Z", q_Z_array[1])
print("q_zeta", q_zeta_array[1])

K = np.array([K_R, K_zeta, K_Z])
B = np.array([B_R, B_T, B_Z])

