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

suffix = "_t3.00"

# Loads in the output variables needed
loadfile = np.load("data_output" + suffix + ".npz")
electron_density_output = loadfile["electron_density_output"]
B_magnitude = loadfile["B_magnitude"]
b_hat_output = loadfile["b_hat_output"]
x_hat_output = loadfile["x_hat_output"]
loadfile.close()

# Loads in input conditions
loadfile = np.load("data_input" + suffix + ".npz")
data_R_coord = loadfile["data_R_coord"]
data_Z_coord = loadfile["data_Z_coord"]
poloidalFlux_grid = loadfile["poloidalFlux_grid"]
loadfile.close()

# Creates the components of the B field in the R, T and Z planes
B_R = B_magnitude * b_hat_output[:, 0]
B_T = B_magnitude * b_hat_output[:, 1]
B_Z = B_magnitude * b_hat_output[:, 2]

# print(len(x_hat_output))
# print(len(electron_density_output))
# print(len(B_R))
# print(len(data_Z_coord))


plt.figure(figsize=(5, 5))
plt.title("Poloidal Plane")
contour_levels = np.linspace(0, 1, 11)
CS = plt.contour(
    data_R_coord,
    data_Z_coord,
    B_T,
    contour_levels,
    vmin=0,
    vmax=1.2,
    cmap="inferno",
)
plt.clabel(
    CS, inline=True, fontsize=10, inline_spacing=-5, fmt="%1.1f", use_clabeltext=True
)
plt.xlim(data_R_coord[0], data_R_coord[-1])
plt.ylim(data_Z_coord[0], data_Z_coord[-1])

plt.xlabel("R / m")
plt.ylabel("Z / m")
