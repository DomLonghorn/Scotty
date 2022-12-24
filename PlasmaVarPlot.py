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

suffix = "_t"

# Loads in the output variables needed
loadfile = np.load("data_output" + suffix + ".npz")
electron_density_output = loadfile["electron_density_output"]
B_magnitude = loadfile["B_magnitude"]
b_hat_output = loadfile["b_hat_output"]
x_hat_output = loadfile["x_hat_output"]
q_R_array = loadfile["q_R_array"]
q_zeta_array = loadfile["q_zeta_array"]
q_Z_array = loadfile["q_Z_array"]
loadfile.close()

# Loads in input conditions
loadfile = np.load("data_input" + suffix + ".npz")
launch_position = loadfile["launch_position"]
data_R_coord = loadfile["data_R_coord"]
data_Z_coord = loadfile["data_Z_coord"]
poloidalFlux_grid = loadfile["poloidalFlux_grid"]
# ne_data_density_array = loadfile["ne_data_density_array"]
loadfile.close()

# Loads in equilibrium data
loadfile = np.load("EqulibriumQuantities" + suffix + ".npz")
B_R_grid = loadfile["B_R_grid"]
B_T_grid = loadfile["B_T_grid"]

# print(B_Z)
## For plotting how the beam propagates from launch to entry
launch_position_X, launch_position_Y, launch_position_Z = find_q_lab_Cartesian(
    launch_position
)
entry_position_X, entry_position_Y, entry_position_Z = find_q_lab_Cartesian(
    np.array([q_R_array[0], q_zeta_array[0], q_Z_array[0]])
)

numberOfDataPoints = np.size(q_R_array)
out_index = numberOfDataPoints

print(B_R_grid)

plt.title("B_T_grid")
CS = plt.contour(
    data_R_coord,
    data_Z_coord,
    np.transpose(B_T_grid.reshape(len(data_R_coord), len(data_Z_coord))),
    cmap="inferno",
)
plt.xlim(data_R_coord[0], data_R_coord[-1])
plt.ylim(data_Z_coord[0], data_Z_coord[-1])
plt.xlabel("R / m")
plt.ylabel("Z / m")
plt.show()


# plt.title("Poloidal Plane")
# contour_levels = np.linspace(0, 1, 11)
# CS = plt.contour(
#     data_R_coord,
#     data_Z_coord,
#     np.transpose(equldensityarray.reshape(len(data_R_coord), len(data_Z_coord))),
#     contour_levels,
#     vmin=0,
#     vmax=1.2,
#     cmap="inferno",
# )
# plt.clabel(
#     CS, inline=True, fontsize=10, inline_spacing=-5, fmt="%1.1f", use_clabeltext=True
# )  # Labels the flux surfaces

# # plt.plot(q_R_array[:out_index], q_Z_array[:out_index], "k")
# # plt.plot([launch_position[0], q_R_array[0]], [launch_position[2], q_Z_array[0]], ":k")
# plt.xlim(data_R_coord[0], data_R_coord[-1])
# plt.ylim(data_Z_coord[0], data_Z_coord[-1])
# plt.xlabel("R / m")
# plt.ylabel("Z / m")
# plt.savefig("Ray_Propagation.jpg")
