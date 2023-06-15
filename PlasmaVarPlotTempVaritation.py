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
import hot_plasma_fun as hpf
import math
from scipy import constants, integrate
import sys

suffix = "_ReferenceColdOMode0.20"


# Loads in the output variables needed
loadfile = np.load("data_output" + suffix + ".npz")
B_magnitude = loadfile["B_magnitude"]
b_hat_output = loadfile["b_hat_output"]
x_hat_output = loadfile["x_hat_output"]
q_R_array = loadfile["q_R_array"]
q_zeta_array = loadfile["q_zeta_array"]
q_Z_array = loadfile["q_Z_array"]
normalised_plasma_freqs=loadfile["normalised_plasma_freqs"],
normalised_gyro_freqs=["normalised_gyro_freqs"],
loadfile.close()

# Loads in input conditions
loadfile = np.load("data_input" + suffix + ".npz")
launch_position = loadfile["launch_position"]
data_R_coord = loadfile["data_R_coord"]
data_Z_coord = loadfile["data_Z_coord"]
poloidalFlux_grid = loadfile["poloidalFlux_grid"]
# ne_data_density_array = loadfile["ne_data_density_array"]
loadfile.close()

loadfile = np.load("analysis_output" + suffix + ".npz")
e_hat_output = loadfile['e_hat_output']
theta_m_output = loadfile['theta_m_output']

# Loads in equilibrium data
loadfile = np.load("Equilibrium_Output" + suffix + ".npz")
B_R_grid = loadfile["B_R_grid"]
B_T_grid = loadfile["B_T_grid"]
B_Z_grid = loadfile["B_Z_grid"]
Electron_Density = loadfile['Electron_Density']
Electron_Density_Grid = loadfile['Electron_Density_Grid']
loadfile.close()

launch_freq_GHz = 55.0
launch_angular_frequency = 2 * math.pi * 10.0**9 * launch_freq_GHz
# print(q_R_array)

resultgrid = np.zeros_like(B_R_grid)
UHR_R_coords = []
UHR_Z_coords = []
CritDensity_R_coords = []
CritDensity_Z_coords = []
Xcutoff_R_coords = []
Xcutoff_Z_coords = []
for i in range(len(B_R_grid)):
    for j in range(len(B_R_grid)):
        B_R = B_R_grid[i][j]
        B_T = B_T_grid[i][j]
        B_Z = B_Z_grid[i][j]
        B_magnitude = np.sqrt(B_R**2 + B_T**2 + B_Z**2)
        Electron_Density = Electron_Density_Grid[i][j]
        Electron_Density = Electron_Density* 10**19

        wce = hpf.wce_func(B_magnitude)
        wpe = hpf.wpe_func(Electron_Density)

        bigX = (wpe**2)/(launch_angular_frequency**2)
        bigY = wce/launch_angular_frequency

        UHRcalc = bigX + bigY**2
        UHRcalc = round(UHRcalc, 2)
        if UHRcalc == 1.0:
            UHR_R_coords.append(data_R_coord[i])
            UHR_Z_coords.append(data_Z_coord[j])
        if round(bigX,2) == 1:
            CritDensity_R_coords.append(data_R_coord[i])
            CritDensity_Z_coords.append(data_Z_coord[j])
        R = 1 - ((bigX)/(1-bigY))
        R = round(R, 3)
        R = abs(R)
        if R <= 0.01:
            Xcutoff_R_coords.append(data_R_coord[i])
            Xcutoff_Z_coords.append(data_Z_coord[j])



# suffix2 = "_XTest11"
# suffix2 = "_CritLaunch2"

# loadfile=np.load("ray_output" + suffix2 + ".npz")
# ray1_R_coords = loadfile['q_R_array_ray']
# ray1_Z_coords = loadfile['q_Z_array_ray']
# K_R_ray1 = loadfile['K_R_array_ray']
# K_Z_ray1 = loadfile['K_Z_array_ray']
# K_zeta_ray1 = loadfile['K_zeta_initial']
# B_R_ray1 = loadfile['B_R_ray']
# B_Z_ray1 = loadfile['B_Z_ray']
# B_T_ray1 = loadfile['B_T_ray']
# ne_ray = loadfile['ne_ray']

# B_magnitude_array = np.sqrt(B_R_ray1**2 + B_T_ray1**2 + B_Z_ray1**2)
# K_magnitude_array = np.sqrt(K_R_ray1**2 + K_Z_ray1**2 + K_zeta_ray1**2)


# theta_m_output = []
# BigX_ray_output = []
# BigY_ray_output = []


RCoords = []
ZCoords = []
temps = []
for i in range(2, 11, 2):
    print(i)
    suffix = '_' + str(i) + 'keVFlag1'
    loadfile=np.load("ray_output" + suffix + ".npz")
    ray_R_coords = loadfile['q_R_array_ray']
    ray_Z_coords = loadfile['q_Z_array_ray']
    RCoords.append(ray_R_coords)
    ZCoords.append(ray_Z_coords)
    if i == 0:
        temps.append("Cold Dispersion Launch")
    else:
         temps.append(str(i)) 
    
ColdDispSuffix = '_ColdLaunch'
loadfile=np.load("ray_output" + ColdDispSuffix + ".npz")
ray_R_coords = loadfile['q_R_array_ray']
ray_Z_coords = loadfile['q_Z_array_ray']

# def checkpolarisation(e_hat):
#     val = np.real(np.dot(e_hat, np.array([0,0,1])))
    
#     val = round(val, 4)
#     if val == 0.0:
#         polarisation = "X-Mode"
#     elif val > 0:
#         polarisation = "O-Mode"
#     else:
#         polarisation = "Error: Polarisation Unknown"
#     return polarisation




plt.figure(figsize=(10, 8))
plt.title("Variation of Te for an O-mode launch, " r'$\theta_{pol}$ = $- 0.2$')
contour_levels = np.linspace(0, 1, 11)
CS = plt.contour(
    data_R_coord,
    data_Z_coord,
    np.transpose(poloidalFlux_grid.reshape(len(data_R_coord), len(data_Z_coord))),
    contour_levels,
    vmin=0,
    vmax=1.2,
    cmap="inferno",
)
plt.clabel(
    CS, inline=True, fontsize=10, inline_spacing=-5, fmt="%1.1f", use_clabeltext=True
)  # Labels the flux surfaces
if UHR_R_coords is not None and UHR_Z_coords is not None:
        plt.scatter(UHR_R_coords, UHR_Z_coords, c='green', s=4, label=r'$X_{UHR}$')
if CritDensity_R_coords is not None and CritDensity_Z_coords is not None:
        plt.scatter(CritDensity_R_coords, CritDensity_Z_coords, c='red', s=4, label=r'$X_{crit}$')
if Xcutoff_R_coords is not None and Xcutoff_Z_coords is not None:
        plt.scatter(Xcutoff_R_coords, Xcutoff_Z_coords, c='blue', s=4, label='x mode Cutoff')

plt.plot(ray_R_coords, ray_Z_coords, label='Cold Dispersion Launch')
for i in range(len(RCoords)):
    plt.plot(RCoords[i], ZCoords[i], label=(temps[i]) + '.0 keV')
plt.plot(RCoords[0][0], ZCoords[0][0], 'kx', label='Launch Point')
plt.xlim(data_R_coord[0], data_R_coord[-1])
plt.ylim(data_Z_coord[0], data_Z_coord[-1])
plt.xlabel("R / m")
plt.ylabel("Z / m")
plt.legend()
plt.show()
# plt.savefig("Ray_Propagation.jpg")
