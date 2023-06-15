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
CutoffDensity = []
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
X_UHR =0.0
suffix2 = '_ColdO'
loadfile=np.load("ray_output" + suffix2 + ".npz")
ray1_R_coords = loadfile['q_R_array_ray']
ray1_Z_coords = loadfile['q_Z_array_ray']
K_R_ray1 = loadfile['K_R_array_ray']
K_Z_ray1 = loadfile['K_Z_array_ray']
K_zeta_ray1 = loadfile['K_zeta_initial']
B_R_ray1 = loadfile['B_R_ray']
B_Z_ray1 = loadfile['B_Z_ray']
B_T_ray1 = loadfile['B_T_ray']
ne_ray = loadfile['ne_ray']
H_output = loadfile['H_output']
theta_pol = loadfile['theta_0']
Temp = loadfile['Te']
loadfile.close()

B_magnitude_array = np.sqrt(B_R_ray1**2 + B_T_ray1**2 + B_Z_ray1**2)
K_magnitude_array = np.sqrt(K_R_ray1**2 + K_Z_ray1**2 + K_zeta_ray1**2)


theta_m_output = []
BigX_ray_output = []
BigY_ray_output = []
for i in range(len(B_R_ray1)):
     
    B_mag = B_magnitude_array[i]
    K_mag = K_magnitude_array[i]

    B_hat = np.array([B_R_ray1[i], B_T_ray1[i], B_Z_ray1[i]])/B_mag
    K_hat = np.array([K_R_ray1[i], K_zeta_ray1 / ray1_R_coords[i], K_Z_ray1[i]])/K_mag

    sin_theta_m_sq = (np.dot(B_hat, K_hat)) ** 2  
    theta_m = np.arcsin(np.sqrt(sin_theta_m_sq))
    theta_m_output.append(theta_m)


    ne = ne_ray[i] * 10**19
    wpe = hpf.wpe_func(ne) 
    wce = hpf.wce_func(B_mag)
    BigY_ray = wce/launch_angular_frequency
    BigX_ray = (wpe**2)/(launch_angular_frequency**2)
    print(np.round(BigX_ray + BigY_ray**2, 1))
    if np.round(BigX_ray + BigY_ray**2, 1) == 1:
         print(BigX_ray)
         print("result!")
         X_UHR = BigX_ray

    BigX_ray_output.append(BigX_ray)
    BigY_ray_output.append(BigY_ray)


theta_m_output = np.array(theta_m_output)
K_perp_array = K_magnitude_array * np.cos(theta_m)
K_para_array = K_magnitude_array * np.sin(theta_m)

N_perp_array1 = (K_perp_array * constants.c)/launch_angular_frequency
N_perp_max_index = np.argmax(N_perp_array1)
# print(N_perp_max_index)
# print(N_perp_array1[498])
# N_perp_array1 = N_perp_array1[0:N_perp_max_index]
# BigX_ray_output = BigX_ray_output[0:N_perp_max_index]

N_para_array1 = (K_para_array * constants.c)/launch_angular_frequency
# N_para_array1 = N_para_array1[0:N_perp_max_index]
theta_pol = np.round(theta_pol, 4)
# X_UHR = 0.0
plt.figure(dpi=120,figsize=(12, 3))
plt.suptitle(r"$\theta_{pol}$ = " + str(round(theta_pol, 4)))
plt.subplot(141)
plt.plot(BigX_ray_output, N_perp_array1, label='O mode launch')
plt.xlabel("X")
plt.ylabel(r'$N_{\perp}$')
plt.title( "Density Profile of O-mode launch")
plt.axvline(x=(X_UHR), color='red', linestyle='--', label=r'$X_{UHR}$')
plt.scatter(BigX_ray_output[0], N_perp_array1[0], color='green', label='Launch Position')
plt.legend()

plt.subplot(142)
plt.plot(BigX_ray_output, N_para_array1, label='O mode launch')
plt.xlabel("X")
plt.ylabel(r'$N_{\parallel}$')
plt.axvline(x=(X_UHR), color='red', linestyle='--', label=r'$X_{UHR}$')
plt.scatter(BigX_ray_output[0], N_para_array1[0], color='green', label='Launch Position')
plt.legend()

suffix3 = '_OModeResonance1'
loadfile=np.load("ray_output" + suffix3 + ".npz")
ray2_R_coords = loadfile['q_R_array_ray']
ray2_Z_coords = loadfile['q_Z_array_ray']



print(ray1_R_coords.shape)
## For plotting how the beam propagates from launch to entry
launch_position_X, launch_position_Y, launch_position_Z = find_q_lab_Cartesian(
    launch_position
)
entry_position_X, entry_position_Y, entry_position_Z = find_q_lab_Cartesian(
    np.array([q_R_array[0], q_zeta_array[0], q_Z_array[0]])
)

numberOfDataPoints = np.size(q_R_array)
out_index = numberOfDataPoints

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



def equilcontourplot(GridToPlot, title, filename="NONE", showfig=False, savefig=True):
    if filename == "NONE":
        filename = title
    plt.title(r'$B_Z$')
    CS = plt.contourf(
    data_R_coord,
    data_Z_coord,
    np.transpose(GridToPlot.reshape(len(data_R_coord), len(data_Z_coord))),
    levels = 100,
    cmap="inferno",
    )
    plt.xlim(data_R_coord[0], data_R_coord[-1])
    plt.ylim(data_Z_coord[0], data_Z_coord[-1])
    plt.xlabel("R / m")
    plt.ylabel("Z / m")
    plt.colorbar()
    if showfig:
        plt.show()
    if savefig:
        plt.savefig(filename + ".png")
    plt.close()

# equilcontourplot(B_R_grid, "B_R grid")
# equilcontourplot(B_T_grid, "B_T_grid")
# equilcontourplot(B_Z_grid, "B_Z_grid")
# equilcontourplot(Electron_Density_Grid, "Electron Density")

# plt.title("Poloidal Plane, Te = 5.0 keV, " r'$\theta_{pol}$ = -$\pi - 0.2$')
plt.subplot(143)
plt.title("Poloidal Plane")
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

# plt.plot(q_R_array[:out_index], q_Z_array[:out_index], "k", label="Ray Trajectory (O mode)")
plt.plot(ray1_R_coords, ray1_Z_coords, "b", label="O Mode Launch")
# plt.plot(ray2_R_coords, ray2_Z_coords, "k", label="O Mode Launch")
plt.plot(ray1_R_coords[0], ray1_Z_coords[0], "rx", label="Launch Position")
# plt.plot([launch_position[0], q_R_array[0]], [launch_position[2], q_Z_array[0]], ":k")
plt.xlim(data_R_coord[0], data_R_coord[-1])
plt.ylim(data_Z_coord[0], data_Z_coord[-1])
plt.xlabel("R / m")
plt.ylabel("Z / m")
plt.legend()


# H_output = H_output[0:N_perp_max_index]
plt.subplot(144)
plt.title('$\log_{10}(|\mathcal{M}|)$')
plt.plot(BigX_ray_output, np.log10(abs(H_output)), label='O mode launch')
plt.scatter(BigX_ray_output[0], np.log10(abs(H_output[0])), color='green', label='Launch Position')
plt.show()
# plt.savefig("Ray_Propagation.jpg")


# plt.figure
# plt.title("Poloidal Plane")
# contour_levels = np.linspace(0, 1, 11)
# CS = plt.contour(
#     data_R_coord,
#     data_Z_coord,
#     np.transpose(poloidalFlux_grid.reshape(len(data_R_coord), len(data_Z_coord))),
#     contour_levels,
#     vmin=0,
#     vmax=1.2,
#     cmap="inferno",
# )
# plt.clabel(
#     CS, inline=True, fontsize=10, inline_spacing=-5, fmt="%1.1f", use_clabeltext=True
# )  # Labels the flux surfaces
# if UHR_R_coords is not None and UHR_Z_coords is not None:
#         plt.scatter(UHR_R_coords, UHR_Z_coords, c='green', s=4, label=r'$X_{UHR}$')
# if CritDensity_R_coords is not None and CritDensity_Z_coords is not None:
#         plt.scatter(CritDensity_R_coords, CritDensity_Z_coords, c='red', s=4, label=r'$X_{crit}$')
# if Xcutoff_R_coords is not None and Xcutoff_Z_coords is not None:
#         plt.scatter(Xcutoff_R_coords, Xcutoff_Z_coords, c='blue', s=4, label='x mode Cutoff')

# # plt.plot(q_R_array[:out_index], q_Z_array[:out_index], "k", label="Ray Trajectory (O mode)")
# plt.plot(ray1_R_coords, ray1_Z_coords, "b", label="X Mode Launch")
# # plt.plot(ray2_R_coords, ray2_Z_coords, "k", label="O Mode Launch")
# plt.plot(ray1_R_coords[0], ray1_Z_coords[0], "rx", label="Launch Position")
# # plt.plot([launch_position[0], q_R_array[0]], [launch_position[2], q_Z_array[0]], ":k")
# plt.xlim(data_R_coord[0], data_R_coord[-1])
# plt.ylim(data_Z_coord[0], data_Z_coord[-1])
# plt.xlabel("R / m")
# plt.ylabel("Z / m")
# plt.legend()
# plt.show()