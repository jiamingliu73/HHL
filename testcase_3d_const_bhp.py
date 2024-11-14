import sys
import numpy as np
import pandas as pd
import collections


def main():

    # input parameters (18)
    phi = 0.18  # porosity
    C_o = 3.5E-6  # oil compressibility, psi-1
    del_x = 1000  # space interval in i, ft
    del_y = 1000  # space interval in j, ft
    del_z = 75  # space interval in k, ft
    numcells_x = 3  # number of grid in x direction
    numcells_y = 3  # number of grid in y direction
    numcells_z = 3  # number of grid in z direction
    k_xy = 15  # permeability in x and y directions, md
    k_z = 15  # permeability in z direction, md
    mu_o = 10  # viscosity, cp
    del_t = 15  # size of time step, days
    t_steps = 5  # number of time steps
    p_i = 6000  # initial pressure, psi
    p_wf = 4000  # bottom hole pressure, psi
    well_cells = [[2, 2, 2]]  # coordinates of well cell
    re = 50  # drainage radius, ft
    rw = 0.25  # wellbore radius, ft

    # parameters that can be assumed constant
    beta_c = 1.127  # unit conversion factor
    B_o = 1  # FVF
    B_o_o = 1  # FVF ref.
    skin = 0  # skin factor

    Q = []
    # transmissibility in x, y, and z directions
    T_ox = beta_c * del_y * del_z * k_xy / (mu_o * B_o * del_x) / 1000
    T_oy = beta_c * del_x * del_z * k_xy / (mu_o * B_o * del_y) / 1000
    T_oz = beta_c * del_x * del_y * k_z / (mu_o * B_o * del_z) / 1000

    total_t = t_steps * del_t  # total time, days
    V_b = del_x * del_y * del_z  # grid volume
    del_t_const = (V_b * phi * C_o) / (0.5615 * mu_o * B_o_o * del_t)  # pressure term constant dt

    X_matrix = np.zeros(shape=(numcells_x * numcells_y * numcells_z, numcells_x * numcells_y * numcells_z))
    b_vector = np.zeros(shape=(1, numcells_x * numcells_y * numcells_z))
    pressure = np.ones(shape=(1, numcells_x * numcells_y * numcells_z))
    pressure = p_i * pressure.T

    p_coef_minus_i = 0
    p_coef_minus_j = 0
    p_coef_minus_k = 0
    # p_coef = 0
    p_coef_plus_i = 0
    p_coef_plus_j = 0
    p_coef_plus_k = 0

    # creating linear system
    for t in range(t_steps):
        flow_rate = []
        well_num = 0
        a = -1
        for k in range(0, numcells_z):
            for j in range(0, numcells_y):
                for i in range(0, numcells_x):
                    a = a + 1
                    if i - 1 < 0:
                        p_coef_minus_i = 0
                    else:
                        p_coef_minus_i = T_ox
                        if k * (numcells_x * numcells_y) + j * numcells_y + (i - 1) >= 0:
                            X_matrix[k * (numcells_x * numcells_y) + j * numcells_y + i][
                                k * (numcells_x * numcells_y) + j * numcells_x + (i - 1)] = p_coef_minus_i

                    if i + 1 < numcells_x:
                        p_coef_plus_i = T_ox
                        if k * (numcells_x * numcells_y) + j * numcells_y + (i + 1) < numcells_x * numcells_y * numcells_z:
                            X_matrix[k * (numcells_x * numcells_y) + j * numcells_y + i][
                                k * (numcells_x * numcells_y) + j * numcells_x + (i + 1)] = p_coef_plus_i
                    else:
                        p_coef_plus_i = 0

                    if j - 1 < 0:
                        p_coef_minus_j = 0
                    else:
                        p_coef_minus_j = T_oy
                        if k * (numcells_x * numcells_y) + j * numcells_y + (i - 2) >= 0:
                            X_matrix[k * (numcells_x * numcells_y) + j * numcells_y + i][
                                k * (numcells_x * numcells_y) + j * numcells_y + (i - 2)] = p_coef_minus_j

                    if j + 1 < numcells_y:
                        p_coef_plus_j = T_oy
                        if k * (numcells_x * numcells_y) + j * numcells_y + (i + 2) < numcells_x * numcells_y * numcells_z:
                            X_matrix[k * (numcells_x * numcells_y) + j * numcells_y + i][
                                k * (numcells_x * numcells_y) + j * numcells_y + (i + 2)] = p_coef_plus_j
                    else:
                        p_coef_plus_j = 0

                    if k - 1 < 0:
                        p_coef_minus_k = 0
                    else:
                        p_coef_minus_k = T_oz
                        if (k - 1) * (numcells_x * numcells_y) + j * numcells_y + i >= 0:
                            X_matrix[k * (numcells_x * numcells_y) + j * numcells_y + i][
                                (k - 1) * (numcells_x * numcells_y) + j * numcells_x + i] = p_coef_minus_k

                    if k + 1 < numcells_z:
                        p_coef_plus_k = T_oz
                        if (k + 1) * (
                                numcells_x * numcells_y) + j * numcells_y + i < numcells_x * numcells_y * numcells_z:
                            X_matrix[k * (numcells_x * numcells_y) + j * numcells_y + i][
                                (k + 1) * (numcells_x * numcells_y) + j * numcells_x + i] = p_coef_plus_k
                    else:
                        p_coef_plus_k = 0

                    X_matrix[k * (numcells_x * numcells_y) + j * numcells_y + i][
                        k * (numcells_x * numcells_y) + j * numcells_x + i] = -(
                            del_t_const + p_coef_minus_i + p_coef_plus_i + p_coef_minus_j + p_coef_plus_j + p_coef_minus_k + p_coef_plus_k)

                    P = pressure[a]
                    b_vector[0][k * (numcells_x * numcells_y) + j * numcells_y + i] = -(del_t_const * P)

                    if [k, j, i] in well_cells:
                        # b_vector[0][k * (numcells_x * numcells_y) + j * numcells_y + i] = b_vector[0][k * (
                        #         numcells_x * numcells_y) + j * numcells_y + i] - (well_rate[well_num])
                        # J = (2 * np.pi * k_xy * del_z) / (mu_o * B_o * (np.log(re / rw) + skin))
                        J = (0.0078 * k_xy * del_z) / (mu_o * B_o * (np.log(re / rw) + 0.5 + skin))
                        if P < p_wf:
                            Qq = 0
                        else:
                            Qq = - J * (P - p_wf)
                        X_matrix[a][a] = X_matrix[a][a] - J
                        b_vector[0][a] = b_vector[0][a] - (J * p_wf)
                        well_num += 1
                        flow_rate.append(Qq)
                        rate = sum(flow_rate)

        pressures = np.linalg.solve(X_matrix, b_vector.T)
        pressure = pressures
        Q.append(rate.tolist())

        production = - del_t * np.sum(Q)

        print(X_matrix)

        print(b_vector)

        print(pressures)

        print(Q)  # flow rate at each time step

        print('Total production is', production, 'in', total_t, 'days')

        data1 = pd.DataFrame(X_matrix)
        data2 = pd.DataFrame(b_vector.T)
        data3 = pd.DataFrame(pressures)
        data4 = pd.DataFrame(Q)

        with pd.ExcelWriter('3dresults.xlsx') as writer:
            data1.to_excel(writer, sheet_name="X_matrix")
            data2.to_excel(writer, sheet_name="b_vector")
            data3.to_excel(writer, sheet_name="Pressures")
            data4.to_excel(writer, sheet_name="Flow rate")


if __name__ == '__main__':
    main()
