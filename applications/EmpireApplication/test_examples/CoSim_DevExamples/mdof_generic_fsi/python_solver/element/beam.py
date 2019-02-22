# Beam element

# import python modules
from scipy. optimize import minimize
from scipy import linalg
import matplotlib.pyplot as plt
import numpy as np


class Beam():

    def __init__(self, properties, rdof, target_freq):

        self.properties = properties
        self.rdof = rdof

        self.rdof_spring = [1, 0]
        self.rdof_beam = [1, 0]
        self.EI = self.optimize(target_freq)
        self.K, self.M, self.B, self.K_big, self.M_big, self.B_big = self.beam(self.EI)

        self.eig_vals, self.eig_vecs_raw, self.eig_freq, self.eig_per = self.eigen_value(
            self.K, self.M)

        self.B, self.B_big = self.damping(2, 0.01)

    def beam(self, EI):

        elem_number = self.properties.levels
        elem_length = self.properties.height / elem_number

        rho = self.properties.density
        material_area = self.properties.length * self.properties.width

        # mass matrix
        # beam
        def mass_matrix():
            # material constant
            material_constant = rho * material_area * elem_length / 420

            # element matrix
            m = np.array(
                [[156, 22 * elem_length, 54, -13 * elem_length],
                 [22 * elem_length, 4 * elem_length * elem_length,
                  13 * elem_length, -3 * elem_length * elem_length],
                 [54, 13 * elem_length, 156, -22 * elem_length],
                 [-13 * elem_length, -3 * elem_length * elem_length, -22 * elem_length, 4 * elem_length * elem_length]])

            # global matrix
            M = np.zeros((2 * elem_number + 2, 2 * elem_number + 2))

            for i in range(elem_number):
                M_temp = np.zeros((2 * elem_number + 2, 2 * elem_number + 2))
                M_temp[2 * i:2 * i + 4, 2 * i:2 * i + 4] = m
                M += material_constant * M_temp

            M_big = M

            # remove the fixed degrees of freedom
            for dof in self.rdof_beam:
                for i in [0, 1]:
                    M = np.delete(M, dof, axis=i)

            # return mass matrix
            return M, M_big

        def stiffness_matrix():
            # stiffness constant
            stiffness_constant = EI / pow(elem_length, 3)

            # element matrix
            k = np.array(
                [[12, 6 * elem_length, -12, 6 * elem_length],
                 [6 * elem_length, 4 * elem_length * elem_length,
                  -6 * elem_length, 2 * elem_length * elem_length],
                    [-12, -6 * elem_length, 12, -6 * elem_length],
                    [6 * elem_length, 2 * elem_length * elem_length, -6 * elem_length, 4 * elem_length * elem_length]])

            # global matrix
            K = np.zeros((2 * elem_number + 2, 2 * elem_number + 2))

            for i in range(elem_number):
                K_temp = np.zeros((2 * elem_number + 2, 2 * elem_number + 2))
                K_temp[2 * i:2 * i + 4, 2 * i:2 * i + 4] = k
                K += stiffness_constant * K_temp

            K_big = K

            # remove the fixed degrees of freedom
            for dof in self.rdof_beam:
                for i in [0, 1]:
                    K = np.delete(K, dof, axis=i)

            # return stiffness matrix
            return K, K_big

        # return matrices
        K, K_big = stiffness_matrix()
        M, M_big = mass_matrix()
        B = np.zeros(K.shape)
        B_big = np.zeros(K.shape)

        return [K, M, B, K_big, M_big, B_big]

    def eigen_value(self, K, M):

        # raw eigenvalues
        eig_vals_raw, eig_vecs_raw = linalg.eigh(K, M)
        # real eigenvalues
        eig_vals = np.sqrt(np.real(eig_vals_raw))
        eig_freq = eig_vals / 2 / np.pi  # in Hz
        eig_per = 1. / eig_freq  # in second

        return[eig_vals, eig_vecs_raw, eig_freq, eig_per]

    def optimize(self, target_freq):

        def objective(param):
            EI = param
            K = self.beam(EI)[0]
            M = self.beam(EI)[1]
            return (self.eigen_value(K, M)[2][0] - target_freq) ** 2

        init = 1000
        result = minimize(objective, init, method='Nelder-Mead')
        return result.x

    def damping(self, nr_modes, damping):

        # Damping ratios
        xi = [damping for x in range(nr_modes)]
        # xi = [damping] * nr_modes

        # LHS assambly
        LHS = np.zeros((nr_modes, nr_modes))  # nModes x nModes matrix

        for i in range(nr_modes):
            for j in range(nr_modes):
                LHS[i][j] = pow(self.eig_vals[i], (2 * j) - 1)

        # RHS assambly
        RHS = np.dot(2, xi)

        # Solve for the for the unknown damping factors
        a = np.linalg.solve(LHS, RHS)
        # Setup Caughey damping matrix
        # CauC = a[0] * self.M + a[1] * self.K + a[2] * \
        #     np.dot(np.dot(self.K, linalg.inv(self.M)), self.K)

        C = np.zeros(len(self.K))

        for i in range(nr_modes):
            C = np.add(C, np.dot(
                self.M, (a[i] * np.asarray(np.matrix(np.dot(linalg.inv(self.M), self.K)) ** i))))
        
        C_big = np.zeros(self.K_big.shape)

        for i in range(nr_modes):
            C_big = np.add(C_big, np.dot(
                self.M_big, (a[i] * np.asarray(np.matrix(np.dot(linalg.inv(self.M_big), self.K_big)) ** i))))
        # x = np.zeros(100)
        # for i in range(100):
        #     x[i] = (self.eig_vals[nr_modes] / 100) * i
        # cita = np.zeros(100)
        # for i in range(100):
        #     cita[i] = 0
        #     for j in range(nr_modes):
        #         cita[i] += 0.5 * a[j] * pow(x[i], (2 * j) - 1)

        # plt.plot(self.eig_vals[0:nr_modes], xi, marker='+')
        # plt.plot(x, cita)
        # plt.ylim([-100, 100])
        # plt.show()

        return C, C_big

    def eigen_value_load(self, mode):

        num_elems = self.properties.levels
        elem_length = self.properties.height / num_elems

        # Sorted inteces
        freq_ind = np.argsort(self.eig_freq)

        # Mass normalization
        [n_row, m_col] = self.eig_vecs_raw.shape
        eig_vecs_norm = np.zeros([n_row, m_col])

        gen_mass_raw = np.zeros(m_col)
        gen_mass_norm = np.zeros(m_col)

        for i in range(len(self.eig_vals)):
            gen_mass_raw[i] = (np.transpose(self.eig_vecs_raw[:, i])).dot(
                self.M).dot(self.eig_vecs_raw[:, i])

            gen_mass_factor = np.sqrt(gen_mass_raw[i])

            eig_vecs_norm[:, i] = self.eig_vecs_raw[:, i] / gen_mass_factor

            gen_mass_norm[i] = (np.transpose(eig_vecs_norm[:, i])).dot(
                self.M).dot(eig_vecs_norm[:, i])

        eigen_form = eig_vecs_norm[:, freq_ind[int(mode)]]
        nodal_disp = eigen_form

        def plot_modes():
            # Extend eigen modes for plotting
            extended_eig_vecs = np.zeros(
                [self.properties.levels + 1, self.properties.levels])
            extended_eig_vecs[1:, :] = eig_vecs_norm[::2, ::2]

            # Plot eigenvalues and shapes
            fig = plt.figure(1)

            # Number of modes to iterate
            modes = 1
            for i in range(modes):
                plt.plot(extended_eig_vecs[:, freq_ind[i]], Z, marker='o', label="mode " + str(i + 1) + ", eigFreq: " + str(np.round(self.eig_freq[freq_ind[i]], 3))
                         + " [Hz], T: " + str(np.round(self.eig_per[freq_ind[i]], 3)) + " [s], omega: " + str(np.round(self.eig_vals[freq_ind[i]], 3)) + " [rad/s]")
                plt.plot(delta, Z, marker='o')

            plt.title("Eigenmodes - mass normalized eigenvectors")
            plt.xlabel("Eigenmode magnitude")
            plt.ylabel("Height [m]")
            plt.legend(loc='best')
            plt.grid(True)
            plt.show()

        return nodal_disp
