# class for torsional bar

# import python modules
from scipy. optimize import minimize
from scipy import linalg
import matplotlib.pyplot as plt
import numpy as np


class TorsionalBar():

    def __init__(self, properties, rdof, target_freq):

        self.properties = properties
        self.rdof = rdof
        self.GJ = self.optimize(target_freq)
        self.K, self.M, self.B, self.K_big, self.M_big, self.B_big = self.torsional_bar(self.GJ)

        self.eig_vals, self.eig_vecs_raw, self.eig_freq, self.eig_per = self.eigen_value(
            self.K, self.M)

        self.B, self.B_big = self.damping(2, 0.01)

    def torsional_bar(self, GJ):

        elem_number = self.properties.levels
        elem_length = self.properties.height / elem_number

        rho = self.properties.density
        material_area = self.properties.length * self.properties.width
        I_0 = self.properties.width ** 2 + self.properties.length ** 2

        # mass matrix
        def mass_matrix():
            # material constant
            material_constant = rho * I_0 * elem_length / 6

            # element matrix
            m = np.array([[2, 1],
                          [1, 2]])

            # global matrix
            M = np.zeros((elem_number + 1, elem_number + 1))

            for i in range(elem_number):
                    M_temp = np.zeros((elem_number + 1, elem_number + 1))
                    M_temp[i:i + 2, i:i + 2] = m
                    M += material_constant * M_temp

            M_big = M

            # remove the fixed degrees of freedom
            for dof in self.rdof:
                for i in [0, 1]:
                    M = np.delete(M, dof, axis=i)

            # return stiffness matrix
            return M, M_big

        # stiffness matrix
        def stiffness_matrix():
            # material constant
            stiffness_constant = GJ / elem_length

            # element matrix
            k = np.array([[1, -1],
                          [-1, 1]])

            # global matrix
            K = np.zeros((elem_number + 1, elem_number + 1))

            for i in range(elem_number):
                    K_temp = np.zeros((elem_number + 1, elem_number + 1))
                    K_temp[i:i + 2, i:i + 2] = k
                    K += stiffness_constant * K_temp

            K_big = K

            # remove the fixed degrees of freedom
            for dof in self.rdof:
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
            GJ = param
            K = self.torsional_bar(GJ)[0]
            M = self.torsional_bar(GJ)[1]
            return (self.eigen_value(K, M)[2][0] - target_freq) ** 2

        init = 1
        result = minimize(objective, init, method='Nelder-Mead')
        return result.x

    def damping(self, nr_modes, daming_ratio):

        # Damping ratios
        xi = [daming_ratio] * nr_modes

        # LHS assambly
        LHS = np.zeros((nr_modes, nr_modes))  # nModes x nModes matrix

        for i in range(nr_modes):
            for j in range(nr_modes):
                LHS[i][j] = pow(self.eig_vals[i], 2 * j - 1)

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

        return C, C_big


    def load_distribution(self, fluid_forces):

        if self.direction == 'X':
            load = fluid_forces[0]
        elif self.direction == 'Y':
            load = fluid_forces[1]
        elif self.direction == 'R':
            load = fluid_forces[2]

        elem_number = self.properties.levels

        def nodal_force(force):
            F = [force / 2, force / 2]
            return F

        nodal_load = list(map(nodal_force, load))

        load = np.zeros(elem_number + 1)

        for i in range(elem_number):
            load_temp = np.zeros(elem_number + 1)
            load_temp[i:i + 2] = nodal_load[i]
            load += load_temp

        # remove the fixed degrees of freedom
        for dof in self.rdof:
            load = np.delete(load, dof, axis=0)

        return load

    def eigen_value_load(self, mode):

        elem_number = self.properties.levels
        elem_length = self.properties.height / elem_number

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
            extended_eig_vecs[1:, :] = eig_vecs_norm

            # Plot eigenvalues and shapes
            fig = plt.figure(1)

            # Number of modes to iterate
            modes = 1
            for i in range(modes):
                plt.plot(extended_eig_vecs[:, freq_ind[i]], Z, marker='o', label="mode " + str(i + 1) + ", eigFreq: " + str(np.round(self.eig_freq[freq_ind[i]], 3))
                         + " [Hz], T: " + str(np.round(self.eig_per[freq_ind[i]], 3)) + " [s], omega: " + str(np.round(self.eig_vals[freq_ind[i]], 3)) + " [rad/s]")

            plt.title("Eigenmodes - mass normalized eigenvectors")
            plt.xlabel("Eigenmode magnitude")
            plt.ylabel("Height [m]")
            plt.legend(loc='best')
            plt.grid(True)
            plt.show()

        return nodal_disp
