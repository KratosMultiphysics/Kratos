import KratosMultiphysics as KM

from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_accelerator import CoSimulationConvergenceAccelerator

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import (
    cs_print_info,
    cs_print_warning,
    SettingsTypeCheck
)
import KratosMultiphysics.CoSimulationApplication.colors as colors

import numpy as np
from copy import deepcopy
from collections import deque


def Create(settings):
    SettingsTypeCheck(settings)
    return IQNILSComplexTwoFieldsConvergenceAccelerator(settings)


class IQNILSComplexTwoFieldsConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    """
    IQN-ILS for two real Kratos fields treated internally as one complex field.
    """

    def __init__(self, settings):
        super().__init__(settings)

        iteration_horizon = self.settings["iteration_horizon"].GetInt()
        timestep_horizon = self.settings["timestep_horizon"].GetInt()

        self.alpha = complex(
            self.settings["alpha_real"].GetDouble(),
            self.settings["alpha_imag"].GetDouble()
        )

        self.R = deque(maxlen=iteration_horizon)
        self.X = deque(maxlen=iteration_horizon)

        self.q = timestep_horizon - 1
        self.v_old_matrices = deque(maxlen=self.q)
        self.w_old_matrices = deque(maxlen=self.q)

        self.V_new = []
        self.W_new = []
        self.V_old = []
        self.W_old = []

        # --------------------------------------------------
        # Residual history output
        # --------------------------------------------------
        self.residual_history = []
        self.coupling_iteration = 0
        self.initial_residual_norm = None

        self.output_file = self.settings["residual_output_file"].GetString()

    def UpdateSolution(self, r_pair, x_pair):
        r_real = np.asarray(r_pair[0], dtype=float)
        r_imag = np.asarray(r_pair[1], dtype=float)

        x_real = np.asarray(x_pair[0], dtype=float)
        x_imag = np.asarray(x_pair[1], dtype=float)

        if r_real.shape != r_imag.shape:
            raise Exception(
                "IQNILSComplexTwoFieldsConvergenceAccelerator: "
                f"residual fields have different sizes: {r_real.shape} and {r_imag.shape}"
            )

        if x_real.shape != x_imag.shape:
            raise Exception(
                "IQNILSComplexTwoFieldsConvergenceAccelerator: "
                f"solution fields have different sizes: {x_real.shape} and {x_imag.shape}"
            )

        if r_real.shape != x_real.shape or r_imag.shape != x_imag.shape:
            raise Exception(
                "IQNILSComplexTwoFieldsConvergenceAccelerator: residual and solution sizes do not match."
            )

        r_complex = r_real + 1j * r_imag
        x_complex = x_real + 1j * x_imag

        # --------------------------------------------------
        # Store residual norm history
        # --------------------------------------------------
        residual_norm = np.linalg.norm(r_complex)

        if self.coupling_iteration == 0:
            relative_residual_norm = 1.0
        else:
            previous_residual_norm = self.residual_history[-1][1]

            if previous_residual_norm > 0.0:
                relative_residual_norm = (
                    residual_norm / previous_residual_norm
                )
            else:
                relative_residual_norm = residual_norm

        self.residual_history.append([
            self.coupling_iteration+1,
            residual_norm,
            relative_residual_norm
        ])

        self.coupling_iteration += 1

        delta_complex = self._ComputeIQNILSUpdate(r_complex, x_complex)

        return np.real(delta_complex), np.imag(delta_complex)

    def _ComputeIQNILSUpdate(self, r, x):
        self.R.appendleft(deepcopy(r))
        self.X.appendleft(x + r)  # r = x_tilde - x

        row = len(r)
        col = len(self.R) - 1
        k = col
        num_old_matrices = len(self.v_old_matrices)

        if len(self.V_old) == 0 and len(self.W_old) == 0:
            if k == 0:
                if self.echo_level > 3:
                    cs_print_info(
                        self._ClassName(),
                        "Doing complex relaxation in first iteration with alpha = ",
                        str(self.alpha)
                    )
                return self.alpha * r

            else:
                if self.echo_level > 3:
                    cs_print_info(self._ClassName(), "Doing complex IQN-ILS extrapolation")
                    cs_print_info(self._ClassName(), "Number of new modes: ", col)

                self.V_new = self._BuildDifferenceMatrix(self.R, col, row)
                V = self.V_new

                if (V.shape[0] < V.shape[1]) and self.echo_level > 0:
                    cs_print_warning(
                        self._ClassName(),
                        ": " + colors.red("WARNING: column number larger than row number!")
                    )

                self.W_new = self._BuildDifferenceMatrix(self.X, col, row)
                W = self.W_new

                delta_r = -self.R[0]

                c = np.linalg.lstsq(V, delta_r, rcond=None)[0]

                delta_x = np.dot(W, c) - delta_r

                return delta_x

        else:
            if k == 0:
                if self.echo_level > 3:
                    cs_print_info(self._ClassName(), "Using complex matrices from previous time steps")
                    cs_print_info(self._ClassName(), "Number of previous matrices: ", num_old_matrices)

                V = self.V_old
                W = self.W_old

                delta_r = -self.R[0]

                c = np.linalg.lstsq(V, delta_r, rcond=None)[0]

                delta_x = np.dot(W, c) - delta_r

                return delta_x

            else:
                if self.echo_level > 3:
                    cs_print_info(self._ClassName(), "Doing complex IQN-ILS extrapolation")
                    cs_print_info(self._ClassName(), "Number of new modes: ", col)
                    cs_print_info(self._ClassName(), "Number of previous matrices: ", num_old_matrices)

                self.V_new = self._BuildDifferenceMatrix(self.R, col, row)
                V = np.hstack((self.V_new, self.V_old))

                if (V.shape[0] < V.shape[1]) and self.echo_level > 0:
                    cs_print_warning(
                        self._ClassName(),
                        ": " + colors.red("WARNING: column number larger than row number!")
                    )

                self.W_new = self._BuildDifferenceMatrix(self.X, col, row)
                W = np.hstack((self.W_new, self.W_old))

                delta_r = -self.R[0]

                c = np.linalg.lstsq(V, delta_r, rcond=None)[0]

                delta_x = np.dot(W, c) - delta_r

                return delta_x

    def _BuildDifferenceMatrix(self, deque_data, col, row):
        mat = np.empty(shape=(col, row), dtype=complex)

        for i in range(col):
            mat[i] = deque_data[i] - deque_data[i + 1]

        return mat.T

    def FinalizeSolutionStep(self):
        # --------------------------------------------------
        # Write residual history for this frequency/step
        # --------------------------------------------------
        if self.residual_history:
            np.savetxt(
                self.output_file,
                np.asarray(self.residual_history),
                delimiter=",",
                header="iteration,residual_norm,relative_residual_norm",
                comments=""
            )

        if len(self.V_new) != 0 and len(self.W_new) != 0:
            self.v_old_matrices.appendleft(self.V_new)
            self.w_old_matrices.appendleft(self.W_new)

        if self.v_old_matrices and self.w_old_matrices:
            self.V_old = np.concatenate(self.v_old_matrices, 1)
            self.W_old = np.concatenate(self.w_old_matrices, 1)

        if self.R and self.X:
            if self.echo_level > 3:
                cs_print_info(self._ClassName(), "Cleaning")
            self.R.clear()
            self.X.clear()

        self.V_new = []
        self.W_new = []

        self.residual_history = []
        self.coupling_iteration = 0
        self.initial_residual_norm = None

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "iteration_horizon"     : 20,
            "timestep_horizon"      : 1,
            "alpha_real"            : 0.01,
            "alpha_imag"            : 0.0,
            "residual_output_file"  : "residual_history_complex_iqnils.csv"
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults