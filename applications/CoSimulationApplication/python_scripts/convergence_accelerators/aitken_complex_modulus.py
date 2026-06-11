import KratosMultiphysics as KM

from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_accelerator import CoSimulationConvergenceAccelerator
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
from copy import deepcopy
from collections import deque


def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return AitkenComplexModulusResidualConvergenceAccelerator(settings)


class AitkenComplexModulusResidualConvergenceAccelerator(CoSimulationConvergenceAccelerator):
    """
    Real Aitken relaxation computed from the component-wise modulus
    of a complex residual.

    Complex residual:
        r_c = r_real + i r_imag

    Scalar residual used for Aitken:
        r_abs = |r_c|

    Real relaxation:
        alpha in R

    Complex update:
        delta_c = alpha * r_c
    """

    def __init__(self, settings):
        super().__init__(settings)

        self.R_abs = deque(maxlen=2)

        self.alpha_old = self.settings["init_alpha"].GetDouble()

        self.init_alpha_abs_max = self.settings["init_alpha_abs_max"].GetDouble()
        self.alpha_abs_max = self.settings["alpha_abs_max"].GetDouble()
        self.alpha_abs_min = self.settings["alpha_abs_min"].GetDouble()

        # Residual history output
        self.write_residual_history = self.settings[
            "write_residual_history"
        ].GetBool()

        self.output_file = self.settings[
            "residual_output_file"
        ].GetString()

        self.residual_history = []
        self.coupling_iteration = 0

        # Previous residual norm
        self.previous_residual_norm = None

    def InitializeSolutionStep(self):

        self.initial_iteration = True

        self.R_abs.clear()

        self.residual_history = []

        self.coupling_iteration = 0

        self.previous_residual_norm = None

    def UpdateSolution(self, r_pair, x_pair):

        r_real = np.asarray(r_pair[0], dtype=float)

        r_imag = np.asarray(r_pair[1], dtype=float)

        if r_real.shape != r_imag.shape:
            raise Exception(
                "AitkenComplexModulusResidualConvergenceAccelerator: "
                f"residual fields have different sizes: "
                f"{r_real.shape} and {r_imag.shape}"
            )

        # Complex residual
        r_complex = r_real + 1j * r_imag

        # Residual used for Aitken
        r_abs = np.abs(r_complex)

        self.R_abs.appendleft(deepcopy(r_abs))

        # Residual norms
        residual_norm = np.linalg.norm(r_complex)

        if self.previous_residual_norm is None:

            relative_residual_norm = 1.0

        else:

            relative_residual_norm = (
                residual_norm / self.previous_residual_norm
                if self.previous_residual_norm > 0.0
                else residual_norm
            )

        self.previous_residual_norm = residual_norm

        # First iteration
        if self.initial_iteration:

            self.initial_iteration = False

            alpha = self.alpha_old

            if abs(alpha) > self.init_alpha_abs_max:

                alpha = np.sign(alpha) * self.init_alpha_abs_max

        # Standard Aitken update
        else:

            r_diff = self.R_abs[0] - self.R_abs[1]

            numerator = np.dot(self.R_abs[1], r_diff)

            denominator = np.dot(r_diff, r_diff)

            if abs(denominator) < 1e-30:

                alpha = self.alpha_old

                if self.echo_level > 0:

                    cs_tools.cs_print_warning(
                        self._ClassName(),
                        f"Aitken denominator close to zero. "
                        f"Reusing alpha = {alpha}"
                    )

            else:

                alpha = -self.alpha_old * numerator / denominator

            # Upper bound
            if abs(alpha) > self.alpha_abs_max:

                alpha = np.sign(alpha) * self.alpha_abs_max

                if self.echo_level > 0:

                    cs_tools.cs_print_warning(
                        self._ClassName(),
                        f"|alpha| reached upper bound: "
                        f"{self.alpha_abs_max}"
                    )

            # Lower bound
            elif abs(alpha) < self.alpha_abs_min and abs(alpha) > 1e-30:

                alpha = np.sign(alpha) * self.alpha_abs_min

                if self.echo_level > 0:

                    cs_tools.cs_print_warning(
                        self._ClassName(),
                        f"|alpha| reached lower bound: "
                        f"{self.alpha_abs_min}"
                    )

        self.alpha_old = alpha

        if self.echo_level > 2:

            cs_tools.cs_print_info(
                self._ClassName(),
                f"Modulus-based real Aitken alpha = {alpha}"
            )

        # Store residual history
        if self.write_residual_history:

            self.residual_history.append([
                self.coupling_iteration + 1,
                residual_norm,
                relative_residual_norm,
                alpha
            ])

        self.coupling_iteration += 1
        
        # Complex update
        delta_complex = alpha * r_complex

        return np.real(delta_complex), np.imag(delta_complex)

    def FinalizeSolutionStep(self):

        if (
            self.write_residual_history
            and len(self.residual_history) > 0
        ):

            np.savetxt(
                self.output_file,
                np.asarray(self.residual_history),
                delimiter=",",
                header=(
                    "iteration,"
                    "residual_norm,"
                    "relative_residual_norm,"
                    "alpha"
                ),
                comments=""
            )

        self.R_abs.clear()

        self.residual_history = []

        self.coupling_iteration = 0

        self.previous_residual_norm = None

    @classmethod
    def _GetDefaultParameters(cls):

        this_defaults = KM.Parameters("""{
            "init_alpha"              : 0.05,
            "init_alpha_abs_max"      : 0.20,
            "alpha_abs_max"           : 2.0,
            "alpha_abs_min"           : 0.0,
            "write_residual_history"  : false,
            "residual_output_file"    : "residual_history_aitken_complex_modulus.csv"
        }""")

        this_defaults.AddMissingParameters(
            super()._GetDefaultParameters()
        )

        return this_defaults