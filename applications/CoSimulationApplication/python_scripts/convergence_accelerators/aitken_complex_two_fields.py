import KratosMultiphysics as KM

from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_convergence_accelerator import CoSimulationConvergenceAccelerator
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
from copy import deepcopy
from collections import deque


def Create(settings):
    cs_tools.SettingsTypeCheck(settings)
    return AitkenComplexTwoFieldsConvergenceAccelerator(settings)


class AitkenComplexTwoFieldsConvergenceAccelerator(CoSimulationConvergenceAccelerator):

    def __init__(self, settings):
        super().__init__(settings)

        self.R = deque(maxlen=2)

        self.alpha_old = complex(
            self.settings["init_alpha_real"].GetDouble(),
            self.settings["init_alpha_imag"].GetDouble()
        )

        self.init_alpha_abs_max = self.settings["init_alpha_abs_max"].GetDouble()
        self.alpha_abs_max = self.settings["alpha_abs_max"].GetDouble()
        self.alpha_abs_min = self.settings["alpha_abs_min"].GetDouble()
        self.max_alpha_phase_deg = self.settings["max_alpha_phase_deg"].GetDouble()

        self.track_residual = self.settings["track_residual"].GetBool()
        self.residual_output_file = self.settings["residual_output_file"].GetString()

        self.residual_history = []
        self.coupling_iteration = 0

    def InitializeSolutionStep(self):
        self.initial_iteration = True
        self.R.clear()
        self.residual_history = []
        self.coupling_iteration = 0

    def _LimitComplexAlphaPhase(self, alpha):
        alpha_abs = abs(alpha)

        if alpha_abs < 1e-30:
            return 0.0 + 0.0j

        phase = np.angle(alpha)
        max_phase = np.deg2rad(self.max_alpha_phase_deg)

        reference_phase = 0.0 if np.real(alpha) >= 0.0 else np.pi

        phase_error = np.angle(
            np.exp(1j * (phase - reference_phase))
        )

        phase_error = np.clip(
            phase_error,
            -max_phase,
            max_phase
        )

        limited_phase = reference_phase + phase_error

        return alpha_abs * np.exp(1j * limited_phase)

    def _TrackResidual(self, residual_norm, alpha):
        if not self.track_residual:
            return

        if self.coupling_iteration == 0:
            relative_residual_norm = 1.0
        else:
            previous_residual_norm = self.residual_history[-1][1]
            relative_residual_norm = (
                residual_norm / previous_residual_norm
                if previous_residual_norm > 0.0
                else residual_norm
            )

        self.residual_history.append([
            self.coupling_iteration + 1,
            residual_norm,
            relative_residual_norm,
            np.real(alpha),
            np.imag(alpha),
            abs(alpha),
            np.rad2deg(np.angle(alpha))
        ])

    def UpdateSolution(self, r_pair, x_pair):
        r_real = np.asarray(r_pair[0], dtype=float)
        r_imag = np.asarray(r_pair[1], dtype=float)

        if r_real.shape != r_imag.shape:
            raise Exception(
                "AitkenComplexTwoFieldsConvergenceAccelerator: "
                f"residual fields have different sizes: {r_real.shape} and {r_imag.shape}"
            )

        r_complex = r_real + 1j * r_imag
        self.R.appendleft(deepcopy(r_complex))

        if self.initial_iteration:
            self.initial_iteration = False

            alpha = self.alpha_old

            if abs(alpha) > self.init_alpha_abs_max:
                alpha = self.init_alpha_abs_max * alpha / abs(alpha)

            alpha = self._LimitComplexAlphaPhase(alpha)

        else:
            r_diff = self.R[0] - self.R[1]

            numerator = np.vdot(self.R[1], r_diff)
            denominator = np.vdot(r_diff, r_diff)

            if abs(denominator) < 1e-30:
                alpha = self.alpha_old

                if self.echo_level > 0:
                    cs_tools.cs_print_warning(
                        self._ClassName(),
                        f"Aitken denominator close to zero. Reusing alpha = {alpha}"
                    )
            else:
                alpha = -self.alpha_old * numerator / denominator

            alpha_abs = abs(alpha)

            if alpha_abs > self.alpha_abs_max:
                alpha = self.alpha_abs_max * alpha / alpha_abs

                if self.echo_level > 0:
                    cs_tools.cs_print_warning(
                        self._ClassName(),
                        f"|alpha| reached upper bound: {self.alpha_abs_max}"
                    )

            elif alpha_abs < self.alpha_abs_min and alpha_abs > 1e-30:
                alpha = self.alpha_abs_min * alpha / alpha_abs

                if self.echo_level > 0:
                    cs_tools.cs_print_warning(
                        self._ClassName(),
                        f"|alpha| reached lower bound: {self.alpha_abs_min}"
                    )

            alpha = self._LimitComplexAlphaPhase(alpha)

        self.alpha_old = alpha

        residual_norm = np.linalg.norm(r_complex)
        self._TrackResidual(residual_norm, alpha)
        self.coupling_iteration += 1

        if self.echo_level > 2:
            cs_tools.cs_print_info(
                self._ClassName(),
                f"Complex two-field Aitken alpha = {alpha}, "
                f"|alpha| = {abs(alpha)}, "
                f"phase = {np.rad2deg(np.angle(alpha))} deg"
            )

        delta_complex = alpha * r_complex

        return np.real(delta_complex), np.imag(delta_complex)

    def FinalizeSolutionStep(self):
        if self.track_residual and self.residual_history:
            np.savetxt(
                self.residual_output_file,
                np.asarray(self.residual_history),
                delimiter=",",
                header=(
                    "iteration,residual_norm,"
                    "relative_residual_norm,"
                    "alpha_real,alpha_imag,"
                    "alpha_abs,alpha_phase_deg"
                ),
                comments=""
            )

        self.R.clear()
        self.residual_history = []
        self.coupling_iteration = 0

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "init_alpha_real"      :  0.05,
            "init_alpha_imag"      :  0.0,
            "init_alpha_abs_max"   :  0.20,
            "alpha_abs_max"        :  2.0,
            "alpha_abs_min"        :  0.0,
            "max_alpha_phase_deg"  :  20.0,
            "track_residual"       :  false,
            "residual_output_file" : "residual_history_complex_aitken.csv"
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults