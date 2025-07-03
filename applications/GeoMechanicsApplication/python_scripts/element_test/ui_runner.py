import traceback
from ui_logger import log_message
from run_simulation import run_simulation

def run_gui_builder(test_type, dll_path, index, material_parameters, input_widgets, cohesion_phi_indices, axes):
    try:

        sigma_init = float(input_widgets["Initial effective cell pressure |σ'ₓₓ|"].get())
        eps_max = float(input_widgets["Maximum Strain |εᵧᵧ|"].get())
        n_steps = float(input_widgets["Number of steps"].get())
        duration = float(input_widgets["Duration"].get())

        if any(val <= 0 for val in [eps_max, n_steps, duration]) or sigma_init < 0:
            raise ValueError("All values must be positive and non-zero.")

        run_simulation(
            test_type=test_type,
            dll_path=dll_path or "",
            index=index,
            material_parameters=material_parameters,
            num_steps=n_steps,
            end_time=duration,
            maximum_strain=eps_max,
            initial_effective_cell_pressure=sigma_init,
            cohesion_phi_indices=cohesion_phi_indices,
            axes=axes
        )

    except Exception:
        log_message("Error during simulation:", "error")
        log_message(traceback.format_exc(), "error")
