from set_triaxial_test import lab_test


def run_triaxial_simulation(dll_path, model_index, parameters, num_steps, end_time, maximum_strain, initial_effective_cell_pressure):
    return lab_test(
        dll_path,
        model_index,
        parameters,
        num_steps,
        end_time,
        maximum_strain,
        initial_effective_cell_pressure
    )