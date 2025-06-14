import os
import shutil
import tempfile
import numpy as np

from material_editor import MaterialEditor
from project_parameter_editor import ProjectParameterEditor
from mdpa_editor import MdpaEditor
from triaxial import TriaxialTestRunner
from triaxial_plots import plot_delta_sigma, plot_volumetric_strain, plot_sigma, plot_p_q, plot_mohr_coulomb_circle


def run_triaxial_simulation(dll_path, index, umat_parameters, num_steps, end_time, maximum_strain, initial_effective_cell_pressure,
                            cohesion_phi_indices=None, axes=None):
    tmp_folder = tempfile.mkdtemp(prefix="triaxial_")

    base_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.join(base_dir, "test_triaxial")
    files_to_copy = [
        "MaterialParameters.json",
        "ProjectParameters.json",
        "triaxial.mdpa"
    ]

    for filename in files_to_copy:
        shutil.copy(os.path.join(src_dir, filename), tmp_folder)

    json_file_path = os.path.join(tmp_folder, "MaterialParameters.json")
    project_param_path = os.path.join(tmp_folder, "ProjectParameters.json")
    mdpa_path = os.path.join(tmp_folder, "triaxial.mdpa")

    material_editor = MaterialEditor(json_file_path)

    if dll_path:
        material_editor.update_material_properties({
            "IS_FORTRAN_UDSM": True,
            "UMAT_PARAMETERS": umat_parameters,
            "UDSM_NAME": dll_path,
            "UDSM_NUMBER": index
        })
        material_editor.set_constitutive_law("SmallStrainUDSM2DPlaneStrainLaw")
    else:
        material_editor.update_material_properties({
            "YOUNG_MODULUS": umat_parameters[0],
            "POISSON_RATIO": umat_parameters[1]
        })
        material_editor.set_constitutive_law("GeoLinearElasticPlaneStrain2DLaw")

    project_editor = ProjectParameterEditor(project_param_path)
    project_editor.update_time_step_properties(num_steps, end_time)
    project_editor.update_end_time_properties(end_time)
    project_editor.update_initial_stress_vector(initial_effective_cell_pressure)

    mdpa_editor = MdpaEditor(mdpa_path)
    mdpa_editor.update_maximum_strain(maximum_strain)
    mdpa_editor.update_initial_effective_cell_pressure(initial_effective_cell_pressure)
    mdpa_editor.update_end_time(end_time)
    mdpa_editor.update_first_timestep(num_steps)

    output_files = [os.path.join(tmp_folder, 'gid_output', "triaxial.post.res")]
    runner = TriaxialTestRunner(output_files, tmp_folder)
    reshaped_values_by_time, vertical_strain, volumetric_strain, von_mises_stress, mean_effective_stresses = runner.run()

    if cohesion_phi_indices:
        c_idx, phi_idx = cohesion_phi_indices
        cohesion = float(umat_parameters[c_idx - 1])
        friction_angle = float(umat_parameters[phi_idx - 1])
    else:
        cohesion = None
        friction_angle = None

    sigma1_list = []
    sigma3_list = []
    eigenvector_list = []
    eigenvalue_list = []
    for time_step, tensors in reshaped_values_by_time.items():
        for sigma in tensors:
            eigenvalues, eigenvectors = np.linalg.eigh(sigma)
            eigenvector_list.append(eigenvectors)
            eigenvalue_list.append(eigenvalues)
            sigma_max = np.max(eigenvalues)
            sigma_min = np.min(eigenvalues)
            sigma1_list.append(sigma_min)
            sigma3_list.append(sigma_max)

    sigma_diff = abs(np.array(sigma1_list) - np.array(sigma3_list))

    if axes:
        for ax in axes:
            ax.clear()

    plot_delta_sigma(axes[0], vertical_strain, sigma_diff)
    plot_volumetric_strain(axes[1], vertical_strain, volumetric_strain)
    plot_sigma(axes[2], sigma1_list, sigma3_list)
    plot_p_q(axes[3], mean_effective_stresses, von_mises_stress)
    plot_mohr_coulomb_circle(axes[4], sigma1_list[-1], sigma3_list[-1], cohesion, friction_angle)
