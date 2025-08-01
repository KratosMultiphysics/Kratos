import os
import shutil
import tempfile
import numpy as np

from material_editor import MaterialEditor
from project_parameter_editor import ProjectParameterEditor
from mdpa_editor import MdpaEditor
from generic_test_runner import GenericTestRunner
from plots import (
    plot_delta_sigma_triaxial, plot_volumetric_vertical_strain_triaxial, plot_principal_stresses_triaxial,
    plot_p_q_triaxial, plot_mohr_coulomb_triaxial,
    plot_strain_stress_direct_shear, plot_principal_stresses_direct_shear,
    plot_p_q_direct_shear, plot_mohr_coulomb_direct_shear
)

def setup_simulation_files(test_type, tmp_folder):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.join(base_dir, f"test_{test_type}")
    for filename in ["MaterialParameters.json", "ProjectParameters.json", "mesh.mdpa"]:
        shutil.copy(os.path.join(src_dir, filename), tmp_folder)
    return (
        os.path.join(tmp_folder, "MaterialParameters.json"),
        os.path.join(tmp_folder, "ProjectParameters.json"),
        os.path.join(tmp_folder, "mesh.mdpa")
    )

def set_material_constitutive_law(json_file_path, dll_path, material_parameters, index):
    editor = MaterialEditor(json_file_path)
    if dll_path:
        editor.update_material_properties({
            "IS_FORTRAN_UDSM": True,
            "UMAT_PARAMETERS": material_parameters,
            "UDSM_NAME": dll_path,
            "UDSM_NUMBER": index
        })
        editor.set_constitutive_law("SmallStrainUDSM2DPlaneStrainLaw")
    else:
        editor.update_material_properties({
            "YOUNG_MODULUS": material_parameters[0],
            "POISSON_RATIO": material_parameters[1]
        })
        editor.set_constitutive_law("GeoLinearElasticPlaneStrain2DLaw")

def set_project_parameters(project_path, num_steps, end_time, initial_stress):
    editor = ProjectParameterEditor(project_path)
    editor.update_property('time_step', end_time / num_steps)
    editor.update_property('end_time', end_time)
    stress_vector = [-initial_stress] * 3 + [0.0]
    editor.update_nested_value("apply_initial_uniform_stress_field", "value", stress_vector)

def set_mdpa(mdpa_path, max_strain, init_pressure, num_steps, end_time, test_type):
    editor = MdpaEditor(mdpa_path)
    editor.update_maximum_strain(max_strain)
    editor.update_end_time(end_time)
    editor.update_first_timestep(num_steps, end_time)
    if test_type == "triaxial":
        editor.update_initial_effective_cell_pressure(init_pressure)
    if test_type == "direct_shear":
        editor.update_middle_maximum_strain(max_strain)

def calculate_principal_stresses(tensors):
    sigma_1, sigma_3 = [], []
    for tensors_at_time in tensors.values():
        for sigma in tensors_at_time:
            eigenvalues, _ = np.linalg.eigh(sigma)
            sigma_1.append(np.min(eigenvalues))
            sigma_3.append(np.max(eigenvalues))
    return sigma_1, sigma_3

def get_cohesion_phi(umat_parameters, indices):
    if indices:
        c_idx, phi_idx = indices
        return float(umat_parameters[c_idx - 1]), float(umat_parameters[phi_idx - 1])
    return None, None

def plot_results(test_type, axes, yy, vol, sigma1, sigma3, shear_xy, shear_strain_xy, mean_stress, von_mises, cohesion, phi):
    if test_type == "triaxial":
        plot_delta_sigma_triaxial(axes[0], yy, abs(np.array(sigma1) - np.array(sigma3)))
        plot_volumetric_vertical_strain_triaxial(axes[1], yy, vol)
        plot_principal_stresses_triaxial(axes[2], sigma1, sigma3)
        plot_p_q_triaxial(axes[3], mean_stress, von_mises)
        plot_mohr_coulomb_triaxial(axes[4], sigma1[-1], sigma3[-1], cohesion, phi)
    elif test_type == "direct_shear":
        plot_strain_stress_direct_shear(axes[0], shear_strain_xy, shear_xy)
        plot_principal_stresses_direct_shear(axes[1], sigma1, sigma3)
        plot_p_q_direct_shear(axes[2], mean_stress, von_mises)
        plot_mohr_coulomb_direct_shear(axes[3], sigma1[-1], sigma3[-1], cohesion, phi)
    else:
        raise ValueError(f"Unsupported test_type: {test_type}")

def run_simulation(test_type, dll_path, index, material_parameters, num_steps, end_time, maximum_strain,
                   initial_effective_cell_pressure, cohesion_phi_indices=None, axes=None):
    tmp_folder = tempfile.mkdtemp(prefix=f"{test_type}_")

    try:
        json_path, project_path, mdpa_path = setup_simulation_files(test_type, tmp_folder)

        set_material_constitutive_law(json_path, dll_path, material_parameters, index)
        set_project_parameters(project_path, num_steps, end_time, initial_effective_cell_pressure)
        set_mdpa(mdpa_path, maximum_strain, initial_effective_cell_pressure, num_steps, end_time, test_type)

        runner = GenericTestRunner([os.path.join(tmp_folder, 'gid_output', "output.post.res")], tmp_folder)
        tensors, yy_strain, vol_strain, von_mises, mean_stress, shear_xy, shear_strain_xy = runner.run()

        sigma_1, sigma_3 = calculate_principal_stresses(tensors)
        cohesion, friction_angle = get_cohesion_phi(material_parameters, cohesion_phi_indices)

        if axes:
            for ax in axes:
                ax.clear()

        plot_results(test_type, axes, yy_strain, vol_strain, sigma_1, sigma_3,
                     shear_xy, shear_strain_xy, mean_stress, von_mises, cohesion, friction_angle)

    finally:
        if os.path.exists(tmp_folder):
            shutil.rmtree(tmp_folder)
