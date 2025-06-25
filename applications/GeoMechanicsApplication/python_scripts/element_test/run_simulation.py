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

def run_simulation(test_type, dll_path, index, umat_parameters, num_steps, end_time, maximum_strain,
                   initial_effective_cell_pressure, cohesion_phi_indices=None, axes=None):
    tmp_folder = tempfile.mkdtemp(prefix=f"{test_type}_")

    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        src_dir = os.path.join(base_dir, f"test_{test_type}")
        files_to_copy = [
            "MaterialParameters.json",
            "ProjectParameters.json",
            "mesh.mdpa"
        ]

        for filename in files_to_copy:
            shutil.copy(os.path.join(src_dir, filename), tmp_folder)

        json_file_path = os.path.join(tmp_folder, "MaterialParameters.json")
        project_param_path = os.path.join(tmp_folder, "ProjectParameters.json")
        mdpa_path = os.path.join(tmp_folder, "mesh.mdpa")

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
        project_editor.update_property('time_step', end_time / num_steps)
        project_editor.update_property('end_time', end_time)
        initial_uniform_stress_value = [-initial_effective_cell_pressure] * 3 + [0.0]
        project_editor.update_nested_value("apply_initial_uniform_stress_field", "value", initial_uniform_stress_value)

        mdpa_editor = MdpaEditor(mdpa_path)
        mdpa_editor.update_maximum_strain(maximum_strain)
        mdpa_editor.update_initial_effective_cell_pressure(initial_effective_cell_pressure)
        mdpa_editor.update_end_time(end_time)
        mdpa_editor.update_first_timestep(num_steps, end_time)

        if test_type == "direct_shear":
            mdpa_editor.update_middle_maximum_strain(maximum_strain)

        output_files = [os.path.join(tmp_folder, 'gid_output', "output.post.res")]
        runner = GenericTestRunner(output_files, tmp_folder)
        (reshaped_values_by_time, vertical_strain, volumetric_strain, von_mises_stress, mean_effective_stresses,
         shear_stress_xy, shear_strain_xy) = runner.run()

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

        if test_type == "triaxial":
            plot_delta_sigma_triaxial(axes[0], vertical_strain, sigma_diff)
            plot_volumetric_vertical_strain_triaxial(axes[1], vertical_strain, volumetric_strain)
            plot_principal_stresses_triaxial(axes[2], sigma1_list, sigma3_list)
            plot_p_q_triaxial(axes[3], mean_effective_stresses, von_mises_stress)
            plot_mohr_coulomb_triaxial(axes[4], sigma1_list[-1], sigma3_list[-1], cohesion, friction_angle)


        elif test_type == "direct_shear":
            plot_strain_stress_direct_shear(axes[0], shear_strain_xy, shear_stress_xy)
            plot_principal_stresses_direct_shear(axes[1], sigma1_list, sigma3_list)
            plot_p_q_direct_shear(axes[2], mean_effective_stresses, von_mises_stress)
            plot_mohr_coulomb_direct_shear(axes[3], sigma1_list[-1], sigma3_list[-1], cohesion, friction_angle)

        else:
            raise ValueError(f"Unsupported test_type: {test_type}")

        return

    finally:
        if os.path.exists(tmp_folder):
            shutil.rmtree(tmp_folder)