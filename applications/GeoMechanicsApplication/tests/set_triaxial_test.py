import os
import json
import numpy as np
import plotly.graph_objects as go
import re
import tkinter as tk
from tkinter import messagebox, filedialog

import KratosMultiphysics as Kratos

from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis
import test_helper


class TriaxialTest:
    def __init__(self, json_file_path):
        self.json_file_path = json_file_path
        self.data = self._read_json()

    def _read_json(self):
        with open(self.json_file_path, 'r') as file:
            return json.load(file)

    def _write_json(self):
        with open(self.json_file_path, 'w') as file:
            json.dump(self.data, file, indent=4)

    def modify_umat_parameters(self, new_value):
        self.data['properties'][0]['Material']['Variables']['UMAT_PARAMETERS'][0] = new_value
        self._write_json()

    def read_umat_parameters(self):
        umat_parameters = self.data['properties'][0]['Material']['Variables']['UMAT_PARAMETERS']
        cohesion = umat_parameters[2]
        friction_angle = umat_parameters[3]
        return cohesion, friction_angle

def run_triaxial_test(output_file_paths):
    os.getcwd()

    # construct parameterfile names of stages to run
    project_path = os.path.join('test_triaxial')
    n_stages = 2
    parameter_file_names = [os.path.join('ProjectParameters_stage' + str(i + 1) + '.json') for i in
                            range(n_stages)]

    # change to project directory
    os.chdir(project_path)

    # setup stages from parameterfiles
    parameters_stages = [None] * n_stages
    for idx, parameter_file_name in enumerate(parameter_file_names):
        with open(parameter_file_name, 'r') as parameter_file:
            parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    stages = [GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

    # execute the stages
    [stage.Run() for stage in stages]

    # # back to working directory
    # os.chdir(currentWorking)

    cauchy_stress_results = []
    displacement_results = []
    strain_results = []

    for output_file_path in output_file_paths:
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_path)

        # Iterate through all results in the output_data
        for result_name, result_items in output_data["results"].items():
            if result_name == "CAUCHY_STRESS_TENSOR":
                for result_item in result_items:
                    time = result_item["time"]
                    values = result_item["values"]
                    cauchy_stress_results.append({
                        "time": time,
                        "values": values
                    })
            elif result_name == "DISPLACEMENT":
                for result_item in result_items:
                    time = result_item["time"]
                    values = result_item["values"]
                    displacement_results.append({
                        "time": time,
                        "values": values
                    })
            elif result_name == "ENGINEERING_STRAIN_TENSOR":
                for result_item in result_items:
                    time = result_item["time"]
                    values = result_item["values"]
                    strain_results.append({
                        "time": time,
                        "values": values
                    })

    reshaped_values_by_time = {}

    # Process Cauchy stress results
    for idx, result in enumerate(cauchy_stress_results):
        time_step = result["time"]
        element_values = result["values"]

        if not element_values:
            continue

        reshaped_values = []
        element = element_values[0]
        sublist = element["value"][0]
        tensor = np.array([
            [sublist[0], sublist[3], sublist[5]],
            [sublist[3], sublist[1], sublist[4]],
            [sublist[5], sublist[4], sublist[2]],
        ])
        reshaped_values.append(tensor)

        reshaped_values_by_time[time_step] = reshaped_values


    volumetric_strains = []
    yy_strains = []

    for idx, result in enumerate(strain_results):
        element_values = result["values"]

        if not element_values:  # Skip empty lists
            continue

        element = element_values[0]
        sublist = element["value"][0]
        epsilon_xx = sublist[0]
        epsilon_yy = sublist[1]
        epsilon_zz = sublist[2]

        epsilon_v = epsilon_xx + epsilon_yy + epsilon_zz

        volumetric_strains.append(epsilon_v)
        yy_strains.append(epsilon_yy)

    return reshaped_values_by_time, yy_strains, volumetric_strains

class MaterialEditor:
    def __init__(self, json_path):
        self.json_path = json_path
        self.data = self._load_json()

    def _load_json(self):
        with open(self.json_path, 'r') as f:
            data = json.load(f)
        return data

    def _update_material_and_save(self, entries: dict):
        variables = self.data["properties"][0]["Material"]["Variables"]
        for key, entry in entries.items():
            # is entry a list
            value = entry
            if isinstance(entry, list):
                value_str = [str(x).strip() for x in entry]
                value = [self._convert_type(x) for x in value_str]
            variables[key] = value

        with open(self.json_path, 'w') as f:
            json.dump(self.data, f, indent=4)
        # messagebox.showinfo("Success", "Material parameters updated successfully!")

    def _convert_type(self, value):
        try:
            if '.' in value or 'e' in value.lower():
                return float(value)
            else:
                return int(value)
        except ValueError:
            return value

class ProjectParameterEditor:
    def __init__(self, json_path):
        self.json_path = json_path
        with open(self.json_path, 'r') as f:
            self.raw_text = f.read()

    def _write_back(self):
        with open(self.json_path, 'w') as f:
            f.write(self.raw_text)

    def update_time_step_properties(self, number_of_step):
        new_time_step = 1.0 / number_of_step # (end_time - start_time) / number_of_steps
        pattern = r'("time_step"\s*:\s*)([0-9eE+.\-]+)'
        replacement = r'\g<1>' + str(new_time_step)
        self.raw_text, count = re.subn(pattern, replacement, self.raw_text)
        # if count == 0:
        #     messagebox.showwarning("Warning", "Could not find 'time_step' to update.")
        # else:
        #     self._write_back()
        #     messagebox.showinfo("Success", f"'time_step' updated to {new_time_step}")

class MdpaEditor:
    def __init__(self, mdpa_path):
        self.mdpa_path = mdpa_path
        with open(self.mdpa_path, 'r') as f:
            self.raw_text = f.read()

    def update_maximum_strain(self, maximum_strain):
        pattern = r'(\s*)\$maximum_strain(\s*)'
        prescribed_displacement = -maximum_strain / 100

        def replacer(match):
            leading_ws = match.group(1)
            trailing_ws = match.group(2)
            return f"{leading_ws}{prescribed_displacement}{trailing_ws}"

        new_text, count = re.subn(pattern, replacer, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", "Could not find '$maximum_strain' to update.")
        else:
            self.raw_text = new_text
        with open(self.mdpa_path, 'w') as f:
            f.write(self.raw_text)
        # messagebox.showinfo("Success", f"'$maximum_strain' replaced with {prescribed_displacement}")


    def update_initial_effective_cell_pressure(self, initial_effective_cell_pressure):
        pattern = r'(\s*)\$initial_effective_cell_pressure(\s*)'

        def replacer(match):
            leading_ws = match.group(1)
            trailing_ws = match.group(2)
            return f"{leading_ws}{initial_effective_cell_pressure}{trailing_ws}"

        new_text, count = re.subn(pattern, replacer, self.raw_text)
        if count == 0:
            messagebox.showwarning("Warning", "Could not find '$initial_effective_cell_pressure' to update.")
        else:
            with open(self.mdpa_path, 'w') as f:
                f.write(new_text)
        #     messagebox.showinfo("Success", f"'$initial_effective_cell_pressure' replaced with {initial_effective_cell_pressure}")

def plot_sigma(sigma_1, sigma_3):
    """
    Plots the principal stresses σ₁ and σ₃.

    Args:
        sigma_1 (list): List of σ₁ values.
        sigma_3 (list): List of σ₃ values.
    """
    fig = go.Figure()

    # Add scatter plot for σ₁ vs σ₃
    fig.add_trace(go.Scatter(
        x=sigma_3,
        y=sigma_1,
        mode='markers',
        marker=dict(size=10, color='blue'),
        name='σ₁ vs σ₃'
    ))

    fig.update_layout(
        title=dict(
            text='σ₁ vs σ₃ Plot',
            x=0.5,  # Center the title
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title='σ₃ (Principal Stress 3) [kN/m²]',
            showline=True,
            autorange='reversed',
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        yaxis=dict(
            title=' σ₁ (Principal Stress 1) [kN/m²]',
            showline=True,
            autorange='reversed',
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        template='plotly_white',
    )
    fig.update_layout(
        xaxis=dict(rangemode='tozero'),
        yaxis=dict(rangemode='tozero'),
    )
    # fig.show()
    return fig

def plot_delta_sigma(vertical_strain, sigma_diff):
    """
    Plots the difference between σ₁ and σ₃ against displacement.

    Args:
        displacement (list): List of displacements.
        diff (list): List of differences between σ₁ and σ₃.
    """
    fig = go.Figure()

    # Add scatter plot for σ₁ vs σ₃
    fig.add_trace(go.Scatter(
        x=vertical_strain,
        y=sigma_diff,
        mode='markers',
        marker=dict(size=10, color='blue'),
        name='σ₁ vs σ₃'
    ))

    fig.update_layout(
        title=dict(
            text='|σ₁ - σ₃| vs Vertical Strain Plot',
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title='Vertical Strain [-]',
            autorange='reversed',
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        yaxis=dict(
            title=' |σ₁ - σ₃| [kN/m²]',
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        template='plotly_white',
    )
    fig.update_layout(
        xaxis=dict(rangemode='tozero'),
        yaxis=dict(rangemode='tozero'),
    )
    return fig

def plot_mohr_coulomb_circle(sigma_1, sigma_3):
    """
    Plots the Mohr-Coulomb circle with σ' on the x-axis and mobilized shear stress on the y-axis.

    Args:
        sigma_1 (float): Maximum principal stress.
        sigma_3 (float): Minimum principal stress.
    """
    # Calculate the center and radius of the Mohr-Coulomb circle
    center = (sigma_1 + sigma_3) / 2
    radius = (sigma_1 - sigma_3) / 2

    # Generate points for the circle
    theta = np.linspace(0, np.pi, 200)
    sigma = center + radius * np.cos(theta)
    tau = radius * np.sin(theta)

    # Generate the dashed line
    phi_rad = np.radians(TriaxialTest(os.path.join('MaterialParameters_stage1.json'))
                         .read_umat_parameters()[1])
    cohesion = TriaxialTest(os.path.join('MaterialParameters_stage1.json')).read_umat_parameters()[0]

    x_line = np.linspace(0, sigma_1, 200)
    y_line = x_line * np.tan(phi_rad) - cohesion

    fig = go.Figure()

    # Add the Mohr-Coulomb circle
    fig.add_trace(go.Scatter(
        x=sigma,
        y=tau,
        mode='lines',
        name='Mohr-Coulomb Circle',
        line=dict(color='blue', width=2)
    ))

    # Add the dashed line for failure line
    fig.add_trace(go.Scatter(
        x=x_line,
        y=y_line,
        mode='lines',
        name='Mobilized Shear Stress = σ\' * tan(ϕ°) + c',
        line=dict(color='red', width=2, dash='dash')
    ))

    fig.update_layout(
        title=dict(
            text="Mohr-Coulomb Circle",
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title="σ' (Effective Stress) [kN/m²]",
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True,
            autorange='reversed'
        ),
        yaxis=dict(
            title="τ (Mobilized Shear Stress) [kN/m²]",
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True,
            autorange='reversed'
        ),
        template='plotly_white',
        xaxis_scaleanchor="y"
    )
    # fig.update_layout(
    #     xaxis=dict(rangemode='tozero'),
    #     yaxis=dict(rangemode='tozero'),
    # )
    return fig

def plot_p_q(p_list, q_list):

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=p_list,
        y=q_list,
        mode='markers+lines',
        name='p\' vs q',
        marker=dict(size=8, color='blue')
    ))

    fig.update_layout(
        title=dict(text="Mean Effective Stress vs Deviatoric Stress", x=0.5),
        xaxis=dict(
            title="p' (Mean Effective Stress) [kN/m²]",
            autorange='reversed',
            showline=True,
            mirror=True,
            linecolor='black'
        ),
        yaxis=dict(
            title="q (Deviatoric Stress) [kN/m²]",
            showline=True,
            mirror=True,
            linecolor='black'
        ),
        template='plotly_white'
    )
    fig.update_layout(
        xaxis=dict(rangemode='tozero'),
        yaxis=dict(rangemode='tozero'),
    )
    return fig

def plot_volumetric_strain(vertical_strain, volumetric_strain):
    """
    Plots volumetric strain against vertical strain.

    Args:
        vertical_strain (list): List of vertical strains.
        volumetric_strain (list): List of volumetric strains.
    """
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=vertical_strain,
        y=volumetric_strain,
        mode='markers',
        marker=dict(size=10, color='blue'),
        name='Volumetric Strain'
    ))

    fig.update_layout(
        title=dict(
            text='Volumetric Strain vs Vertical Strain Plot',
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title='Vertical Strain [−]',
            autorange='reversed',
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        yaxis=dict(
            title='Volumetric Strain [−]',
            autorange='reversed',
            showline=True,
            linewidth=2,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            tickcolor='black',
            ticklen=5,
            mirror=True
        ),
        template='plotly_white',
    )
    fig.update_layout(
        xaxis=dict(rangemode='tozero'),
        yaxis=dict(rangemode='tozero'),
    )
    return fig

def lab_test(dll_path, index, umat_parameters, time_step, maximum_strain, initial_effective_cell_pressure):
    json_file_path = 'test_triaxial/MaterialParameters_stage1.json'
    material_editor = MaterialEditor(json_file_path)
    material_editor._update_material_and_save({"UMAT_PARAMETERS": umat_parameters,
                                               "UDSM_NAME": dll_path,
                                               "UDSM_NUMBER": index})

    project_param_path = 'test_triaxial/ProjectParameters_stage2.json'
    project_editor = ProjectParameterEditor(project_param_path)
    project_editor.update_time_step_properties(time_step)

    mdpa_path = 'test_triaxial/triaxial.mdpa'
    mdpa_editor = MdpaEditor(mdpa_path)
    mdpa_editor.update_maximum_strain(maximum_strain)
    mdpa_editor.update_initial_effective_cell_pressure(initial_effective_cell_pressure)


# List of output files to process
    output_files = [
        os.path.join('gid_output', "triaxial_Stage_2.post.res")]

    reshaped_values_by_time, vertical_strain, volumetric_strain = run_triaxial_test(output_files)

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

    p_list = []
    q_list = []

    for eigenvalues in eigenvalue_list:
        sigma1, sigma2, sigma3 = np.sort(eigenvalues)
        p = (sigma1 + sigma2 + sigma3) / 3
        q = np.sqrt(0.5 * ((sigma1 - sigma2) ** 2 + (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2))

        p_list.append(p)
        q_list.append(q)
    sigma_diff = abs(np.array(sigma1_list) - np.array(sigma3_list))

    # plot_sigma(sigma1_list, sigma3_list)
    # plot_delta_sigma(vertical_strain, sigma_diff)
    # plot_mohr_coulomb_circle(sigma1_list[-1], sigma3_list[-1])
    # plot_p_q(p_list, q_list)
    # plot_volumetric_strain(vertical_strain, volumetric_strain)

    figs = [
        plot_delta_sigma(vertical_strain, sigma_diff),
        plot_volumetric_strain(vertical_strain, volumetric_strain),
        plot_sigma(sigma1_list, sigma3_list),
        plot_p_q(p_list, q_list),
        plot_mohr_coulomb_circle(sigma1_list[-1], sigma3_list[-1]),
    ]
    return figs

if __name__ == "__main__":
    figs = lab_test("../MohrCoulomb64.dll", ["10000", "0.3", "0.0", "30.0", "0.0", "0.0", "1.0", "0.3"],
             1, 0.01, 20, 100)
