import os
import json
import numpy as np
import plotly.graph_objects as go
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

    reshaped_values_by_time = {}

    # Process Cauchy stress results
    for idx, result in enumerate(cauchy_stress_results):
        time_step = result["time"]
        element_values = result["values"]

        if not element_values:  # Skip empty lists
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

    # Process displacement results
    node_2_displacements = []
    for result in displacement_results:
        element_values = result["values"]

        # Find the entry for node:2
        for element in element_values:
            if element["node"] == 2:
                # Extract the second array (index 1) from the value
                if len(element["value"]) > 1:  # Ensure there is a second array
                    node_2_displacements.append(element["value"][1])
                break  # Exit loop once node:2 is found

    return reshaped_values_by_time, node_2_displacements

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
    fig.show()

def plot_delta_sigma(displacement, diff):
    """
    Plots the difference between σ₁ and σ₃ against displacement.

    Args:
        displacement (list): List of displacements.
        diff (list): List of differences between σ₁ and σ₃.
    """
    fig = go.Figure()

    # Add scatter plot for σ₁ vs σ₃
    fig.add_trace(go.Scatter(
        x=displacement,
        y=diff,
        mode='markers',
        marker=dict(size=10, color='blue'),
        name='σ₁ vs σ₃'
    ))

    fig.update_layout(
        title=dict(
            text='|σ₁ - σ₃| vs displacement Plot',
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title='Displacement [m]',
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
    fig.show()

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
            # range=[0, center + 2.5 * radius],
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
            # range=[0, 1.5 * radius],
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
    fig.show()

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
    fig.show()


class MaterialEditorApp:
    def __init__(self, master, json_path):
        self.master = master
        self.master.title("Material Parameters Editor")
        self.master.geometry("800x900")  # Bigger window
        self.json_path = json_path
        self.entries = {}

        self._load_json()
        self._build_gui()

    def _load_json(self):
        with open(self.json_path, 'r') as f:
            self.data = json.load(f)
        self.variables = self.data["properties"][0]["Material"]["Variables"]

    def _build_gui(self):
        canvas = tk.Canvas(self.master)
        scrollbar = tk.Scrollbar(self.master, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        row = 0

        # Helper frame for UMAT info
        helper_text = (
            "UMAT_PARAMETERS Mapping for MohrCoulomb64.dll:\n"
            "  1 : E        → Young's modulus\n"
            "  2 : Nu(ν)    → Poisson's ratio (unloading/reloading)\n"
            "  3 : C        → Cohesion\n"
            "  4 : Phi      → Friction angle (°)\n"
            "  5 : Psi      → Dilation angle (°)\n"
            "  6 : Tens     → Allowable tensile stress\n"
            "  7 : Yield    → Yield function index (1 = Mohr Coulomb)\n"
            "  8 : Nu_undr  → Undrained Poisson's ratio"
        )

        help_frame = tk.LabelFrame(scrollable_frame, text="UMAT_PARAMETERS Help", padx=10, pady=5)
        help_label = tk.Label(help_frame, text=helper_text, justify="left", font=("Courier", 10))
        help_label.pack(anchor="w")
        help_frame.grid(row=row, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        row += 1

        tk.Label(scrollable_frame, text="Edit Parameters Below", font=("Arial", 14, "bold")).grid(
            row=row, column=0, columnspan=2, pady=5)
        row += 1

        for key, value in self.variables.items():
            tk.Label(scrollable_frame, text=key).grid(row=row, column=0, sticky="w", padx=10, pady=5)

            entry = tk.Entry(scrollable_frame, width=80 if isinstance(value, list) else 30)
            entry.insert(0, ', '.join(map(str, value)) if isinstance(value, list) else str(value))
            entry.grid(row=row, column=1, padx=10, pady=5)
            self.entries[key] = entry
            row += 1

        tk.Button(scrollable_frame, text="Save Changes", command=self._save_changes,
                  bg="green", fg="white", padx=10, pady=5).grid(row=row, column=0, columnspan=2, pady=20)

    def _save_changes(self):
        for key, entry in self.entries.items():
            value_str = entry.get().strip()
            if value_str.lower() == "true":
                value = True
            elif value_str.lower() == "false":
                value = False
            elif "," in value_str:
                try:
                    value = [self._convert_type(v) for v in value_str.split(",")]
                except ValueError:
                    messagebox.showerror("Error", f"Invalid list for {key}")
                    return
            else:
                try:
                    value = self._convert_type(value_str)
                except ValueError:
                    messagebox.showerror("Error", f"Invalid value for {key}")
                    return

            self.variables[key] = value

        with open(self.json_path, 'w') as f:
            json.dump(self.data, f, indent=4)

        messagebox.showinfo("Success", "Material parameters updated successfully!")
        self.master.destroy()

    def _convert_type(self, value):
        try:
            if '.' in value or 'e' in value.lower():
                return float(value)
            else:
                return int(value)
        except ValueError:
            return value


def launch_material_editor(json_path=None):
    root = tk.Tk()

    if not json_path:
        json_path = filedialog.askopenfilename(
            title="Select MaterialParameters JSON",
            filetypes=[("JSON files", "*.json")]
        )

    if json_path:
        app = MaterialEditorApp(root, json_path)
        root.mainloop()
    else:
        print("No file selected. Exiting.")
        exit()

if __name__ == "__main__":
    json_file_path = 'test_triaxial/MaterialParameters_stage1.json'

    launch_material_editor(json_file_path)

    # List of output files to process
    output_files = [
        os.path.join('gid_output', "triaxial_Stage_2.post.res")]

    reshaped_values_by_time, node_2_displacements = run_triaxial_test(output_files)

    displacement_list=[]
    for i in range (len(node_2_displacements)):
        displacement = np.linalg.norm(np.array(node_2_displacements[i]))
        displacement_list.append(displacement)

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


    vector_list = eigenvector_list
    value_list = eigenvalue_list
    plot_sigma(sigma1_list, sigma3_list)

    diff = abs(np.array(sigma1_list) - np.array(sigma3_list))
    plot_delta_sigma(displacement_list, diff)
    plot_mohr_coulomb_circle(sigma1_list[-1], sigma3_list[-1])
    plot_p_q(p_list, q_list)

