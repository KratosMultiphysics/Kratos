import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from material_editor import MaterialEditor
from set_triaxial_test import run_triaxial_simulation
from ui_plot_manager import render_plots

def build_ui_from_model(root, parent_frame, dll_path, model_dict):
    global fig, axes, canvas

    def _create_input_fields(model_var, param_frame, button_frame, plot_frame):
        for widget in param_frame.winfo_children():
            widget.destroy()
        for widget in button_frame.winfo_children():
            widget.destroy()

        entry_widgets = {}
        triaxial_entry_widgets = {}

        default_values = {
            "1. E": "10000", "2. n_ur": "0.3", "3. c'": "0.0", "4. f_peak": "30.0",
            "5. y_peak": "0.0", "6. s_t, cut-off": "0.0", "7. yield function (MC=1 DP=2 MNC=3 MN=4)": "1",
            "8. n_un (UMAT)": "0.3", "YOUNG_MODULUS": "10000", "POISSON_RATIO": "0.3"
        }

        selected_model = model_var.get()
        index = model_dict["model_name"].index(selected_model)
        params = model_dict["param_names"][index]

        ttk.Label(param_frame, text="Soil Input Parameters", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
        for param in params:
            row = ttk.Frame(param_frame)
            row.pack(fill="x", padx=10, pady=2)
            ttk.Label(row, text=param).pack(side="left", padx=5)
            entry = ttk.Entry(row)
            entry.insert(0, default_values.get(param, "N/A"))
            entry.pack(side="left", fill="x", expand=True)
            entry_widgets[param] = entry

        triaxial_frame = ttk.Frame(param_frame, padding="10")
        triaxial_frame.pack(fill="x", pady=10)

        ttk.Label(triaxial_frame, text="Triaxial Input Data", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
        for label, unit, default in [
            ("Initial effective cell pressure |σ'₃₃|", "kN/m²", "100"),
            ("Maximum Strain |εᵧᵧ|", "%", "10"),
            ("Number of steps", "", "100"),
            ("Duration", "s", "1.0")
        ]:
            row = ttk.Frame(triaxial_frame)
            row.pack(fill="x", padx=10, pady=2)
            ttk.Label(row, text=label).pack(side="left", padx=5)
            entry = ttk.Entry(row)
            entry.insert(0, default)
            entry.pack(side="left", fill="x", expand=True, padx=5)
            ttk.Label(row, text=unit).pack(side="left", padx=5)
            triaxial_entry_widgets[label] = entry

        def _gather_and_validate_inputs():
            try:
                umat_params = [e.get() for e in entry_widgets.values()]
                eps_max = float(triaxial_entry_widgets["Maximum Strain |εᵧᵧ|"].get())
                sigma_init = float(triaxial_entry_widgets["Initial effective cell pressure |σ'₃₃|"].get())
                n_steps = float(triaxial_entry_widgets["Number of steps"].get())
                duration = float(triaxial_entry_widgets["Duration"].get())

                if n_steps <= 0:
                    raise ValueError("Number of steps must be greater than zero.")
                if duration <= 0:
                    raise ValueError("Duration must be greater than zero.")
                if eps_max <= 0:
                    raise ValueError("Maximum strain must be greater than zero.")
                if sigma_init < 0:
                    raise ValueError("Initial effective cell pressure cannot be negative.")

                return umat_params, eps_max, sigma_init, n_steps, duration
            except ValueError as e:
                messagebox.showerror("Invalid Input", str(e))
                return None

        def _run_simulation():
            result = _gather_and_validate_inputs()
            if not result:
                return
            umat_params, eps_max, sigma_init, n_steps, duration = result

            if dll_path:
                index = model_dict["model_name"].index(model_var.get()) + 1
                figs = run_triaxial_simulation(dll_path, index, umat_params, n_steps, duration, eps_max, sigma_init)
            else:
                editor = MaterialEditor("test_triaxial/MaterialParameters.json")
                entries = {
                    "YOUNG_MODULUS": umat_params[0],
                    "POISSON_RATIO": umat_params[1]
                }
                editor._update_material_and_save(entries)
                figs = run_triaxial_simulation("", -2, [], n_steps, duration, eps_max, sigma_init)

            render_plots(figs, axes, canvas)

        ttk.Button(button_frame, text="Run Calculation", command=_run_simulation).pack(pady=5)

    model_var = tk.StringVar(root)
    model_var.set(model_dict["model_name"][0])

    dropdown_frame = ttk.Frame(parent_frame, padding="10", width=700, height=100)
    dropdown_frame.pack_propagate(False)
    dropdown_frame.pack(side="left", fill="y", padx=10, pady=10)

    ttk.Label(dropdown_frame, text="Select a Model:", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
    model_menu = ttk.Combobox(dropdown_frame, textvariable=model_var, values=model_dict["model_name"], state="readonly")
    model_menu.pack(side="top", fill="x", expand=True, padx=5)

    param_frame = ttk.Frame(dropdown_frame, padding="10")
    param_frame.pack(fill="both", expand=True, pady=10)

    button_frame = ttk.Frame(dropdown_frame, padding="10")
    button_frame.pack(fill="x", pady=10)

    plot_frame = ttk.Frame(parent_frame, padding="5", width=800, height=600)
    plot_frame.pack(side="right", fill="both", expand=True, padx=5, pady=5)

    width_per_ax, height_per_ax = 6, 5
    nrows, ncols = 3, 2
    fig = plt.figure(figsize=(ncols * width_per_ax, nrows * height_per_ax))
    gs = GridSpec(nrows, ncols, figure=fig, wspace=0.4, hspace=0.6)
    axes = [fig.add_subplot(gs[i]) for i in range(5)]
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.get_tk_widget().pack(fill="both", expand=True)
    if len(fig.axes) == 6:
        fig.delaxes(fig.axes[5])

    model_var.trace("w", lambda *args: _create_input_fields(model_var, param_frame, button_frame, plot_frame))
    _create_input_fields(model_var, param_frame, button_frame, plot_frame)
