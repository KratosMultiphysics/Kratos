import tkinter as tk
from tkinter import ttk, scrolledtext
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from set_triaxial_test import run_triaxial_simulation
from ui_plot_manager import render_plots
from ui_logger import init_log_widget, log_message, clear_log
import traceback
from ui_udsm_parser import input_parameters_format_to_unicode

def build_ui_from_model(root, parent_frame, dll_path, model_dict):
    global fig, axes, canvas

    def create_soil_input_fields(param_frame, params, units, default_values):
        entry_widgets = {}
        ttk.Label(param_frame, text="Soil Input Parameters", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
        for i, param in enumerate(params):
            unit = units[i] if i < len(units) else ""
            row = ttk.Frame(param_frame)
            row.pack(fill="x", padx=10, pady=2)
            ttk.Label(row, text=param).pack(side="left", padx=5)
            entry = ttk.Entry(row)
            entry.insert(0, default_values.get(param, "Please enter a number"))
            entry.pack(side="left", fill="x", expand=True)
            ttk.Label(row, text=unit).pack(side="left", padx=5)
            entry_widgets[param] = entry
        return entry_widgets

    def create_triaxial_fields(param_frame):
        triaxial_entry_widgets = {}
        frame = ttk.Frame(param_frame, padding="10")
        frame.pack(fill="x", pady=10)
        ttk.Label(frame, text="Triaxial Input Data", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
        for label, unit, default in [
            ("Initial effective cell pressure |σ'₃₃|", "kN/m²", "100"),
            ("Maximum Strain |εᵧᵧ|", "%", "10"),
            ("Number of steps", "", "100"),
            ("Duration", "s", "1.0")
        ]:
            row = ttk.Frame(frame)
            row.pack(fill="x", padx=10, pady=2)
            ttk.Label(row, text=label).pack(side="left", padx=5)
            entry = ttk.Entry(row)
            entry.insert(0, default)
            entry.pack(side="left", fill="x", expand=True, padx=5)
            ttk.Label(row, text=unit).pack(side="left", padx=5)
            triaxial_entry_widgets[label] = entry
        return triaxial_entry_widgets

    def gather_and_validate(entry_widgets, triaxial_widgets):
        try:
            umat_params = [e.get() for e in entry_widgets.values()]
            eps_max = float(triaxial_widgets["Maximum Strain |εᵧᵧ|"].get())
            sigma_init = float(triaxial_widgets["Initial effective cell pressure |σ'₃₃|"].get())
            n_steps = float(triaxial_widgets["Number of steps"].get())
            duration = float(triaxial_widgets["Duration"].get())

            if n_steps <= 0 or duration <= 0 or eps_max <= 0 or sigma_init < 0:
                raise ValueError in log_message("All values must be positive and non-zero.", "error")

            return umat_params, eps_max, sigma_init, n_steps, duration
        except ValueError as error:
            # messagebox.showerror("Invalid Input", str(e))
            log_message(f"Input error: {error}", "error")
            return None

    def run_simulation(entry_widgets, triaxial_widgets, model_var, mohr_checkbox, cohesion_var, phi_var):
        clear_log()
        log_message("Starting triaxial calculation...", "info")
        log_message("Validating input...", "info")
        root.update_idletasks()

        result = gather_and_validate(entry_widgets, triaxial_widgets)
        if not result:
            return

        umat_params, eps_max, sigma_init, n_steps, duration = result
        log_message("Calculating...", "info")
        root.update_idletasks()

        if mohr_checkbox.get():
            c_idx = int(cohesion_var.get())
            phi_idx = int(phi_var.get())
            cohesion_phi_indices = (c_idx, phi_idx)
        else:
            cohesion_phi_indices = None

        try:
            if dll_path:
                index = model_dict["model_name"].index(model_var.get()) + 1
                figs = run_triaxial_simulation(dll_path, index, umat_params, n_steps, duration, eps_max, sigma_init,
                                               cohesion_phi_indices=cohesion_phi_indices)
            else:
                figs = run_triaxial_simulation(
                    dll_path="",
                    index=-2,
                    umat_parameters=[float(umat_params[0]), float(umat_params[1])],
                    num_steps=n_steps,
                    end_time=duration,
                    maximum_strain=eps_max,
                    initial_effective_cell_pressure=sigma_init,
                    cohesion_phi_indices=cohesion_phi_indices
                )
            render_plots(figs, axes, canvas)
            log_message("Simulation completed successfully.", "info")

        except Exception as e:
            traceback_str = traceback.format_exc()
            log_message("An error occurred during simulation:", "error")
            log_message(traceback_str, "error")

    def _create_input_fields(model_var, param_frame, button_frame):
        for w in param_frame.winfo_children() + button_frame.winfo_children():
            w.destroy()

        raw_defaults  = {
            "1. E": "10000", "2. n_ur": "0.3", "3. c'": "0.0", "4. f_peak": "30.0",
            "5. y_peak": "0.0", "6. s_t, cut-off": "0.0", "7. yield function (MC=1 DP=2 MNC=3 MN=4)": "1",
            "8. n_un (UMAT)": "0.3", "YOUNG_MODULUS": "10000", "POISSON_RATIO": "0.3"
        }

        selected_model = model_var.get()
        index = model_dict["model_name"].index(selected_model)
        params = model_dict["param_names"][index]
        units = model_dict.get("param_units", [[]])[index]

        default_values = {
            input_parameters_format_to_unicode(k): v for k, v in raw_defaults.items()
        }
        entry_widgets = create_soil_input_fields(param_frame, params, units, default_values)
        triaxial_widgets = create_triaxial_fields(param_frame)

        ttk.Button(button_frame, text="Run Calculation",
                   command=lambda: run_simulation(entry_widgets, triaxial_widgets, model_var, mohr_checkbox,
                                                  cohesion_var, phi_var)).pack(pady=5)

        mohr_checkbox = tk.BooleanVar()
        cohesion_var = tk.StringVar(value="3")
        phi_var = tk.StringVar(value="4")

        def toggle_mohr_options():
            if mohr_checkbox.get():
                c_label.pack()
                c_dropdown.pack()
                phi_label.pack()
                phi_dropdown.pack()
            else:
                c_label.pack_forget()
                c_dropdown.pack_forget()
                phi_label.pack_forget()
                phi_dropdown.pack_forget()

        ttk.Checkbutton(param_frame, text="Mohr-Coulomb Model", variable=mohr_checkbox,
                        command=toggle_mohr_options).pack(anchor="w", padx=10, pady=(10, 5))

        c_label = ttk.Label(param_frame, text="Cohesion Index (1-based)")
        c_dropdown = ttk.Combobox(param_frame, textvariable=cohesion_var,
                                  values=[str(i+1) for i in range(len(params))], state="readonly")

        phi_label = ttk.Label(param_frame, text="Friction Angle Index (1-based)")
        phi_dropdown = ttk.Combobox(param_frame, textvariable=phi_var,
                                    values=[str(i+1) for i in range(len(params))], state="readonly")


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

    log_frame = ttk.Frame(dropdown_frame, padding="5")
    log_frame.pack(fill="x", padx=10, pady=(0, 10))

    ttk.Label(log_frame, text="Log Output:", font=("Arial", 10, "bold")).pack(anchor="w")
    log_widget = scrolledtext.ScrolledText(log_frame, height=6, width=40, state="disabled", wrap="word", font=("Courier", 9))
    log_widget.pack(fill="x", expand=False)
    init_log_widget(log_widget)


    plot_frame = ttk.Frame(parent_frame, padding="5", width=800, height=600)
    plot_frame.pack(side="right", fill="both", expand=True, padx=5, pady=5)

    fig = plt.figure(figsize=(12, 15))
    gs = GridSpec(3, 2, figure=fig, wspace=0.4, hspace=0.6)
    axes = [fig.add_subplot(gs[i]) for i in range(5)]
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.get_tk_widget().pack(fill="both", expand=True)
    if len(fig.axes) == 6:
        fig.delaxes(fig.axes[5])

    model_var.trace("w", lambda *args: _create_input_fields(model_var, param_frame, button_frame))
    _create_input_fields(model_var, param_frame, button_frame)
