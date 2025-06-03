import tkinter as tk
from tkinter import ttk, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from set_triaxial_test import run_triaxial_simulation
from ui_plot_manager import render_plots
from ui_logger import init_log_widget, log_message, clear_log
from ui_udsm_parser import input_parameters_format_to_unicode
import traceback

class GeotechTestUI:
    def __init__(self, root, parent_frame, test_name, dll_path, model_dict):
        self.root = root
        self.parent = parent_frame
        self.test_name = test_name
        self.dll_path = dll_path
        self.model_dict = model_dict

        self.model_var = tk.StringVar(root)
        self.model_var.set(model_dict["model_name"][0])

        self._init_frames()
        self._init_dropdown_section()
        self._init_plot_canvas()
        self._create_input_fields()

    def _init_frames(self):
        self.left_frame = ttk.Frame(self.parent, padding="10", width=700)
        self.left_frame.pack_propagate(False)
        self.left_frame.pack(side="left", fill="y", padx=10, pady=10)

        self.dropdown_frame = ttk.Frame(self.left_frame)
        self.dropdown_frame.pack(fill="x")

        self.param_frame = ttk.Frame(self.left_frame, padding="10")
        self.param_frame.pack(fill="both", expand=True, pady=10)

        self.button_frame = ttk.Frame(self.left_frame, padding="10")
        self.button_frame.pack(fill="x", pady=10)

        self.log_frame = ttk.Frame(self.left_frame, padding="5")
        self.log_frame.pack(fill="x", padx=10, pady=(0, 10))

    def _init_plot_canvas(self):
        self.plot_frame = ttk.Frame(self.parent, padding="5", width=800, height=600)
        self.plot_frame.pack(side="right", fill="both", expand=True, padx=5, pady=5)

        self.fig = plt.figure(figsize=(12, 15))
        self.gs = GridSpec(3, 2, figure=self.fig, wspace=0.4, hspace=0.6)
        self.axes = [self.fig.add_subplot(self.gs[i]) for i in range(5)]
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        if len(self.fig.axes) == 6:
            self.fig.delaxes(self.fig.axes[5])

    def _init_dropdown_section(self):
        ttk.Label(self.dropdown_frame, text="Select a Model:", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
        model_menu = ttk.Combobox(self.dropdown_frame, textvariable=self.model_var, values=self.model_dict["model_name"], state="readonly")
        model_menu.pack(side="top", fill="x", expand=True, padx=5)
        self.model_var.trace("w", lambda *args: self._create_input_fields())

        ttk.Label(self.log_frame, text="Log Output:", font=("Arial", 10, "bold")).pack(anchor="w")
        self.log_widget = scrolledtext.ScrolledText(self.log_frame, height=6, width=40, state="disabled", wrap="word", font=("Courier", 9))
        self.log_widget.pack(fill="x", expand=False)
        init_log_widget(self.log_widget)

    def _create_input_fields(self):
        for w in self.param_frame.winfo_children() + self.button_frame.winfo_children():
            w.destroy()

        index = self.model_dict["model_name"].index(self.model_var.get())
        params = self.model_dict["param_names"][index]
        units = self.model_dict.get("param_units", [[]])[index]

        raw_defaults = {
            "1. E": "10000", "2. n_ur": "0.3", "3. c'": "0.0", "4. f_peak": "30.0",
            "5. y_peak": "0.0", "6. s_t, cut-off": "0.0", "7. yield function (MC=1 DP=2 MNC=3 MN=4)": "1",
            "8. n_un (UMAT)": "0.3", "YOUNG_MODULUS": "10000", "POISSON_RATIO": "0.3"
        }
        default_values = {
            input_parameters_format_to_unicode(k): v for k, v in raw_defaults.items()
        }

        self.entry_widgets = self._create_entries(self.param_frame, "Soil Input Parameters", params, units, default_values)

        self.mohr_checkbox = tk.BooleanVar()
        self.cohesion_var = tk.StringVar(value="3")
        self.phi_var = tk.StringVar(value="4")

        self._create_mohr_options(params)

        self.triaxial_widgets = self._create_entries(self.param_frame, "Triaxial Input Data", [
            "Initial effective cell pressure |σ'₃₃|", "Maximum Strain |εᵧᵧ|", "Number of steps", "Duration"
        ], ["kN/m²", "%", "", "s"], {k: v for k, v in zip([
            "Initial effective cell pressure |σ'₃₃|", "Maximum Strain |εᵧᵧ|", "Number of steps", "Duration"
        ], ["100", "10", "100", "1.0"])})

        ttk.Button(self.button_frame, text="Run Calculation", command=self._run_simulation).pack(pady=5)

    def _create_entries(self, frame, title, labels, units, defaults):
        widgets = {}
        ttk.Label(frame, text=title, font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
        for i, label in enumerate(labels):
            unit = units[i] if i < len(units) else ""
            row = ttk.Frame(frame)
            row.pack(fill="x", padx=10, pady=2)
            ttk.Label(row, text=label).pack(side="left", padx=5)
            entry = ttk.Entry(row)
            entry.insert(0, defaults.get(label, ""))
            entry.pack(side="left", fill="x", expand=True)
            ttk.Label(row, text=unit).pack(side="left", padx=5)
            widgets[label] = entry
        return widgets

    def _create_mohr_options(self, params):
        mohr_row = ttk.Frame(self.param_frame)
        mohr_row.pack(fill="x", padx=10, pady=5)

        checkbox = ttk.Checkbutton(mohr_row, text="Mohr-Coulomb Model", variable=self.mohr_checkbox,
                                   command=self._toggle_mohr_options)
        checkbox.pack(side="left")

        self.c_label = ttk.Label(mohr_row, text="Cohesion Index (1-based)")
        self.c_dropdown = ttk.Combobox(mohr_row, textvariable=self.cohesion_var,
                                       values=[str(i+1) for i in range(len(params))], state="readonly", width=10)

        self.phi_label = ttk.Label(mohr_row, text="Friction Angle Index (1-based)")
        self.phi_dropdown = ttk.Combobox(mohr_row, textvariable=self.phi_var,
                                         values=[str(i+1) for i in range(len(params))], state="readonly", width=10)

    def _toggle_mohr_options(self):
        widgets = [self.c_label, self.c_dropdown, self.phi_label, self.phi_dropdown]
        if self.mohr_checkbox.get():
            for w in widgets:
                w.pack(side="left", padx=5)
        else:
            for w in widgets:
                w.pack_forget()

    def _run_simulation(self):
        clear_log()
        log_message("Starting triaxial calculation...", "info")
        self.root.update_idletasks()

        try:
            umat_params = [e.get() for e in self.entry_widgets.values()]
            eps_max = float(self.triaxial_widgets["Maximum Strain |εᵧᵧ|"].get())
            sigma_init = float(self.triaxial_widgets["Initial effective cell pressure |σ'₃₃|"].get())
            n_steps = float(self.triaxial_widgets["Number of steps"].get())
            duration = float(self.triaxial_widgets["Duration"].get())

            if any(val <= 0 for val in [eps_max, n_steps, duration]) or sigma_init < 0:
                raise ValueError("All values must be positive and non-zero.")

            cohesion_phi_indices = None
            if self.mohr_checkbox.get():
                cohesion_phi_indices = (int(self.cohesion_var.get()), int(self.phi_var.get()))

            index = self.model_dict["model_name"].index(self.model_var.get()) + 1 if self.dll_path else -2
            figs = run_triaxial_simulation(
                dll_path=self.dll_path or "",
                index=index,
                umat_parameters=[float(x) for x in umat_params],
                num_steps=n_steps,
                end_time=duration,
                maximum_strain=eps_max,
                initial_effective_cell_pressure=sigma_init,
                cohesion_phi_indices=cohesion_phi_indices
            )
            render_plots(figs, self.axes, self.canvas)
            log_message("Simulation completed successfully.", "info")

        except Exception as e:
            log_message("An error occurred during simulation:", "error")
            log_message(traceback.format_exc(), "error")
