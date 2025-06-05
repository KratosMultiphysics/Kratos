import tkinter as tk
from tkinter import ttk, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading

from set_triaxial_test import run_triaxial_simulation
from ui_plot_manager import render_plots
from ui_logger import init_log_widget, log_message, clear_log
from ui_udsm_parser import input_parameters_format_to_unicode
import traceback


MAX_STRAIN_LABEL = "Maximum Strain |εᵧᵧ|"
INIT_PRESSURE_LABEL = "Initial effective cell pressure |σ'ₓₓ|"
STRESS_INC_LABEL = "Stress increment |σ'ᵧᵧ|"
NUM_STEPS_LABEL = "Number of steps"
DURATION_LABEL = "Duration"
FL2_UNIT_LABEL = "kN/m²"
SECONDS_UNIT_LABEL = "s"
PERCENTAGE_UNIT_LABEL = "%"
WITHOUT_UNIT_LABEL = ""


class GeotechTestUI:
    def __init__(self, root, parent_frame, test_name, dll_path, model_dict, external_widgets=None):
        self.root = root
        self.parent = parent_frame
        self.test_name = test_name
        self.dll_path = dll_path
        self.model_dict = model_dict

        self.model_var = tk.StringVar(root)
        self.model_var.set(model_dict["model_name"][0])
        self.current_test = tk.StringVar(value="Triaxial")

        self._init_frames()
        self._init_dropdown_section()
        self._init_plot_canvas()
        self._create_input_fields()

        self.is_running = False
        self.external_widgets = external_widgets if external_widgets else []

    def _start_simulation_thread(self):
        if self.is_running:
            return
        self.is_running = True
        self._freeze_gui()
        threading.Thread(target=self._run_simulation, daemon=True).start()

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

        self.test_selector_frame = ttk.Frame(self.param_frame, padding="5")
        self.test_selector_frame.pack(fill="x", pady=(10, 5))

        self.test_buttons = {}
        for test_name in ["Triaxial", "Oedometer", "Direct Shear"]:
            btn = tk.Button(
                self.test_selector_frame,
                text=test_name,
                font=("Arial", 8, "bold"),
                width=10,
                height=10,
                relief="raised",
                command=lambda name=test_name: self._switch_test(name)
            )
            btn.pack(side="left", padx=5, pady=5)
            self.test_buttons[test_name] = btn

        self.test_input_frame = ttk.Frame(self.param_frame, padding="10")
        self.test_input_frame.pack(fill="both", expand=True)

        self._switch_test("Triaxial")

        self.run_button = ttk.Button(self.button_frame, text="Run Calculation", command=self._start_simulation_thread)
        self.run_button.pack(pady=5)

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
        self.mohr_frame = ttk.Frame(self.param_frame)
        self.mohr_frame.pack(fill="x", padx=10, pady=5)

        self.mohr_checkbox_widget = ttk.Checkbutton(
            self.mohr_frame,
            text="Mohr-Coulomb Model",
            variable=self.mohr_checkbox,
            command=self._toggle_mohr_options
        )
        self.mohr_checkbox_widget.pack(side="left")

        self.c_label = ttk.Label(self.mohr_frame, text="Cohesion Index (1-based)")
        self.c_dropdown = ttk.Combobox(self.mohr_frame, textvariable=self.cohesion_var,
                                       values=[str(i+1) for i in range(len(params))], state="readonly", width=10)

        self.phi_label = ttk.Label(self.mohr_frame, text="Friction Angle Index (1-based)")
        self.phi_dropdown = ttk.Combobox(self.mohr_frame, textvariable=self.phi_var,
                                         values=[str(i+1) for i in range(len(params))], state="readonly", width=10)

    def _toggle_mohr_options(self):
        widgets = [self.c_label, self.c_dropdown, self.phi_label, self.phi_dropdown]
        if self.mohr_checkbox.get():
            log_message("Mohr-Coulomb model is selected.", "info")
            for w in widgets:
                w.pack(side="left", padx=5)
        else:
            for w in widgets:
                w.pack_forget()

    def _switch_test(self, test_name):
        clear_log()
        self.current_test.set(test_name)

        for w in self.test_input_frame.winfo_children():
            w.destroy()

        if test_name == "Triaxial":
            self.triaxial_widgets = self._create_entries(
                self.test_input_frame,
                "Triaxial Input Data",
                [INIT_PRESSURE_LABEL, MAX_STRAIN_LABEL, NUM_STEPS_LABEL, DURATION_LABEL],
                [FL2_UNIT_LABEL, PERCENTAGE_UNIT_LABEL, WITHOUT_UNIT_LABEL, SECONDS_UNIT_LABEL],
                {INIT_PRESSURE_LABEL: "100", MAX_STRAIN_LABEL: "10",
                 NUM_STEPS_LABEL: "100", DURATION_LABEL: "1.0"}
            )
        elif test_name == "Oedometer":
            self.oedometer_widgets = self._create_entries(
                self.test_input_frame,
                "Oedometer Input Data",
                [DURATION_LABEL, STRESS_INC_LABEL, NUM_STEPS_LABEL],
                [FL2_UNIT_LABEL, FL2_UNIT_LABEL, WITHOUT_UNIT_LABEL],
                {DURATION_LABEL: "1.0", STRESS_INC_LABEL: "100", NUM_STEPS_LABEL: "100"}
            )
        elif test_name == "Direct Shear":
            self.shear_widgets = self._create_entries(
                self.test_input_frame,
                "Direct Shear Input Data",
                [INIT_PRESSURE_LABEL, MAX_STRAIN_LABEL, NUM_STEPS_LABEL, DURATION_LABEL],
                [FL2_UNIT_LABEL, PERCENTAGE_UNIT_LABEL, "", SECONDS_UNIT_LABEL],
                {INIT_PRESSURE_LABEL: "100", MAX_STRAIN_LABEL: "10",
                 NUM_STEPS_LABEL: "100", DURATION_LABEL: "1.0"}
            )

        log_message(f"{test_name} test selected.", "info")

    def _run_simulation(self):
        clear_log()
        try:
            log_message("Starting calculation... Please wait...", "info")
            log_message("Validating input...", "info")
            self.root.update_idletasks()

            umat_params = [e.get() for e in self.entry_widgets.values()]

            if self.current_test.get() != "Triaxial":
                raise NotImplementedError(f"{self.current_test.get()} simulation not yet implemented.")

            eps_max = float(self.triaxial_widgets[MAX_STRAIN_LABEL].get())
            sigma_init = float(self.triaxial_widgets[INIT_PRESSURE_LABEL].get())
            n_steps = float(self.triaxial_widgets[NUM_STEPS_LABEL].get())
            duration = float(self.triaxial_widgets[DURATION_LABEL].get())

            if any(val <= 0 for val in [eps_max, n_steps, duration]) or sigma_init < 0:
                raise ValueError("All values must be positive and non-zero.")

            cohesion_phi_indices = None
            if self.mohr_checkbox.get():
                cohesion_phi_indices = (int(self.cohesion_var.get()), int(self.phi_var.get()))

            index = self.model_dict["model_name"].index(self.model_var.get()) + 1 if self.dll_path else -2

            log_message("Calculating...", "info")
            self.root.update_idletasks()

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
            self.root.after(0, lambda: render_plots(figs, self.axes, self.canvas))
            log_message("Simulation completed successfully.", "info")

        except Exception:
            log_message("An error occurred during simulation:", "error")
            log_message(traceback.format_exc(), "error")

        finally:
            self.root.after(0, self._unfreeze_gui)
            self.is_running = False

    def _enable_run_button(self):
            self.run_button.config(state="normal")
            self.is_running = False

    def _set_widget_state(self, parent, state):
        for child in parent.winfo_children():
            if isinstance(child, (ttk.Entry, ttk.Combobox, tk.Button, ttk.Button, tk.Checkbutton, ttk.Checkbutton)):
                child.configure(state=state)
            elif isinstance(child, scrolledtext.ScrolledText):
                child.config(state=state if state == "normal" else "disabled")
            elif isinstance(child, (ttk.Frame, tk.Frame)):
                self._set_widget_state(child, state)

        for widget in self.external_widgets:
            widget.configure(state=state)

    def _freeze_gui(self):
        self._set_widget_state(self.left_frame, "disabled")

    def _unfreeze_gui(self):
        self._set_widget_state(self.left_frame, "normal")
        self.run_button.config(state="normal")
