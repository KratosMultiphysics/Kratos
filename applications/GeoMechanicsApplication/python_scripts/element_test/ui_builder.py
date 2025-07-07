import os
import math
import tkinter as tk
from tkinter import ttk, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import threading
import traceback

from ui_runner import run_gui_builder
from ui_logger import init_log_widget, log_message, clear_log
from ui_udsm_parser import input_parameters_format_to_unicode

from ui_labels import (
    TRIAXIAL, DIRECT_SHEAR,
    MAX_STRAIN_LABEL, INIT_PRESSURE_LABEL, STRESS_INC_LABEL, NUM_STEPS_LABEL, DURATION_LABEL,
    FL2_UNIT_LABEL, SECONDS_UNIT_LABEL, PERCENTAGE_UNIT_LABEL, WITHOUT_UNIT_LABEL
)

class GeotechTestUI:
    def __init__(self, root, parent_frame, test_name, dll_path, model_dict, external_widgets=None):
        self.root = root
        self.parent = parent_frame
        self.test_name = test_name
        self.dll_path = dll_path
        self.model_dict = model_dict
        self.is_linear_elastic = model_dict["model_name"][0].lower() == "linear elastic model"

        self.model_var = tk.StringVar(root)
        self.model_var.set(model_dict["model_name"][0])
        self.current_test = tk.StringVar(value=TRIAXIAL)

        self._init_frames()

        self.plot_frame = ttk.Frame(self.parent, padding="5", width=800, height=600)
        self.plot_frame.pack(side="right", fill="both", expand=True, padx=5, pady=5)

        self.is_running = False
        self.external_widgets = external_widgets if external_widgets else []

        self._init_dropdown_section()
        self._create_input_fields()

    def _start_simulation_thread(self):
        if self.is_running:
            return
        self.is_running = True
        self._disable_gui()
        threading.Thread(target=self._run_simulation, daemon=True).start()

    def _init_frames(self):
        self.left_frame = ttk.Frame(self.parent, padding="10", width=615)
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

    def _init_plot_canvas(self, num_plots):
        self._destroy_existing_plot_canvas()

        self.fig = plt.figure(figsize=(12, 8), dpi=100)
        rows = math.ceil(math.sqrt(num_plots))
        cols = math.ceil(num_plots / rows)

        self.gs = GridSpec(rows, cols, figure=self.fig, wspace=0.4, hspace=0.6)
        self.axes = [self.fig.add_subplot(self.gs[i]) for i in range(num_plots)]
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def _init_dropdown_section(self):
        ttk.Label(self.dropdown_frame, text="Select a Model:", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
        self.model_menu = ttk.Combobox(self.dropdown_frame, textvariable=self.model_var, values=self.model_dict["model_name"], state="readonly")
        self.model_menu.pack(side="top", fill="x", expand=True, padx=5)
        self.model_var.trace("w", lambda *args: self._create_input_fields())

        if self.is_linear_elastic:
            self.model_menu.configure(state="disabled")
        else:
            self.model_menu.configure(state="readonly")

        ttk.Label(self.log_frame, text="Log Output:", font=("Arial", 10, "bold")).pack(anchor="w")
        self.log_widget = scrolledtext.ScrolledText(self.log_frame, height=6, width=40, state="disabled", wrap="word", font=("Courier", 9))
        self.log_widget.pack(fill="x", expand=False)

        self.log_widget.bind("<Key>", lambda e: "break")
        self.log_widget.bind("<Button-1>", lambda e: "break")
        self.log_widget.bind("<FocusIn>", lambda e: self.root.focus())

        init_log_widget(self.log_widget)

    def _create_input_fields(self):
        for w in self.param_frame.winfo_children() + self.button_frame.winfo_children():
            w.destroy()

        index = self.model_dict["model_name"].index(self.model_var.get())
        params = self.model_dict["param_names"][index]
        units = self.model_dict.get("param_units", [[]])[index]

        default_values = {}
        self.entry_widgets = self._create_entries(self.param_frame, "Soil Input Parameters", params, units, default_values)

        self.mohr_checkbox = tk.BooleanVar()
        self.cohesion_var = tk.StringVar(value="3")
        self.phi_var = tk.StringVar(value="4")
        self._create_mohr_options(params)
        if self.is_linear_elastic:
            self.mohr_checkbox_widget.configure(state="disabled")

        self.test_selector_frame = ttk.Frame(self.param_frame, padding="5")
        self.test_selector_frame.pack(fill="x", pady=(10, 5))

        self.test_buttons = {}
        self.test_images = {
            TRIAXIAL: tk.PhotoImage(file=os.path.join(os.path.dirname(__file__), "assets", "triaxial.png")),
            DIRECT_SHEAR: tk.PhotoImage(file=os.path.join(os.path.dirname(__file__), "assets", "direct_shear.png"))
        }

        for test_name in [TRIAXIAL, DIRECT_SHEAR]:
            btn = tk.Button(
                self.test_selector_frame,
                text=test_name,
                image=self.test_images[test_name],
                compound="top",
                font=("Arial", 10, "bold"),
                width=100,
                height=100,
                relief="raised",
                command=lambda name=test_name: self._switch_test(name)
            )
            btn.pack(side="left", padx=5, pady=5)
            self.test_buttons[test_name] = btn

        self.test_input_frame = ttk.Frame(self.param_frame, padding="10")
        self.test_input_frame.pack(fill="both", expand=True)

        self._switch_test(TRIAXIAL)

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
                                       values=[str(i+1) for i in range(len(params))], state="readonly", width=2)

        self.phi_label = ttk.Label(self.mohr_frame, text="Friction Angle Index (1-based)")
        self.phi_dropdown = ttk.Combobox(self.mohr_frame, textvariable=self.phi_var,
                                         values=[str(i+1) for i in range(len(params))], state="readonly", width=2)

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

        for name, button in self.test_buttons.items():
            if name == test_name:
                button.config(relief="sunken", bg="SystemButtonFace", state="normal")
            else:
                button.config(relief="raised", bg="SystemButtonFace", state="normal")

        for w in self.test_input_frame.winfo_children():
            w.destroy()

        if test_name == TRIAXIAL:
            self._init_plot_canvas(num_plots=5)
            ttk.Label(self.test_input_frame, text="Triaxial Input Data", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=(5, 0))
            self._add_test_type_dropdown(self.test_input_frame)
            self.triaxial_widgets = self._create_entries(
                self.test_input_frame,
                "",
                [INIT_PRESSURE_LABEL, MAX_STRAIN_LABEL, NUM_STEPS_LABEL, DURATION_LABEL],
                [FL2_UNIT_LABEL, PERCENTAGE_UNIT_LABEL, WITHOUT_UNIT_LABEL, SECONDS_UNIT_LABEL],
                {INIT_PRESSURE_LABEL: "100", MAX_STRAIN_LABEL: "20",
                 NUM_STEPS_LABEL: "100", DURATION_LABEL: "1.0"}
            )

        elif test_name == DIRECT_SHEAR:
            self._init_plot_canvas(num_plots=4)
            ttk.Label(self.test_input_frame, text="Direct Simple Shear Input Data", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=(5, 0))
            self._add_test_type_dropdown(self.test_input_frame)
            self.shear_widgets = self._create_entries(
                self.test_input_frame,
                "",
                [INIT_PRESSURE_LABEL, MAX_STRAIN_LABEL, NUM_STEPS_LABEL, DURATION_LABEL],
                [FL2_UNIT_LABEL, PERCENTAGE_UNIT_LABEL, WITHOUT_UNIT_LABEL, SECONDS_UNIT_LABEL],
                {INIT_PRESSURE_LABEL: "100", MAX_STRAIN_LABEL: "20",
                 NUM_STEPS_LABEL: "100", DURATION_LABEL: "1.0"}
            )

        log_message(f"{test_name} test selected.", "info")


    def _add_test_type_dropdown(self, parent):
        ttk.Label(parent, text="Type of Test:", font=("Arial", 10, "bold")).pack(anchor="w", padx=5, pady=(5, 2))

        self.test_type_var = tk.StringVar(value="Drained")
        self.test_type_menu = ttk.Combobox(
            parent,
            textvariable=self.test_type_var,
            values=["Drained"],
            state="readonly",
            width=12
        )
        self.test_type_menu.pack(anchor="w", padx=10, pady=(0, 10))

        self.test_type_menu.bind("<<ComboboxSelected>>")

    def _run_simulation(self):
        try:
            log_message("Starting calculation... Please wait...", "info")
            log_message("Validating input...", "info")
            self.root.update_idletasks()

            material_params = [e.get() for e in self.entry_widgets.values()]

            cohesion_phi_indices = None
            if not self.is_linear_elastic and self.mohr_checkbox.get():
                cohesion_phi_indices = (int(self.cohesion_var.get()), int(self.phi_var.get()))

            index = self.model_dict["model_name"].index(self.model_var.get()) + 1 if self.dll_path else None
            test_type = self.current_test.get()

            log_message("Calculating...", "info")
            self.root.update_idletasks()

            if test_type == TRIAXIAL:
                run_gui_builder(
                    test_type="triaxial",
                    dll_path=self.dll_path or "",
                    index=index,
                    material_parameters=[float(x) for x in material_params],
                    input_widgets=self.triaxial_widgets,
                    cohesion_phi_indices=cohesion_phi_indices,
                    axes=self.axes
                )

            elif test_type == DIRECT_SHEAR:
                run_gui_builder(
                    test_type="direct_shear",
                    dll_path=self.dll_path or "",
                    index=index,
                    material_parameters=[float(x) for x in material_params],
                    input_widgets=self.shear_widgets,
                    cohesion_phi_indices=cohesion_phi_indices,
                    axes=self.axes
                )

            self.canvas.draw()
            log_message(f"{test_type} test completed successfully.", "info")

        except Exception:
            log_message("An error occurred during simulation:", "error")
            log_message(traceback.format_exc(), "error")

        finally:
            self.root.after(0, self._enable_gui)
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

    def _disable_gui(self):
        self._set_widget_state(self.left_frame, "disabled")

    def _enable_gui(self):
        self._set_widget_state(self.left_frame, "normal")
        self.run_button.config(state="normal")

        if self.is_linear_elastic:
            self.mohr_checkbox_widget.configure(state="disabled")
            self.model_menu.configure(state="disabled")
        else:
            self.model_menu.configure(state="readonly")

    def _destroy_existing_plot_canvas(self):
        if hasattr(self, "plot_frame") and self.plot_frame.winfo_exists():
            for widget in self.plot_frame.winfo_children():
                widget.destroy()
        self.fig = None
        self.canvas = None
        self.axes = []
