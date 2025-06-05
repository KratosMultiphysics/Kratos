import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from ui_builder import GeotechTestUI
from ui_udsm_parser import udsm_parser


SELECT_DLL = "Select DLL File"
LINEAR_ELASTIC = "Linear Elastic Model"

def create_menu():
    root = tk.Tk()
    root.title("Triaxial Test")
    root.state('zoomed')
    root.resizable(True, True)

    top_frame = ttk.Frame(root, padding="10")
    top_frame.pack(side="top", fill="x")

    main_frame = ttk.Frame(root)
    main_frame.pack(side="top", fill="both", expand=True)

    def load_dll():
        dll_path = filedialog.askopenfilename(title=SELECT_DLL, filetypes=[("DLL files", "*.dll")])
        if not dll_path:
            messagebox.showerror("Error", "No DLL file selected.")
            return

        try:
            model_dict = udsm_parser(dll_path)
        except Exception as e:
            messagebox.showerror("DLL Error", f"Failed to parse DLL: {e}")
            return

        for widget in main_frame.winfo_children():
            widget.destroy()

        GeotechTestUI(root, main_frame, test_name="Triaxial", dll_path=dll_path, model_dict=model_dict,
                      external_widgets=[model_source_menu])

    def load_linear_elastic():
        model_dict = {
            "model_name": [LINEAR_ELASTIC],
            "num_params": [2],
            "param_names": [["YOUNG MODULUS", "POISSON RATIO"]],
            "param_units": [["kN/m²", "–"]]
        }

        for widget in main_frame.winfo_children():
            widget.destroy()

        GeotechTestUI(root, main_frame, test_name="Triaxial", dll_path=None, model_dict=model_dict,
                      external_widgets=[model_source_menu])

    def handle_model_source_selection(event):
        choice = model_source_var.get()
        if choice == SELECT_DLL:
            load_dll()
        elif choice == LINEAR_ELASTIC:
            load_linear_elastic()

    model_source_var = tk.StringVar(value="Select Model Source")
    model_source_menu = ttk.Combobox(
        top_frame,
        textvariable=model_source_var,
        values=[SELECT_DLL, LINEAR_ELASTIC],
        state="readonly"
    )
    model_source_menu.bind("<<ComboboxSelected>>", handle_model_source_selection)
    model_source_menu.pack(side="left", padx=5)

    root.mainloop()

if __name__ == "__main__":
    create_menu()