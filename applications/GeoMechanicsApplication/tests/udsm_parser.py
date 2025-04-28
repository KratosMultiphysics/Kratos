import ctypes
import pefile
import string
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import plotly.graph_objects as go
from set_triaxial_test import lab_test
from matplotlib.gridspec import GridSpec
import numpy as np

import re

def clean_c_buffer(char_buffer):
    try:
        raw = char_buffer.raw.split(b'\x00')[0]
        decoded = raw.decode("utf-8", errors="ignore")
        # Remove leading numbers and spaces (like "0    ")
        decoded = re.sub(r"^\s*\d+\s+", "", decoded)
        # Strip non-display chars
        decoded_clean = re.sub(r"[^\w\s\-\(\)\[\]\.,':=+/]", "", decoded)
        return decoded_clean.strip()
    except Exception:
        return "<?>"  # fallback

def find_symbol_in_dll(dll_path, dll_lib, symbol_name):
    try:
        pe = pefile.PE(dll_path)
        symbol_name_lower = symbol_name.lower()

        # Iterate over all exported symbols
        for exp in pe.DIRECTORY_ENTRY_EXPORT.symbols:
            name = exp.name
            if name and name.decode().lower() == symbol_name_lower:
                return getattr(dll_lib, name.decode())  # Return the correct symbol name

        return None
    except Exception as e:
        print(f"Error reading DLL: {e}")
        return None

def get_model_count(getmodelcount):
    try:
        getmodelcount.argtypes = (ctypes.POINTER(ctypes.c_int),)
        getmodelcount.restype = None
        result = ctypes.c_int()
        getmodelcount(ctypes.byref(result))
    except:
        raise Exception("Function GetModelCount not found in DLL")
    return result.value

def get_model_name(getmodelname, model_no):
    BUFFER_SIZE = 256
    char_buffer = ctypes.create_string_buffer(BUFFER_SIZE)
    getmodelname.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_char_p, ctypes.c_long)
    getmodelname.restype = None
    getmodelname(ctypes.byref(ctypes.c_int(model_no)), char_buffer, ctypes.c_long(256))
    return clean_c_buffer(char_buffer)

def get_param_count(getparamcount, model_no):
    getparamcount.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
    getparamcount.restype = None
    result = ctypes.c_int()
    getparamcount(ctypes.byref(ctypes.c_int(model_no)), ctypes.byref(result))
    return result.value

def get_param_name(getparamname, model_no, param_no):
    BUFFER_SIZE = 256
    char_buffer = ctypes.create_string_buffer(BUFFER_SIZE)
    getparamname.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                             ctypes.c_char_p, ctypes.c_long)
    getparamname.restype = None
    getparamname(ctypes.byref(ctypes.c_int(model_no)), ctypes.byref(ctypes.c_int(param_no)),
                 char_buffer, ctypes.c_long(256))
    return clean_c_buffer(char_buffer)

def get_state_var_count(getstatevarcount, model_no):
    getstatevarcount.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
    getstatevarcount.restype = None
    result = ctypes.c_int()
    getstatevarcount(ctypes.byref(ctypes.c_int(model_no)), ctypes.byref(result))
    return result.value

def get_state_var_name(getstatevarname, model_no, state_var_no):
    BUFFER_SIZE = 256
    char_buffer = ctypes.create_string_buffer(BUFFER_SIZE)
    getstatevarname.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                ctypes.c_char_p, ctypes.c_long)
    getstatevarname.restype = None
    getstatevarname(ctypes.byref(ctypes.c_int(model_no)), ctypes.byref(ctypes.c_int(state_var_no)),
                    char_buffer, ctypes.c_long(256))
    return clean_c_buffer(char_buffer)

def get_material_files(material_directory):
    # Get the list of files in the material directory
    material_files = os.listdir(material_directory)
    return material_files
def udsm_parser(dll_path):
    """
    This function parses a user defined soil model (UDSM) from a DLL file.
    It extracts the model name, number of parameters, and the parameter names.
    """
    import ctypes

    # Load the DLL
    dll_lib = ctypes.CDLL(dll_path, winmode=0)

    # Get the function pointers
    getmodelcount = find_symbol_in_dll(dll_path, dll_lib, "getmodelcount")
    getmodelname = find_symbol_in_dll(dll_path, dll_lib, "getmodelname")
    getparamcount = find_symbol_in_dll(dll_path, dll_lib, "getparamcount")
    getparamname = find_symbol_in_dll(dll_path, dll_lib, "getparamname")
    #getstatevarcount = find_symbol_in_dll(dll_path, dll_lib, "getstatevarcount")
    #getstatevarname = find_symbol_in_dll(dll_path, dll_lib, "getstatevarname")

    # Close the DLL
    model_count = get_model_count(getmodelcount)

    # create library object with the structure that contains the parameters from each of the models in the file
    model_dict = {"model_name": [], "num_params": [], "param_names": []}
    for model_no in range(1,model_count+1):
        model_dict["model_name"].append(get_model_name(getmodelname, model_no))
        num_params = get_param_count(getparamcount, model_no)
        model_dict["num_params"].append(num_params)
        param_names = []
        for param in range(1,num_params+1):
            param_names.append(get_param_name(getparamname, model_no, param))
        model_dict["param_names"].append(param_names)

    return model_dict

def create_menu():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    umat_entry_parameters = []

    # Prompt user to select a DLL file
    dll_path = filedialog.askopenfilename(title="Select DLL File", filetypes=[("DLL files", "*.dll")])
    if not dll_path:
        messagebox.showerror("Error", "No DLL file selected.")
        return

    # Parse the UDSM from the DLL file
    model_dict = udsm_parser(dll_path)

    # Create a new window for the menu
    menu_window = tk.Toplevel()
    menu_window.title("Triaxial Test")
    menu_window.state('zoomed')  # Maximize the window
    menu_window.resizable(True, True)

    # Global variables for the plot canvas and axes
    global fig, axes, canvas
    fig = None
    axes = None
    canvas = None


# Create a frame for the dropdown menu
    dropdown_frame = ttk.Frame(menu_window, padding="10")
    dropdown_frame.pack(side="left", fill="y", padx=10, pady=10)

    ttk.Label(dropdown_frame, text="Select a Model:", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
    model_var = tk.StringVar(menu_window)
    model_var.set(model_dict["model_name"][0])  # Set default value
    model_menu = ttk.Combobox(dropdown_frame, textvariable=model_var, values=model_dict["model_name"], state="readonly")
    model_menu.pack(side="top", fill="x", expand=True, padx=5)

    # Create a frame for the parameters
    param_frame = ttk.Frame(dropdown_frame, padding="10")
    param_frame.pack(fill="both", expand=True, pady=10)

    button_frame = ttk.Frame(menu_window, padding="10")
    button_frame.pack(side="bottom", fill="x", padx=10, pady=10)

    def update_parameters(*args):
        # Clear the parameter frame
        for widget in param_frame.winfo_children():
            widget.destroy()

        # Dictionary to store entry widgets
        entry_widgets = {}
        triaxial_entry_widgets = {}

        # Define default values for parameters
        default_values = {
            "1. E": "10000",  # Example default value
            "2. n_ur": "0.3",
            "3. c'": "0.0",
            "4. f_peak": "30.0",
            "5. y_peak": "0.0",
            "6. s_t, cut-off": "0.0",
            "7. yield function (MC=1 DP=2 MNC=3 MN=4)": "1",
            "8. n_un (UMAT)": "0.3"
        }

        # Get the selected model and its parameters
        selected_model = model_var.get()
        index = model_dict["model_name"].index(selected_model)
        params = model_dict["param_names"][index]

        # Add a label for "Triaxial Input Data"
        parameter_frame = ttk.Frame(param_frame, padding="10")
        parameter_frame.pack(fill="x", pady=10)
        ttk.Label(parameter_frame, text="Soil Input Parameters", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)
        # Add text input fields for each parameter
        for param in params:
            param_row = ttk.Frame(param_frame)
            param_row.pack(fill="x", padx=10, pady=2)

            ttk.Label(param_row, text=param, font=("Arial", 10)).pack(side="left", padx=5)
            entry = ttk.Entry(param_row)
            entry.pack(side="left", fill="x", expand=True)

            # Set default value if available
            default_value = default_values.get(param, "N/A")  # Use "N/A" if no default value is defined
            entry.insert(0, default_value)

            entry_widgets[param] = entry  # Store the entry widget

        # Add a label for "Triaxial Input Data"
        triaxial_frame = ttk.Frame(param_frame, padding="10")
        triaxial_frame.pack(fill="x", pady=10)

        ttk.Label(triaxial_frame, text="Triaxial Input Data", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=5)

        triaxial_inputs = [
            ("Initial effective cell pressure |σ'\u2093\u2093|", "kN/m²", "100"),
            ("Maximum Strain |ε\u1d67\u1d67|", "%", "10"),
            ("Number of steps", "", "100")
        ]

        for label_text, unit, default_value in triaxial_inputs:
            input_row = ttk.Frame(triaxial_frame)
            input_row.pack(fill="x", padx=10, pady=2)

            ttk.Label(input_row, text=label_text, font=("Arial", 10)).pack(side="left", padx=5)
            entry = ttk.Entry(input_row)
            entry.pack(side="left", fill="x", expand=True, padx=5)

            # Set default value
            entry.insert(0, default_value)

            ttk.Label(input_row, text=unit, font=("Arial", 10)).pack(side="left", padx=5)
            triaxial_entry_widgets[label_text] = entry  # Store the entry widget

        # Save the entry widgets globally for access in run_calculation
        global parameter_entries, input_entries, model_index
        parameter_entries = entry_widgets
        input_entries = triaxial_entry_widgets
        model_index = index + 1
    def run_calculation():
        # Retrieve values from the entry widgets
        parameters = [entry.get() for key, entry in parameter_entries.items()]
        try:
            initial_effective_stress = float(input_entries["Initial effective cell pressure |σ'\u2093\u2093|"].get())
            maximum_strain = float(input_entries["Maximum Strain |ε\u1d67\u1d67|"].get())
            time_step = float(input_entries["Number of steps"].get())
        except ValueError:
            messagebox.showerror("Error", "Invalid input for 'Triaxial Input Data'. Please enter numeric values.")
            return

        figs = lab_test(dll_path, model_index, parameters, time_step, maximum_strain, initial_effective_stress)
        for i, ax in enumerate(axes):
            ax.clear()

            if i < len(figs):  # Make sure there's a plotly fig for this subplot
                fig = figs[i]

                # Extract x and y labels from the Plotly figure
                x_label = fig.layout.xaxis.title.text if fig.layout.xaxis.title else None
                y_label = fig.layout.yaxis.title.text if fig.layout.yaxis.title else None

                # Extract all traces and plot them
                for trace in fig.data:
                    if isinstance(trace, go.Scatter):
                        x = trace.x
                        y = trace.y
                        label = trace.name if hasattr(trace, 'name') else None
                        ax.plot(x, y, label=label)

                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)

                ax.invert_xaxis()
                if i in [1, 2, 4]:  # Flip y-axis for fig 1 and 3
                    ax.invert_yaxis()
                    ax.set_ylim(0, np.min(y)*1.2)  # Ensure x-axis starts from 0
                else:
                    ax.set_ylim(0, np.max(y)*1.2)

                ax.set_xlim(0, np.min(x)*1.2)  # Ensure x-axis starts from 0


            titles = ["Delta Sigma", "Volumetric Strain", "Sigma Plot", "p-q Plot", "Mohr-Coulomb Circle"]
            ax.set_title(titles[i] if i < len(titles) else f"Plot {i + 1}")
            ax.legend()

        canvas.draw()
    run_button = ttk.Button(button_frame, text="Run Calculation", command=run_calculation)
    run_button.pack(pady=5)

    # Bind the dropdown menu to update parameters on selection change
    model_var.trace("w", update_parameters)

    # Initialize the parameters for the default model
    update_parameters()

    # Create a frame for the plots
    plot_frame = ttk.Frame(menu_window, padding="5")
    plot_frame.pack(side="right", fill="both", expand=True, padx=5, pady=5)

    # Desired size per subplot
    width_per_ax = 4
    height_per_ax = 3

    # Grid layout
    nrows, ncols = 3, 2
    figsize = (ncols * width_per_ax, nrows * height_per_ax)

    # Define grid layout with spacing
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(nrows, ncols, figure=fig, wspace=0.4, hspace=0.6)  # Adjust wspace and hspace for spacing

    # Create subplots with individual aspect ratios
    axes = []
    for i in range(nrows * ncols):
        if i < 5:  # Only create 5 subplots
            ax = fig.add_subplot(gs[i])
            ax.plot([0, 1, 2, 3], [0, 0, 0, 0])  # Placeholder data
            # ax.set_title(f"Plot {i + 1}")
            # ax.set_aspect('equal') #if i % 2 == 0 else 'equal')  # Example: alternate aspect ratios
            axes.append(ax)

    # Hide the 6th subplot
    if len(axes) == 6:
        fig.delaxes(axes[5])  # or axes[5].set_visible(False)

    # Embed the matplotlib figure in the tkinter frame
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack(fill="both", expand=True)

    menu_window.mainloop()


if __name__ == "__main__":

    create_menu()
    # Example usage
    #dll_path = r"c:\Users\nuttall\OneDrive - Stichting Deltares\Desktop\Bugs\ColumnWetDry\udsm.dll"
    #model_name, num_params, param_names = udsm_parser(dll_path)

