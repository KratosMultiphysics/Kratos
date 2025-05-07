import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from ui_builder import build_ui_from_model
from ui_udsm_parser import udsm_parser


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
        dll_path = filedialog.askopenfilename(title="Select DLL File", filetypes=[("DLL files", "*.dll")])
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

        build_ui_from_model(root, main_frame, dll_path, model_dict)

    select_dll_button = ttk.Button(top_frame, text="Select DLL File", command=load_dll)
    select_dll_button.pack(side="left")

    root.mainloop()

if __name__ == "__main__":
    create_menu()