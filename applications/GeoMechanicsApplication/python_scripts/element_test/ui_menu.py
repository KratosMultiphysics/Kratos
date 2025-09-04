import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk, scrolledtext, Menu
from platformdirs import user_data_dir
from pathlib import Path
from ui_builder import GeotechTestUI
from ui_udsm_parser import udsm_parser
from ui_labels import APP_TITLE, APP_VERSION, APP_NAME, APP_AUTHOR, SELECT_UDSM, LINEAR_ELASTIC, FONT_SEGOE_UI

import ctypes
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("deltares.ElementTestSuite.ui")

data_dir = Path(user_data_dir(APP_NAME, APP_AUTHOR))
data_dir.mkdir(parents=True, exist_ok=True)

LICENSE_FLAG_PATH = data_dir / "license_accepted.flag"


def show_license_agreement():
    license_text = """
    Pre-Release Software Licensing Agreement for testing of Pre-Release Software
    ----------------------------------------------------------------------------

    THE PARTIES
    1. STICHTING DELTARES, having its registered office and place of business in
    Delft (The Netherlands) at Boussinesqweg 1, listed in the Commercial Register under
    number 41146461, hereinafter called “Deltares “, in this present matter legally
    represented by C. van den Kieboom, Sales Officer Deltares Software Centre; and
    2. The end user, hereinafter called "Licensee"
    
    WHEREAS:
    • Deltares owns the intellectual property of the computer program being developed by
    Deltares, as described below, hereinafter referred to as “Pre-Release Software”;
    • The Pre-Release Software is under development by Deltares and is not fit for any use
    besides internal review, testing and evaluation;
    • Licensee wishes to acquire a non-exclusive and non-transferable license, without the
    right of sub-licensing, to use internally within the organization of Licensee, the
    PreRelease Software for review, testing and evaluation of functionality of the Pre-Release
    Software;
    • Deltares grants Licensee a Pre Release Software Licensing Agreement for review,
    testing and evaluation of Pre Release Software under the following conditions.
    
    AGREE AS FOLLOWS:
    
    Article 1 Definitions
    ---------------------
    Agreement This agreement, including the Appendices.
    Pre Release Software The computer program or the computer programs under development as
    described below under “Description Pre Release Software” and the Documentation and later
    versions of said computer program and Documentation.
    Documentation 	The manual or manuals and other documents that correspond to the
    Pre Release Software.
    
    Article 2 License
    -----------------
    Deltares grants Licensee until 1st of November 2022 starting at 1st of July 2022 and
    without a monetary charge, a non-exclusive, non-transferable, non sub-licensable license
    to use for the purpose and with the limitations specified in Article 3, the
    Pre Release Software, as described below under “Description Pre Release Software”
    for the purpose of research 
    
    Article 3 Use of Pre Release Software and the Documentation
    -----------------------------------------------------------
    1. Licensee shall only be authorised to use the Pre Release Software within its own
    organisation, for its own use and for review, testing and evaluation of functionality of
    the Pre-Release Software. Licensee shall not be permitted to make any other use of
    Pre Release Software, use the software for levee assignment or design, or to make
    available or grant access to Pre Release Software to any third party, subject to article 5
    paragraph 1 under b.
    2. Licensee shall not be authorised to modify and/or adjust Pre Release Software and/or
    to (otherwise) carry out alterations to it and/or to integrate Pre Release Software in
    other software, unless and only in so far as Licensee has obtained express written
    permission to that effect in advance from Deltares.
    3. Licensee shall not be authorised to (have others) copy Pre Release Software in any
    manner whatsoever or to (have others) multiply it (in any other way), except for backup
    purposes.
    
    Article 4 Intellectual Property Rights, Ownership
    -------------------------------------------------
    Deltares owns the copyright to Pre Release Software. With the Agreement Deltares only
    grants Licensee the license rights in connection with Pre Release Software as described
    in Article 2. Licensee accepts that with this license granted no transfer of ownership
    whatsoever, including the transfer of the intellectual property rights, is made
    to Licensee.
    
    Article 5 Confidentiality
    -------------------------
    1. Licensee shall keep confidential Pre Release Software which Licensee has obtained
    and/or obtains, in any manner, from Deltares under or in connection with the
    Agreement.
    This obligation shall at any rate include:
    a. Treating of Pre Release Software confidentially;
    b. releasing Pre Release Software solely to those employees of Licensees or a third party
    acting on behalf of Licensee under the conditions of this Agreement who require access to
    Pre Release Software, whereby Licensee will oblige these employees and third parties to
    the same confidentiality as Licensee;
    c. the non-disclosure of information and/or data related to the Agreement to third parties
    and/or refraining from making such information and/or data public in any other way without
    the prior express and written consent of Deltares, to be obtained for each separate event.
    d. using information and/or data obtained solely for the purposes for which they were
    obtained.
    2. Licensee's obligation of confidentiality referred to in Article 5.1 shall not apply to
    information and/or data that were already at Licensee's free disposal, or were part of
    the public domain, or were already included in generally accessible literature at the
    time when they were obtained by Licensee, or that were obtained by Licensee from a
    third party or third parties who was or were free to disclose the relevant information
    and/or data and who had not obtained the information and/or data from Deltares.
    3. The provisions in this article and article 8 shall remain in full force after
    termination of
    the Agreement, as set forth in Article 7, as well.
    
    Article 6 No guarantee, no warrantee
    ------------------------------------
    The Pre -Release Software is under development by Deltares. Licensee acknowledges and
    agrees that the Pre Release Software may contain bugs, defects, errors and may not be
    expected to function fully upon installation nor that results obtained with the
    Pre Release Software are correct or of proper quality. Licensee also acknowledges
    that Deltares is under no obligation to correct any bugs, defects, or errors in
    the Pre Release Software or to otherwise support or maintain the Pre Release Software. 
    
    Article 7 Duration, Termination
    -------------------------------
    1. The Agreement concluded for a period until 1st of November 2022 starting at
    1st of July 2022, subject to termination in accordance with the provisions of
    article 7.2.
    2. Parties are entitled without cause being required, to terminate the Agreement in
    writing with immediate effect, without judicial intervention being required.
    3. In the event of termination of the Agreement, Licensee shall immediately return to
    Deltares all copies of Pre Release Software.
    4. In addition, in the event of termination of the Agreement, Licensee shall furthermore
    immediately cease using (additional copies of) Pre Release Software and delete Pre
    Release Software from (all) their computer(s).
    5. "In writing" as mentioned in this article shall also be a fax message but not an e-mail
    message.
    
    Article 8 Liability
    -------------------
    1. Licensee agrees that Deltares (including its personnel and non-employees who (have)
    undertake(n) activities for Deltares) shall not be responsible to Licensee for any
    loss-of, direct, indirect, incidental, special or consequential damages arising out
    of the license agreement or the use of Pre Release Software, to the extent permitted by
    Netherlands law.
    2. Licensee shall indemnify, hold harmless and defend Deltares against any action
    brought by a third party against Deltares to the extent that such a claim is connected to
    the use of Pre Release Software by Licensee and/or third parties at whose disposal the
    Licensee has placed Pre Release Software in accordance with this Agreement and/or
    these results or to whom he has otherwise made them known the results, including use
    of the results of use by Licensee and/or third parties and the installation of the Pre
    Release Software by Licensee.
    
    Article 9 Other provisions
    --------------------------
    1. Changes in and/or deviations to the Agreement are valid only if they are explicitly
    agreed between the parties in writing.
    2. The parties are not allowed to assign any rights and/or obligations under the
    Agreement, entirely or in part, to third parties without the prior written consent of the
    other party.
    3. Any disputes arising from the Agreement or from agreements arising there from, shall
    be submitted solely to the competent court of The Hague.
    4. This Agreement and all the agreements arising there from are governed exclusively by
    Netherlands law 
    """

    license_window = tk.Toplevel()
    license_window.title("Pre-Release License Agreement")
    license_window.geometry("800x600")
    license_window.grab_set()
    license_window.protocol("WM_DELETE_WINDOW", lambda: os._exit(0))

    tk.Label(license_window, text="Please review and accept the agreement to continue.",
             font=("Arial", 12, "bold"), pady=10).pack()

    text_area = scrolledtext.ScrolledText(license_window, wrap="word", font=("Courier", 10))
    text_area.insert("1.0", license_text)
    text_area.config(state="disabled")
    text_area.pack(expand=True, fill="both", padx=10, pady=10)

    button_frame = tk.Frame(license_window)
    button_frame.pack(pady=10)

    def accept():
        try:
            with open(LICENSE_FLAG_PATH, "w") as f:
                f.write("ACCEPTED")
        except Exception as e:
            messagebox.showerror("Error", f"Could not save license acceptance: {e}")
            os._exit(1)
        license_window.destroy()

    def decline():
        messagebox.showinfo("Exit", "You must accept the license agreement to use this software.")
        os._exit(0)

    tk.Button(button_frame, text="Accept", width=15, command=accept).pack(side="left", padx=10)
    tk.Button(button_frame, text="Decline", width=15, command=decline).pack(side="right", padx=10)

def show_about_window():
    about_win = tk.Toplevel()
    about_win.title("About")
    about_win.geometry("500x400")
    about_win.resizable(False, False)
    about_win.grab_set()

    tk.Label(about_win, text=APP_TITLE, font=(FONT_SEGOE_UI, 14, "bold")).pack(pady=(20, 5))
    tk.Label(about_win, text=APP_VERSION, font=(FONT_SEGOE_UI, 12)).pack(pady=(0, 5))
    tk.Label(about_win, text="Powered by:", font=(FONT_SEGOE_UI, 12)).pack(pady=(0, 5))

    image_frame = tk.Frame(about_win)
    image_frame.pack(pady=10)

    try:
        path1 = os.path.join(os.path.dirname(__file__), "assets", "kratos.png")
        path2 = os.path.join(os.path.dirname(__file__), "assets", "deltares.png")

        photo1 = tk.PhotoImage(file=path1)
        photo2 = tk.PhotoImage(file=path2)

        label1 = tk.Label(image_frame, image=photo1)
        label1.image = photo1
        label1.pack(pady=2)

        label2 = tk.Label(image_frame, image=photo2)
        label2.image = photo2
        label2.pack(pady=15)

    except Exception:
        tk.Label(about_win, text="[One or both images could not be loaded]", fg="red").pack()

    tk.Label(about_win, text="Contact: kratos@deltares.nl", font=(FONT_SEGOE_UI, 12)).pack(pady=(0, 2))
    tk.Button(about_win, text="Close", command=about_win.destroy).pack(pady=10)


def create_menu():
    root = tk.Tk()

    menubar = Menu(root)
    root.config(menu=menubar)

    file_menu = Menu(menubar, tearoff=0)
    file_menu.add_command(label="Exit", command=lambda: root.quit())
    menubar.add_cascade(label="File", menu=file_menu)

    about_menu = Menu(menubar, tearoff=0)
    about_menu.add_command(label="License", command=show_license_agreement)
    about_menu.add_command(label="About", command=show_about_window)
    menubar.add_cascade(label="Help", menu=about_menu)

    if not os.path.exists(LICENSE_FLAG_PATH):
        show_license_agreement()

    try:
        icon_path = os.path.join(os.path.dirname(__file__), "assets", "icon.ico")
        root.iconbitmap(default=icon_path)
    except Exception as e:
        print(f"Could not set icon: {e}")

    root.title(f"{APP_TITLE} - {APP_VERSION}")
    root.state('zoomed')
    root.resizable(True, True)

    top_frame = ttk.Frame(root, padding="10")
    top_frame.pack(side="top", fill="x")

    main_frame = ttk.Frame(root)
    main_frame.pack(side="top", fill="both", expand=True)

    def load_dll():
        dll_path = filedialog.askopenfilename(title=SELECT_UDSM, filetypes=[("DLL files", "*.dll")])
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
        if choice == SELECT_UDSM:
            load_dll()
        elif choice == LINEAR_ELASTIC:
            load_linear_elastic()

    model_source_var = tk.StringVar(value="Select Model Source")
    model_source_menu = ttk.Combobox(
        top_frame,
        textvariable=model_source_var,
        values=[SELECT_UDSM, LINEAR_ELASTIC],
        state="readonly"
    )
    model_source_menu.bind("<<ComboboxSelected>>", handle_model_source_selection)
    model_source_menu.pack(side="left", padx=5)

    def on_close():
        root.quit()
        root.destroy()
        os._exit(0)

    root.protocol("WM_DELETE_WINDOW", on_close)
    root.mainloop()

if __name__ == "__main__":
    create_menu()
