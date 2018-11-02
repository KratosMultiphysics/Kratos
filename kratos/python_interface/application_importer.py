from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os.path
import sys
import KratosMultiphysics
from KratosMultiphysics import Logger


def ImportApplication(application, application_name, application_folder, caller, mod_path=None):
    Globals = KratosMultiphysics.KratosGlobals
    Kernel = Globals.Kernel
    main_caller = Globals.AuthorizedCaller
    applications_root = Globals.ApplicationsRoot
    # caller and main_caller are generated from the output of inspect.stack().
    # In particular position [1] is the name of the file containing the call.
    # Note that position [0] (a frame instance) could also be used for the check,
    # but can return false if both calls are made from the python interpreter
    if main_caller[1] != caller[1]:
        msg = "\n***\n*    Python file " + str(caller[1]) + "\n*    requires " + str(application_name) + "\n*    Please import it from your main Python script, " +str(main_caller[1]+"\n*    If your main Python script is already importing it, you may need to re-order the imports"+'\n***')
        # print caller
        # print main_caller
        raise RuntimeError(msg)
    elif application_name not in Globals.RequestedApplications:  # This check is possibly redundant, as Python won't import the same module twice
        Logger.PrintInfo("", "Importing    " + application_name)
        # Add application to dictionary of registered applications
        Globals.RequestedApplications[application_name] = application
        # Add python scrips folder to path
        application_path = os.path.join(applications_root, application_folder)
        python_path = os.path.join(application_path, 'python_scripts')
        sys.path.append(python_path)
        # Add constitutive laws python scrips folder to path
        constitutive_laws_path = os.path.join(python_path, 'constitutive_laws')
        sys.path.append(constitutive_laws_path)

        if mod_path is not None: # optional for backwards compatibility
            mod_path.append(python_path)

        # Add application to kernel
        Kernel.ImportApplication(application)
        # Dynamic renumbering of variables to ensure consistency
        Kernel.Initialize()
        for iterName, iterApplication in list(Globals.RequestedApplications.items()):
            # print("Initializing",iterName)
            Kernel.InitializeApplication(iterApplication)
