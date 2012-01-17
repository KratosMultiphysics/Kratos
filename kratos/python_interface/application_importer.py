import os.path
import sys
import KratosMultiphysics

def ImportApplication(application,application_name,application_folder,caller):
  Globals = KratosMultiphysics.KratosGlobals
  Kernel = Globals.Kernel
  main_caller = Globals.AuthorizedCaller
  applications_root = Globals.ApplicationsRoot
  # caller and main_caller are generated from the output of inspect.stack(). 
  # In particular position [1] is the name of the file containing the call.
  # Note that position [0] (a frame instance) could also be used for the check,
  # but can return false if both calls are made from the python interpreter
  if main_caller[1] != caller[1]:
    msg = "Python file "+str(caller[1])+" requires "+str(application_name)+"\n    Please import it from your main Python script, "+str(main_caller[1])
    #print caller
    #print main_caller
    raise RuntimeError(msg)
  elif Globals.RequestedApplications.has_key(application_name) == False: # This check is possibly redundant, as Python won't import the same module twice
    print "Importing "+application_name
    # Add application to dictionary of registered applications
    Globals.RequestedApplications[application_name] = application
    # Add python scrips folder to path
    application_path = os.path.join(applications_root,application_folder)
    python_path = os.path.join(application_path,'python_scripts')
    sys.path.append(python_path)
    # Add application to kernel
    Kernel.AddApplication(application)
    # Dynamic renumbering of variables to ensure consistency
    Kernel.Initialize()
    for iterName,iterApplication in Globals.RequestedApplications.items():
      print "Initializing ",iterName
      Kernel.InitializeApplication(iterApplication)
