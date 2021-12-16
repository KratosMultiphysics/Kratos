import KratosMultiphysics as KM

# Other imports
from importlib import import_module

class NeuralNetworkProcessFactory(object):
    """
    This class is the factory for the neural network processes.
    """
    def __init__(self):
        pass

    def ConstructListOfProcesses( self, process_list, model = None ):
        constructed_processes = []
        if process_list.IsArray():
            processes = process_list
            size = processes.size()
        else:
            processes = [process_list]
            size = len(processes)
        for i in range(0,size):
            item = processes[i]
            if not item.Has("python_module"):
                KM.Logger.PrintWarning("Your list of processes: ", processes)
                raise NameError('"python_module" must be defined in your parameters. Check all your processes')

            # python-script that contains the process
            python_module_name = item["python_module"].GetString()

            if item.Has("kratos_module"): # for Kratos-processes
                """Location of the process in Kratos; e.g.:
                - KratosMultiphysics
                - KratosMultiphysics.FluidDynamicsApplication
                """

                kratos_module_name = item["kratos_module"].GetString()
                if not kratos_module_name.startswith("KratosMultiphysics"):
                    kratos_module_name = "KratosMultiphysics." + kratos_module_name

                full_module_name = kratos_module_name + "." + python_module_name
                python_module = import_module(full_module_name)
                if kratos_module_name == "KratosMultiphysics.NeuralNetworkApplication":
                    p = python_module.Factory(item)
                else:
                    p = python_module.Factory(item, model)
                constructed_processes.append( p )

            else: # for user-defined processes
                python_module = import_module(python_module_name)
                p = python_module.Factory(item)
                constructed_processes.append( p )

        return constructed_processes