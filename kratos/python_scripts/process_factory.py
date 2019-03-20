import KratosMultiphysics

################# Please do not change the following class:
class KratosProcessFactory(object):
    def __init__(self, Model):
        self.Model = Model #Model is a place

    def ConstructListOfProcesses( self, process_list ):
        constructed_processes = []
        for i in range(0,process_list.size()):
            item = process_list[i]
            if not item.Has("python_module"):
                KratosMultiphysics.Logger.PrintWarning("Your list of processes: ", process_list)
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

                try:
                    full_module_name = kratos_module_name + "." + python_module_name
                    python_module = __import__(full_module_name,fromlist=[python_module_name])
                except ImportError: # old import mechanism for backwards-compatibility
                    __import__(kratos_module_name)
                    python_module = __import__(python_module_name)

                p = python_module.Factory(item, self.Model)
                constructed_processes.append( p )

            else: # for user-defined processes
                python_module = __import__(python_module_name)
                p = python_module.Factory(item, self.Model)
                constructed_processes.append( p )

        return constructed_processes


########## here we generate the common kratos processes --- IMPLEMENTED IN C++ ###################
def Factory(settings, Model):
    if(settings["process_name"].GetString() == "ApplyConstantScalarValueProcess"):
        model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
        return KratosMultiphysics.ApplyConstantScalarValueProcess(model_part, settings["Parameters"])

    elif(settings["process_name"].GetString() == "ApplyConstantVectorValueProcess"):
        model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
        return KratosMultiphysics.ApplyConstantVectorValueProcess(model_part, settings["Parameters"])

    raise Exception("Process name not found ",)
