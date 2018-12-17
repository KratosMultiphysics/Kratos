from KratosMultiphysics import *

#################33please do not change the following class
class KratosProcessFactory(object):
    def __init__(self, Model):
        self.Model = Model #Model is a place

    def ConstructListOfProcesses( self, process_list ):
        constructed_processes = []
        for i in range(0,process_list.size()):
            item = process_list[i]
            # The kratos_module is the application where the script must be loaded. ex. KratosMultiphysics.StructuralMechanicsApplication
            if(item.Has("kratos_module")):
                kratos_module = __import__(item["kratos_module"].GetString())
            # The python_module is the actual scrpt that must be launch
            if(item.Has("python_module")):
                python_module = __import__(item["python_module"].GetString())
                p = python_module.Factory(item, self.Model)
                constructed_processes.append( p )
            else: # Otherwise an error is thrown
                Logger.PrintWarning("Your list of processes: ", process_list)
                raise NameError("python_module must be defined in your parameters. Check all your processes")


        return constructed_processes


########## here we generate the common kratos processes --- IMPLEMENTED IN C++ ###################
def Factory(settings, Model):
    if(settings["process_name"].GetString() == "ApplyConstantScalarValueProcess"):
        model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
        return ApplyConstantScalarValueProcess(model_part, settings["Parameters"])

    elif(settings["process_name"].GetString() == "ApplyConstantVectorValueProcess"):
        model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
        return ApplyConstantVectorValueProcess(model_part, settings["Parameters"])

    raise Exception("process name not found ",)
