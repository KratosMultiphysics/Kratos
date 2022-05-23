
def GetProcessParameters(list_of_processes, name):
    for processes in list_of_processes.values():
        for process in processes:
            if process['python_module'].GetString() == name:
                return process['Parameters']

def GetModelerParameters(list_of_modelers, name):
    for modeler in list_of_modelers:
        if modeler['modeler_name'].GetString() == name:
            return modeler['Parameters']
