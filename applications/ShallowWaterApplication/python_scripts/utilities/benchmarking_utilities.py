
def GetProcessParameters(list_of_processes, name):
    for processes in list_of_processes.values():
        for process in processes:
            if process['python_module'].GetString() == name:
                return process['Parameters']

def GetModelerParameters(list_of_modelers, name):
    for modeler in list_of_modelers:
        if modeler['modeler_name'].GetString() == name:
            return modeler['Parameters']

def KeepOnlyThisProcess(list_of_processes, module):
    for list_name, proc_list in list_of_processes.items():
        for proc in proc_list:
            if proc['python_module'].GetString() == module:
                keep_proc = proc.Clone()
                keep_list = list_name
        list_of_processes.RemoveValue(list_name)
    list_of_processes.AddEmptyArray(keep_list)
    list_of_processes[keep_list].Append(keep_proc)
