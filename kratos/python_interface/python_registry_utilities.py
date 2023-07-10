import KratosMultiphysics

def RegisterAll(PythonModuleName, PythonRegistryListModule):
    RegisterModelersList(PythonModuleName, PythonRegistryListModule)
    RegisterOperationsList(PythonModuleName, PythonRegistryListModule)
    RegisterProcessesList(PythonModuleName, PythonRegistryListModule)
    RegisterStagesList(PythonModuleName, PythonRegistryListModule)
    RegisterOrchestratorsList(PythonModuleName, PythonRegistryListModule)

def RegisterModelersList(PythonModuleName, PythonRegistryListModule):
    __CheckRegistryListIsInModule(PythonRegistryListModule, "modelers")
    __RegisterItemList(PythonModuleName, PythonRegistryListModule.python_modelers_to_be_registered, "Modelers")

def RegisterOperationsList(PythonModuleName, PythonRegistryListModule):
    __CheckRegistryListIsInModule(PythonRegistryListModule, "operations")
    __RegisterItemList(PythonModuleName, PythonRegistryListModule.python_operations_to_be_registered, "Operations")

def RegisterProcessesList(PythonModuleName, PythonRegistryListModule):
    __CheckRegistryListIsInModule(PythonRegistryListModule, "processes")
    __RegisterItemList(PythonModuleName, PythonRegistryListModule.python_processes_to_be_registered, "Processes")

def RegisterStagesList(PythonModuleName, PythonRegistryListModule):
    __CheckRegistryListIsInModule(PythonRegistryListModule, "stages")
    __RegisterItemList(PythonModuleName, PythonRegistryListModule.python_stages_to_be_registered, "Stages")

def RegisterOrchestratorsList(PythonModuleName, PythonRegistryListModule):
    __CheckRegistryListIsInModule(PythonRegistryListModule, "orchestrators")
    __RegisterItemList(PythonModuleName, PythonRegistryListModule.python_orchestrators_to_be_registered, "Orchestrators")

def __CheckRegistryListIsInModule(PythonRegistryListModule, ListKeyword):
    list_variable_name = f"python_{ListKeyword}_to_be_registered"
    if not list_variable_name in dir(PythonRegistryListModule):
        err_msg = f"Trying to register '{ListKeyword}' in Python registry but the corresponding list is not provided in '{PythonRegistryListModule.__name__}'."
        raise Exception(err_msg)

def __RegisterItemList(PythonModuleName, ItemsList, ItemKeyword):
    for item in ItemsList:
        class_name = item.split('.')[-1]
        module_name = f'{PythonModuleName}.{".".join(item.split(".")[:-1])}'

        registry_entry_point = f"{ItemKeyword}.{PythonModuleName}.{class_name}"
        registry_entry_data = {
            "ClassName" : f"{class_name}",
            "ModuleName" : f"{module_name}"
        }
        KratosMultiphysics.Registry.AddItem(registry_entry_point, registry_entry_data)