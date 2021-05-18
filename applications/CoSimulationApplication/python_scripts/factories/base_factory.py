from importlib import import_module

def Create(settings, args, folder, module_name=None):
    """This function creates and returns a class from a module based on the use input
    First the requested module is searched in CoSimulation
    If it is not part of CoSimulation then it is attempted to be imported directly from PYTHONPATH
    """
    if module_name is None: # if the module name was not specified separately then take it from the settings
        module_name = settings["type"].GetString()
    full_module_path = ".".join([folder, module_name])

    try:
        imported_module = import_module(full_module_path)
    except ImportError:
        try:
            imported_module = import_module(module_name)
        except ImportError:
            raise ImportError('Module "{}" could neither be imported from CoSimulation nor from PYTHONPATH'.format(module_name))

    return imported_module.Create(settings, *args)
