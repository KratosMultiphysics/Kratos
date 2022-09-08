# Importing the import module in order to be able to select a base class
from importlib import import_module

# Importing the base factory
from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateSolverWrapper(settings, models, solver_name):
    """This function creates and returns the Wrapper for the Solver used for CoSimulation"""
    use_base_type = False if not settings.Has("base_type") else (True if settings["base_type"].Has("module") else False)
    if use_base_type:
        base_type_settings = settings["base_type"]
        base_type_module = base_type_settings["module"].GetString()
        full_module_path = ".".join(["KratosMultiphysics.CoSimulationApplication", base_type_module])
        try:
            base_type_imported_module = import_module(full_module_path)
        except ImportError:
            try:
                base_type_imported_module = import_module(base_type_module)
            except ImportError:
                raise ImportError('Base type for coupled solver "{}" could neither be imported from CoSimulation nor from PYTHONPATH'.format(base_type_module))
        if base_type_settings.Has("class"):
            base_type_class = base_type_settings["class"].GetString()
        else:
            module_name = base_type_module.split(".")[-1]
            import string
            base_type_class = string.capwords(module_name.replace("_", " ")).replace(" ", "")
        base_type = getattr(base_type_imported_module, base_type_class)
        return base_factory.Create(settings, [models, solver_name, base_type], "KratosMultiphysics.CoSimulationApplication")
    else:
        return base_factory.Create(settings, [models, solver_name], "KratosMultiphysics.CoSimulationApplication")