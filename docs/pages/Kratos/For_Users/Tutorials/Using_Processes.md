---
title: Using processes to customize a simulation
keywords: 
tags: [Using-processes-to-customize-a-simulation.md]
sidebar: kratos_for_users
summary: 
---

# Using `Kratos::Process` to customize a simulation
Kratos offers many ways of customizing the behavior of simulations. One way is to use [`Kratos::Process`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/processes/process.h), either from 
* KratosCore
* Kratos-Application
* user defined 

Adding a process in a simulation is straight-forward, it has to be added to the list of processes in the `ProjectParameters.json`.

## Adding processes in the simulation
In the following it is shown how to add processes from different locations (KratosCore, Kratos-Application or user-defined) to a simulation in `ProjectParameters.json`. Further information can be found [in the description of the AnalysisStage](Common-Python-Interface-of-Applications-for-Users#analysisstage-usage).

```js
"processes" : {
    "from_kratos_core" : [{
        "python_module" : "assign_vector_variable_process",
        "kratos_module" : "KratosMultiphysics",
        "Parameters"    : {
            # ...
        }
    }],
    "from_kratos_application" : [{
        "python_module" : "compute_body_fitted_drag_process",
        "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
        "Parameters"    : {
            # ...
        }
    }],
    "user_defined" : [{
        "python_module" : "user_defined_process_script",
        "Parameters"    : {
            # ...
        }
    }]
}
```


_Note_: In order for the user-defined process to work, it has to be available on the `PYTHONPATH` of the system. This can e.g. be the current working directory.

## Creating a user-defined process

Following is an example containing the functions a python-process can implement. This can be used as a basis for implementing a user-defined process.

In order to see where the functions are being called, one can check the [AnalysisStage](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/analysis_stage.py). 

Since most applications have adapted the AnalysisStage from the KratosCore, it is very easy to use a different application, one does not have to re-learn how the things are working.
```python
import KratosMultiphysics as KM

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DummyPythonProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class DummyPythonProcess(KM.Process):
    """This class is a dummy-process that shows how the functions that can be implemented
    in order to customize the behavior

    Public member variables:
    model -- the container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KM.Process.__init__(self) # calling the baseclass constructor

        default_settings = KM.Parameters("""{
            "default_setting_bool"     : true,
            "default_setting_int"      : 5,
            "default_setting_double"   : -1.45,
            "default_setting_string"   : "dummy",
            "default_setting_array"    : [1.0, 2.5, 4.8],
            "default_setting_subparam" : {
                "level_1" : 123
            }
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        # Alternative:
        # settings.RecursivelyValidateAndAssignDefaults(default_settings)

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def Check(self):
        """ This method verifies that the input is correct

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed just before the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeOutputStep(self):
        """ This method is executed before writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed after writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalize(self):
        """ This method is executed after the computations, at the end of the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass
```
_Note_: Unused functions can be deleted, it is not necessary to have an empty implementation.