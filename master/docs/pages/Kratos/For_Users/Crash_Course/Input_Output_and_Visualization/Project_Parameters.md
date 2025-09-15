---
title: Reading ProjectParameters
keywords: 
tags: [Python Script Tutorial Reading ProjectParameters Project Parameters]
sidebar: kratos_for_users
summary: 
---

The *Kratos* `Parameters` object is a container based on the well known *JavaScript Object Notation* (**JSON**) standard. Even though it can contain any type of key-value information, in *Kratos* it is used to contain configuration settings for solvers, processes or utilities. 

In this tutorial, the use of the **JSON** format together with the *Kratos* `Parameters` class is reviewed using a standard *Kratos* simulation configuration file (`ProjectParameters.json`), in this case coming from the  `FluidDynamicsApplication`, as example.

# Setup
First of all we need to create a python file with following code to import the *Kratos*:

```python
from KratosMultiphysics import *
```

# Reading a **JSON** file
In this subsection we will try to parse the `ProjectParameters.json` file to construct the *Kratos* Parameters object. The `ProjectParameters.json` file reads as follows


```json
{
    "problem_data"                     : {
        "problem_name"    : "parameters_tutorial",
        "model_part_name" : "MainModelPart",
        "domain_size"     : 2,
        "parallel_type"   : "OpenMP",
        "echo_level"      : 0,
        "start_time"      : 0.0,
        "end_time"        : 45
    },
    "output_configuration"             : {
        "result_file_configuration" : {
            "gidpost_flags"       : {
                "GiDPostMode"           : "GiD_PostBinary",
                "WriteDeformedMeshFlag" : "WriteDeformed",
                "WriteConditionsFlag"   : "WriteConditions",
                "MultiFileFlag"         : "SingleFile"
            },
            "file_label"          : "time",
            "output_control_type" : "step",
            "output_frequency"    : 1,
            "body_output"         : true,
            "node_output"         : false,
            "skin_output"         : false,
            "plane_output"        : [],
            "nodal_results"       : ["VELOCITY","PRESSURE"],
            "gauss_point_results" : []
        },
        "point_data_configuration"  : []
    },
    "restart_options"                  : {
        "SaveRestart"      : "False",
        "RestartFrequency" : 0,
        "LoadRestart"      : "False",
        "Restart_Step"     : 0
    },
    "solver_settings"                  : {
        "solver_type"                 : "Monolithic",
        "model_import_settings"       : {
            "input_type"     : "mdpa",
            "input_filename" : "parameters_tuto"
        },
        "echo_level"                  : 0,
        "compute_reactions"           : false,
        "dynamic_tau"                 : 1.0,
        "oss_switch"                  : 0,
        "maximum_iterations"          : 10,
        "relative_velocity_tolerance" : 0.001,
        "absolute_velocity_tolerance" : 1e-5,
        "relative_pressure_tolerance" : 0.001,
        "absolute_pressure_tolerance" : 1e-5,
        "volume_model_part_name"      : "Parts_Fluid",
        "skin_parts"                  : ["AutomaticInlet2D_Inlet","Outlet2D_Outlet","NoSlip2D_No_Slip_Walls","NoSlip2D_No_Slip_Cylinder"],
        "no_skin_parts"               : [],
        "time_stepping"               : {
            "automatic_time_step" : false,
            "time_step"           : 0.1
        }
    },
    "initial_conditions_process_list"  : [],
    "boundary_conditions_process_list" : [{
        "python_module" : "apply_inlet_process",
        "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
        "Parameters"    : {
            "model_part_name" : "AutomaticInlet2D_Inlet",
            "variable_name"   : "VELOCITY",
            "modulus"         : "6*y*(1-y)*sin(pi*t*0.5)",
            "direction"       : "automatic_inwards_normal",
            "interval"        : [0,1]
        }
    },{
        "python_module" : "apply_inlet_process",
        "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
        "Parameters"    : {
            "model_part_name" : "AutomaticInlet2D_Inlet",
            "variable_name"   : "VELOCITY",
            "modulus"         : "6*y*(1-y)",
            "direction"       : "automatic_inwards_normal",
            "interval"        : [1,"End"]
        }
    },{
        "python_module" : "apply_outlet_process",
        "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
        "Parameters"    : {
            "model_part_name"    : "Outlet2D_Outlet",
            "variable_name"      : "PRESSURE",
            "constrained"        : true,
            "value"              : 0.0,
            "hydrostatic_outlet" : false,
            "h_top"              : 0.0
        }
    },{
        "python_module" : "apply_noslip_process",
        "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
        "Parameters"    : {
            "model_part_name" : "NoSlip2D_No_Slip_Walls"
        }
    },{
        "python_module" : "apply_noslip_process",
        "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
        "Parameters"    : {
            "model_part_name" : "NoSlip2D_No_Slip_Cylinder"
        }
    }],
    "gravity"                          : [{
        "python_module" : "assign_vector_by_direction_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"  : "AssignVectorByDirectionProcess",
        "Parameters"    : {
            "model_part_name" : "Parts_Fluid",
            "variable_name"   : "BODY_FORCE",
            "modulus"         : 0.0,
            "constrained"     : false,
            "direction"       : [0.0,-1.0,0.0]
        }
    }],
    "auxiliar_process_list"            : []
}
```

and can be parsed to construct a *Kratos* Parameters object with the next two lines of code

```python
json_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters(json_file.read())
```

To visualize the content of the *ProjectParameters* we need to make the *ProjectParameters* printable. This can be done by calling the **PrettyPrintJsonString** method from the *ProjectParameters* object in this way

```python
print(ProjectParameters.PrettyPrintJsonString())
```

# Working with the *Kratos* Parameters
Once we have parsed the `ProjectParameters.json` file, we can start to check, get and edit its information. At this point it is interesting to mention that the Parameters works in a similar manner that a common Python dictionary does. For instance we can extract the solver settings required by the Python solvers by doing

```python
solver_settings = ProjectParameters["solver_settings"]
```

Similarly, we can do the same operation for a list entry. For instance we can iterate through the entire list of boundary conditions settings as is done below. Note that we have used the method `size` to obtain the length of the iterated list.

```python
for i in range(ProjectParameters["boundary_conditions_process_list"].size()):
    boundary_condition_settings = ProjectParameters["boundary_conditions_process_list"][i]
```

Complementary, we can check if any field exists before trying to retrieve its value with the `Has` method. This method returns a boolean variable with value `True` if the field exists and false otherwise.

```python
ProjectParameters.Has("output_configuration")
```

To get or modify the value of any field, the *Kratos* Parameters incorporates the `Get` and `Set` methods, which are particularized for all the variable types. Some example of its usage can be found in the lines below.

```python
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()
domain_size = ProjectParameters["problem_data"]["domain_size"].GetInt()
ProjectParameters["problem_data"]["end_time"].SetDouble(20.0)
ProjectParameters["problem_data"]["model_part_name"].SetString("NewMainModelPart")
```

Since this is a basic tutorial on the use of the *Kratos* Parameters object, only the basic features have been described.  For more advanced operations check the **JSON** configuration file tutorial in ([here]pages/(How-to-write-a-JSON-configuration-file) )

**Next** [Reading ModelPart From Input File](Reading_Input)<br>
**Prev** [Hello Kratos](../Hello_World)