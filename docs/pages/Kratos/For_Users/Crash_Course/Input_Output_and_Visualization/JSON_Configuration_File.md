---
title: JSON Configuration File
keywords: 
tags: [Write JSON Configuration File]
sidebar: kratos_for_users
summary: 
---

## Overview

As complexity increases one needs more and more flexibility in “configuring” a simulation. This is typically achieved by defining a "configuration file" which organizes all of the "user instructions" to be used in the configuration of the solver. 

Some requirements for this configuration are:

*    Must allow nesting of parameters (tree-like structure)
*    Should allow using default values
*    Must allow some form of type checking
*    Should permit to echo all of the defaults in a single place 

In order to meet such requirements Kratos introduces a ``Parameters`` object, which provides a thin wrapper to JSON strings.

## What is JSON?

<div style="text-align:center;">
    <img style="width:120px;" src="https://upload.wikimedia.org/wikipedia/commons/c/c9/JSON_vector_logo.svg">
</div>

In computing, [JSON](https://en.wikipedia.org/wiki/JSON) is an open-standard format that uses human-readable text to transmit data objects consisting of attribute–value pairs. It is the most common data format used for asynchronous browser/server communication, largely replacing XML, and is used by AJAX. JSON is a language-independent data format.

## Use of JSON

The typical structure of a JSON file is like the following: 

```json
{
    "problem_data"             : {
        "problem_name"    : "minimal_structure",
        "start_time"      : 0.0,
        "end_time"        : 1.0,
        "echo_level"      : 0
    },
    "solver_settings"          : {
        "solver_type"                        : "static",
        "echo_level"                         : 0,
        "analysis_type"                      : "non_linear",
        "model_part_name" : "Structure",
        "domain_size"     : 2,
        "model_import_settings"              : {
            "input_type"     : "mdpa",
            "input_filename" : "minimal_structure"
        },
        "time_stepping"  : {
            "time_step"       : 1.1
        },
        "line_search"                        : false,
        "convergence_criterion"              : "residual_criterion",
        "displacement_relative_tolerance"    : 0.0001,
        "displacement_absolute_tolerance"    : 1e-9,
        "residual_relative_tolerance"        : 0.0001,
        "residual_absolute_tolerance"        : 1e-9,
        "max_iteration"                      : 10,
        "linear_solver_settings"             : {
            "solver_type" : "Super_LU",
            "scaling"     : false,
            "verbosity"   : 0
        },
        "problem_domain_sub_model_part_list" : ["Parts_Parts_Auto2"],
        "processes_sub_model_part_list"      : ["DISPLACEMENT_Displacement_Auto1","SelfWeight2D_Self_weight_Auto1"],
        "rotation_dofs"                      : false
    },
    "processes" : {
        "constraints_process_list" : [{
            "implemented_in_file"   : "impose_vector_value_by_components_process",
            "implemented_in_module" : "KratosMultiphysics",
            "help"                  : "This process fixes the selected components of a given vector variable",
            "process_name"          : "ImposeVectorValueByComponentsProcess",
            "Parameters"            : {
                "mesh_id"         : 0,
                "model_part_name" : "DISPLACEMENT_Displacement_Auto1",
                "variable_name"   : "DISPLACEMENT",
                "is_fixed_x"      : true,
                "is_fixed_y"      : true,
                "is_fixed_z"      : true,
                "value"           : [0.0,0.0,0.0]
            }
        }],
        "loads_process_list"       : [{
            "implemented_in_file"   : "process_factory",
            "implemented_in_module" : "KratosMultiphysics",
            "check"                 : "DirectorVectorNonZero direction",
            "help"                  : "This process ",
            "process_name"          : "ApplyConstantVectorValueProcess",
            "Parameters"            : {
                "mesh_id"         : 0,
                "model_part_name" : "SelfWeight2D_Self_weight_Auto1",
                "variable_name"   : "VOLUME_ACCELERATION",
                "factor"          : 9.8,
                "direction"       : [10.0,0.0,0.0]
            }
        }]
   },    
   "output_processes" : {
        "gid_output" : [{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "Structure.computing_domain",
                "output_name"            : "minimal_structure",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"       : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"          : "step",
                        "output_control_type" : "step",
                        "output_frequency"    : 4,
                        "body_output"         : true,
                        "node_output"         : false,
                        "skin_output"         : false,
                        "plane_output"        : [],
                        "nodal_results"       : ["DISPLACEMENT","REACTION"],
                        "gauss_point_results" : ["VON_MISES_STRESS"]
                    },
                    "point_data_configuration"  : []
                }
            }
        }]
    }
}
```

You can check your JSON file usin a JSON validator, [like](http://jsonlint.com/) .

Such text shall be read to an ``input_string`` which is then parsed and stored within a ``Parameters`` object

this is done in python as:

```python
input_string = open("ProjectParameters.json",'r').read()
settings = Parameters( input_string )
``` 

The settings object now behaves effectively as a combination of lists and dictionaries, for example one can obtain the name of the problem to be solved as 

```python
 name = settings["problem_data"]["problem_name"].GetString()
```

An important aspect of this is that the user is expected to ask specifically for the type to be provided For the case at hand (the name) one is asking for the item ``settings["problem_data"]["problem_name"]`` to be provided as a string by calling the ``GetString()`` method. 

Similarly one can call: 

* ``GetInt()/SetInt(...)``
* ``GetDouble()/SetDouble(...)``
* ``GetBool()/SetBool(...)``
* ``GetString()/SetString(...) ``

The type can be queried by the functions:

* ``IsInt()``
* ``IsDouble()``
* ``IsBool()``
* ``IsString()``
* ``IsArray()``
* ``IsSubParameter()`` 

In the case of an array item the function ``size()`` is provided, and the different entries can be accessed by the [] operator. Note that each entry may have a different type and the ``GetXXX()`` method shall be used even for array entries. 

New values can be added to the ``settings`` by doing:

```python
settings.AddEmptyValue("item_name")
settings["item_name"].SetDouble(1.0) #this would set the new item to the value of 1.0
``` 

Entire parameter lists can be addeded by the function:

```python
settings.AddValue("subparameter_name", other_parameters_object)
``` 

At any point one can obtain a formatted echo of the ``settings`` by doing:

```python
print(settings.PrettyPrintJsonString())
``` 

A compact version, without any extra character (useful for example for serialization) can be obtained by doing:

```python
settings.WriteJsonString()
``` 

## Style Conventions for JSON
Information about style conventions for the JSON file can be found in the [Style Guide](../../Conventions/Style_Guide).
