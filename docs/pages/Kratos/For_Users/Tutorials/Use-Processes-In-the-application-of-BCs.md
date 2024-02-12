---
title: Processes In the application of BCs
keywords: 
tags: [How-to-use-Processes-In-the-application-of-BCs.md]
sidebar: kratos_for_users
summary: 
---

# Description of Main Processes for the application of boundary conditions

Kratos implements a few processes for the application of boundary conditions, which make an attempt to cover the most common user cases.
The list of such key processes is the following:

- AssignScalarVariableProcess
- AssignVectorVariableProcess
- AssignVectorByDirectionProcess

let's focus on each of them:

# AssignScalarVariableProcess
It is designed to apply boundary conditions to all of the nodes on a given submodelpart.
The process is constructed by providing a "Parameters" object, initialized by a string in json format.

The minimal code needed to use the process from the Project parameters is the following block of code:

    {
        "python_module" : "assign_scalar_variable_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"  : "AssignScalarVariableProcess",
        "Parameters"    : {
            "model_part_name" : "Main",
            "variable_name"   : "DISPLACEMENT_X",
            "value"           : "sqrt(x**2+y**2)*t",
        }
    }
    
of which the first three parameters are fixed and needed for the automatic construction of the process. 

The user can specify the behaviour of the module by adjusting the content of the Paramters section of the STRING

- **model_part_name** -- identifies the submodelpart on which the process will be applied
- **variable_name** allows to specify the variable to be prescribed
- **value** is the place to specify the value to be applied. Such field can be either a double number or a simple mathematical expression written in terms of the variables "x y z t". If a function is prescribed, python math notation will be used.

If no other options are prescribed the variable will be "fixed" (that is its value constrained not to change). Such status will be mantained for the entire simulation.

More options are however available to tune the behaviour of the process. Namely the user can prescribe a more complete input parameter of the type

    {
        "python_module"   : "assign_scalar_variable_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"          : "AssignScalarVariableProcess",
        "Parameters"            : {
            "model_part_name" : "Main",
            "variable_name"   : "DISPLACEMENT_X",
            "interval"        : [0.0, 5.0],
            "constrained"     : true,
            "value"           : "sqrt(x**2+y**2)*t",
            "local_axes"      :{
                "origin" : [0.0, 0.0, 0.0],
                "axes"  : [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
            }
        }
    }
    
the optional extra parameters have the following meaning
- **interval** is the time span during the simulation process in which the variable will be prescribed. Fixity will be imposed (if required, see next items) only during such time interval. Similarly the value will be applied only during the time period prescribed by the interval variable.
- **constrained** allows controlling wether the value is to be fixed, or if its value is simply to be changed (without Fixing it)
- **local_axes** allows to prescribe the system of coordinates in which the function is to be evaluated. The parameters "x y z" to be used in the evaluation of the mathematical expressions, will in this case be evaluated with respect to an origin located at **origin** considerin the basis formed by the prescribed **axes**. That is, if (xo,yo,z0) are the coordinates of the origin and V1,V2,V3 are the basis vector, the x actually used in the formula will be computed as x = dot(x-xo, V1)

an important feature is that in the case in which "constrained"=True the value will be fixed at the beginning of the step and be left free at the end of the step. This is needed to ensure that the fixity is only applied within the interval. 

**TODO** the possibility of considering input as a table will be added in the future

# AssignVectorVariableProcess
Such process attempts to generalize the AssignScalarVariableProcess to the case of vector values. The process uses internally the scalar version and hence inherits all of its capabilities as well as, when possible, the same interface
Some difference are however present in the application of the component values and fixity.

A quite general example of application of the process is as follows:

    {
        "python_module"   : "assign_vector_variable_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"          : "AssignVectorProcess",
        "Parameters"            : {
                "model_part_name"      : "Main",
                "variable_name"        : "DISPLACEMENT",
                "interval"             : [25.0, "End"],
                "constrained"          : [true,true,false],
                "value"                : [null, "x+y*t", 10.0],
                "local_axes"           : {}
            }
    }
    
The most remarkable difference from the scalar case is that now **value** is a tuple, in which each component can be either a function or a constant value. The value "null" (without parenthesis) has in this context a special meaning: it means that the null component should be completely ignored in the application of the process.
Two formats are admitted for **constrained**. A vector in the form [true,true,false], implies that the component x and y should be fixed and the z component should be only applied without fixing its value. Alternatively the user can use a single bool value for the 3 components, thus telling that the fixity is the same for all the three components and is controlled by the input value (which defaults to true)


# AssignVectorByDirectionProcess
This is a process similar to the previous, in which however the 3 components are set at once by providing direction and modulus.

    {
        "python_module"   : "assign_vector_by_direction_process",
        "kratos_module" : "KratosMultiphysics",
        "process_name"          : "AssignVectorByDirectionProcess",
        "Parameters"            : {
                "model_part_name"      : "Main",
                "variable_name"        : "VELOCITY",
                "interval"             : [25.0, "End"],
                "modulus"              : "sqrt(abs(x*y))",
                "constrained"          : true,
                "direction"            : [0.0, 1.0, 1.0],
                "local_axes"           : {}
            }
    }
    
- **direction** indicates the orientation in which the variable should be applied. Note that such direction will be normalized internally, so its norm must be non zero.
- **modulus** admits a function or a constant value and tells the modulus of the quantity to be applied
- **constrained** only bool values are admissible. Constraints at once (or frees at once) all the components

# Deprecated Processes
This list shows which of the above processes replace deprecated processes:
* ImposeScalarValueProcess => AssignScalarVariableProcess
* ImposeVectorValueByComponentsProcess => AssignVectorVariableProcess
* ImposeVectorValueByDirectionProcess => AssignVectorByDirectionProcess
* ApplyConstantValueProcess => AssignScalarVariableProcess
* ApplyCustomFunctionProcess => AssignScalarVariableProcess / AssignVectorVariableProcess
* ApplyVariableValueProcess => AssignScalarVariableProcess / AssignVectorVariableProcess