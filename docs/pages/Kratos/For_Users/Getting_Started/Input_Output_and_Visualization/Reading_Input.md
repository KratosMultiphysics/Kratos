---
title: Reading ModelPart From Input File
keywords: 
tags: [Python Script Tutorial Reading ModelPart Input File]
sidebar: kratos_for_users
summary: 
---

The `ModelPart` represents an arbitrary part of the `Model` to be simulated and stores the mesh and additional data for it. Most of the *Kratos* routines take a `ModelPart` as their argument. So always is necessary to create and fill a `ModelPart`. In this tutorial we will describe how to create and fill a model part from given input file.

# Setup
First of all we need to create a python file with following code to import the *Kratos*:

```python
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
```
{: data-lang="Python"}

# Creating a ModelPart
To create a `ModelPart`, one has to create a new `Model` first, and then to call its constructor passing the `ModelPart``s name as its argument:

```python
this_model = Model()
fluid_model_part = this_model.CreateModelPart("FluidPart")
```
{: data-lang="Python"}

You can print the fluid_model_part:

```python
>>> print(fluid_model_part)
-FluidPart- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 0
    Current solution step index : 0

    Mesh 0 :
        Number of Nodes      : 0
        Number of Properties : 0
        Number of Elements   : 0
        Number of Conditions : 0
```
{: data-lang="Python Output"}

It can be seen that the `ModelPart` is empty and has the buffer size equal to 1. This means that no history of the nodal solution steps variables will be saved.

The next step is to define the variables we want to store as historical variables in nodes of this `ModelPart` as follow:

```python
fluid_model_part.AddNodalSolutionStepVariable(VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
```
{: data-lang="Python"}

# Reading ModelPart File
The input file of the *Kratos* has `.mdpa` (stand for ModelPart) and contains the properties, nodes, elements, conditions and initial values. A convenient way to create this file is to use the interface prepared for *GiD* pre and post processor. [Here]pages/(Input-data) you can find more information about the input file. Here we assume that the `Cylinder.mdpa` input file is already created using *GiD*:

For reading the `.mdpa` file first we have to create a `ModelPartIO` object passing the input file path/name to its constructor:

```python
fluid_model_part_io = ModelPartIO("path/to/file/example")
```

**NOTE:** the file name is used without the `.mdpa` extension!

Then use this *IO* object to read the file and store the mesh and data in `ModelPart`:

```python
fluid_model_part_io.ReadModelPart(fluid_model_part)
```

And printing again the `ModelPart`:

```console
 >>> print(fluid_model_part)
 StructurePart model part
-FluidPart- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 5
    Current solution step index : 0

    Mesh 0 :
        Number of Nodes      : 3072
        Number of Properties : 2
        Number of Elements   : 5778
        Number of Conditions : 366

    -NoSlip2D_No_Slip_Cylinder- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes      : 126
            Number of Properties : 0
            Number of Elements   : 0
            Number of Conditions : 126
    -Parts_Fluid- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes      : 3072
            Number of Properties : 0
            Number of Elements   : 5778
            Number of Conditions : 0
    -AutomaticInlet2D_Inlet- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes      : 21
            Number of Properties : 0
            Number of Elements   : 0
            Number of Conditions : 20
    -Outlet2D_Outlet- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes      : 21
            Number of Properties : 0
            Number of Elements   : 0
            Number of Conditions : 20
    -NoSlip2D_No_Slip_Walls- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes      : 202
            Number of Properties : 0
            Number of Elements   : 0
            Number of Conditions : 200
```
{: data-lang="Python Output"}

# Setting the Buffer Size
If we need to store the values of the nodal solution step variables in previous time steps, we must modify the buffer size **AFTER** defining the historical variables:

```python
fluid_model_part.SetBufferSize(2)
```

This would store the values for two previous steps in addition to the current ones.

**Next** [Writing Output File](Writing_Output)<br>
**Prev** [Reading ProjectParameters](Project_Parameters)