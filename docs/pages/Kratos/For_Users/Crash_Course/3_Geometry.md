---
title: Kratos Input Files and IO
keywords: 
tags: [Kratos-input-files-and-IO.md]
sidebar: kratos_for_users
summary: 
---

## 1. The Model and ModePart

The following section is a short version of a more detailed [tutorial](Python-Script-Tutorial:-Reading-ModelPart-From-Input-File).

### 1.1 Create A Model and ModelPart

As you have seen before the Model and the modelpart are the data structures responsible for storing the geometrical information about your simulation.

A `ModelPart` has to be created via a `Model` object, which can contain several `ModelParts`. Right after creating the empty `ModelPart`, the variables needed for the following calculations have to be added. The empty `ModelPart` is then filled using the `ModelPartIO` that reads the information from an .mdpa file. In general, the analysis object takes care of these steps, especially because it knows which variables to add.

Generaly, you will not have to deal with the lecture of a `ModelPart`, but here you will do it directly in your python script. 

If the `.mdpa` file contains application dependent elements, the corresponding Kratos application has to be imported. In our structural example, the elements are from the `StructuralMechanicsApplication`. 

Lets create a small python script `tutorial_modelpart.py` with the following lines:

```python
import KratosMultiphysics.StructuralMechanicsApplication

model = KratosMultiphysics.Model()
model_part_1 = model.CreateModelPart("MyModelPart1")
model_part_2 = model.CreateModelPart("MyModelPart2")
```

This is similar to what we had in the first part, we import an application and we create a `Model`, but we have added a new line that calls a method from the `Model` called `CreateModelPart`, which accepts a name as an argument.

We can now print both this model and modelpart to see what they contain:

```python
# Print the model
print(model)
```

Which will print the information of the whole `model`, containing both `model_part_1` and `model_part_2`:

```bash
-MyModelPart1- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 0
    Current solution step index : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0

-MyModelPart2- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 0
    Current solution step index : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0 
```

And we can also print the info of a single `ModelPart`:

```python
# Print the model_part_1
print(model_part_1)
```

Which will only print the information of the `model_part_1`:

```bash
-MyModelPart1- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 0
    Current solution step index : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0
```

### 1.2 Create a Database

We have created a `Model` with various `ModelPart`, but those are essentially empty. In order fro them to be usefull we will need fill it with the gemoetrical information and the physical data.

The first we must do is to tell our modelparts which phsical data will be storing inside. In Kratos, for performance, this needs to be known before reading the actual geometry and is set individually for every different `ModelPart`. Will explain it further in the data management section. 

For now, lets imagine that we want that our model holds the `DISPLACEMENT` of our nodes over time, we can tell it to the model with this instruction:

```python
# Adding variables BEFORE reading the .mdpa
model_part_1.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
```

There are many [variables](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python/add_containers_to_python.cpp) that we can use in the core and in applications (ex: [structural_mechanics_application.h](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/StructuralMechanicsApplication/custom_python/structural_mechanics_python_application.cpp)), and if none of them suits your needs you can always deifne your own.

Be aware that this step of of crucial important, specially if you pretemd to read a mpda file that containts entity info. If you fail to tell the model_part which variable are available, you will se an error like this while trying to read or run your simulation:

```bash
```

As stated in the Basics section, the `Solver` is nomraly in charge of reading the `ModelPart`, and the reason for that is in fact that the `Solver` knows which of those variables are used in the simulation.

### 1.3 Read a .mdpa file

We have initialized our model and modelpart and assigned the variables we will use. Its now time to actually read the modelpart.

This can be done through the `ModelPartIO` class. You only have to specigy the name of the file **without the .mdpa extension** and then call the read function with the model_part you want to be the one holding the file data. In this case we will go with `KratosWorkshop2019_high_rise_building_CSM` and `model_part_1`

```python
# Create a ModelPartIO and use the Read function
model_part_io = KratosMultiphysics.ModelPartIO("KratosWorkshop2019_high_rise_building_CSM")
model_part_io.ReadModelPart(model_part_1)
```

You can of course print the `model_part_1` to check that the lecture has been done sucessfully:

```python
# Print the model_part_1
print(model_part_1)
```

Will now yield:
```bash
-MyModelPart1- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 3
    Current solution step index : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 273
        Number of Properties  : 1
        Number of Elements    : 228
        Number of Conditions  : 82
        Number of Constraints : 0

    -DISPLACEMENT_Ground- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 7
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 0
            Number of Constraints : 0
    -LineLoad2D_InterfaceStructure- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 83
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 82
            Number of Constraints : 0
    -Parts_Structure- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 273
            Number of Properties  : 0
            Number of Elements    : 228
            Number of Conditions  : 0
            Number of Constraints : 0
```

### 1.4 Manipulating the ModelPart

While you can read whole files, it is also possible to assign its entities manually. 

A modelpart is generally constructed using `Nodes`, `Elements` and `Conditions`, and those entities can be created directly in python.

Let's try to create a simple square using triangular elements. First we need to create the nodes, we will need 4 of them and we will add them to the `model_part_2` with the `CreateNewNode(ID, X, Y, Z)` function:

```python
node_1 = model_part_2.CreateNewNode(1, 0.00000, 0.00000, 0.00000)
node_2 = model_part_2.CreateNewNode(2, 1.00000, 0.00000, 0.00000)
node_3 = model_part_2.CreateNewNode(3, 1.00000, 1.00000, 0.00000)
node_4 = model_part_2.CreateNewNode(4, 1.00000, 0.00000, 0.00000)
```

Once we have all nodes, we can add a couple of elements.

The signature for the function that creates elements is `CreateElement(Name, ID, [LIST OF NODES ID], PROPERTIES)`.

This function is a bit more complicated as Kratos has different concepts for the class `Element` which holds the phisical properties, and the `Geometry` which holds the spatical information. If you explore Kratos files you will se that thipically this is shown in the name, for example:

- `Element2D2N`: Means that we are using the generic `Element`, with two spatial coordinates `2D`, with two nodes `2N` (a line)
- `LaplacianElement3D4N`: Means that we are using a `LaplacianElement`, with three spatial coordinates `3D`, with 3 nodes `4N` (a teth)

Since we want a triangle and we don't care about the Z coordinates, we can use an `Element2D3N`

The properties are just another database that holds information but which is shared for all entitites of the modelpart. For the properties field we will just use `model_part_2.GetProperties(1)`. This will create empty properties for us. 

```python
element_1 = model_part_2.CreateNewElement("Element2D3N", 1, [1,3,2], model_part_2.GetProperties(1))
element_2 = model_part_2.CreateNewElement("Element2D3N", 2, [1,4,3], model_part_2.GetProperties(1))
```

⚠️ Be mindful when asigning the list of nodes, as the order will determine the normal of the element if its a surface.

We can check that out modelpart has the newly added info by printing it, as typical:

```bash
-MyModelPart2- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 0
    Current solution step index : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 4
        Number of Properties  : 1
        Number of Elements    : 2
        Number of Conditions  : 0
        Number of Constraints : 0
```

### 1.5 Submodelparts

Aside for the possibility of a `Model` to contain several `ModelParts`, we also have the hability to create subdivisions inside a `ModelPart` which are called `SubModelPart`. 

Each submodelpart represents a subset of the entitis of a given `ModelPart` that for whatever reason are interesant to keep grouped. For example, the elements and nodes that belong to an inlet of a wind tunnel can be a submodelpart of the geometry representinf the wind tunnel.

It is very important to notice that: 
1) While the different modelParts in a model a essentially different objects, submodelparts are only sets of existing entities. New entitites will not be created when added to a submodelpart. 
2) Entitites can belong to different submodelparts
3) Submodelpart can have other submodelparts, hence, providing you with a mechanism to heriarchicly divide your modelparts. but an entity that belong to a given submodelpart will also belong to all its parents.

## 2. Output
The `vtk_output` block in the ProjectParameters.json gives you an impression on the potential settings for the output. Here you will create just a minimal version of it.

In this part of the tutorial you will create a minimal configuration of a VTK output process. 

```python
from vtk_output_process import VtkOutputProcess
vtk_output_configuration = KratosMultiphysics.Parameters("""{
        "model_part_name"        : \""""+this_model_part.Name+"""\",
        "output_sub_model_parts" : false,
        "nodal_solution_step_data_variables" : ["DISPLACEMENT"]
    }""")

vtk_output = VtkOutputProcess(this_model, vtk_output_configuration)
```

The output process is usually called at defined places inside the analysis. In order to use it, several functions need to be called in the right order.

```python
vtk_output.ExecuteInitialize()
vtk_output.ExecuteBeforeSolutionLoop()
vtk_output.ExecuteInitializeSolutionStep()
vtk_output.PrintOutput()
vtk_output.ExecuteFinalizeSolutionStep()
vtk_output.ExecuteFinalize()
```
