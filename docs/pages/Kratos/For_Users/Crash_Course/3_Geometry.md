---
title: Kratos Input Files and IO
keywords: 
tags: [Kratos-input-files-and-IO.md]
sidebar: kratos_for_users
summary: 
---

# 1. The Model and ModePart

As previously discussed, the `Model` and `ModelPart` are essential data structures in Kratos, responsible for storing the geometrical and structural information for your simulation.

## 1.1 Create A Model and ModelPart

A `ModelPart` has to be created via a `Model` object, which can contain several `ModelParts`. Right after creating the empty `ModelPart`, the variables needed for the following calculations have to be added. The empty `ModelPart` is then filled using the `ModelPartIO` that reads the information from an .mdpa file.

If the `.mdpa` file contains application dependent elements, the corresponding Kratos application has to be imported. In our structural example, the elements are from the `StructuralMechanicsApplication`. 

 Typically, these steps are handled by the `Solver` object, which knows which variables and components to load based on the simulation configuration. However, to understand the process more directly, let’s create a simple Python script, `tutorial_modelpart.py`, that demonstrates how to read a ModelPart manually:

```python
import KratosMultiphysics.StructuralMechanicsApplication

model = KratosMultiphysics.Model()
model_part_1 = model.CreateModelPart("MyModelPart1")
model_part_2 = model.CreateModelPart("MyModelPart2")
```

In this example, we are importing an application, creating a `Model`, and calling a method, `CreateModelPart`, which takes a name as an argument to create a `ModelPart`. This additional step allows us to build a dedicated section in the model to store `nodes`, `elements`, and related data structures, specific to our simulation needs.

To understand how these objects are structured and what information they contain, we can add print statements to examine the contents of both the `Model` and `ModelPart`:

```python
# Print the model
print(model)
```

will print the information of the whole `model`, containing both `model_part_1` and `model_part_2`:

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

will only print the information of the `model_part_1`:

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

## 1.2 Create a Database

We have created a `Model` with various `ModelPart`, but those are essentially empty. In order fro them to be usefull we will need fill it with the gemoetrical information and the physical data.

The first essential step in setting up `ModelParts` in Kratos is to specify the type of physical data they will store. In Kratos, the structure and type of data stored in a `ModelPart` must be defined in advance for performance reasons. This setup allows Kratos to allocate resources efficiently, improving computational speed and memory management.

Each `ModelPart` is configured to store specific physical quantities, known as `Variables` such as `DISPLACEMET`, `VELOCITY`, or `FORCE`. This configuration is required before reading in the geometry or other model data, and it must be performed separately for each `ModelPart` you create.

For now, lets imagine that we want that our model holds the `DISPLACEMENT` of our nodes over time, we can tell it to the model with this instruction:

```python
# Adding variables BEFORE reading the .mdpa
model_part_1.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
```

There are many [variables](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python/add_containers_to_python.cpp) that we can use in the core and in applications (ex: [structural_mechanics_application.h](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/StructuralMechanicsApplication/custom_python/structural_mechanics_python_application.cpp)), and if none of them suits your needs you can always deifne your own.

Be aware that this step is of crucial important, specially if you pretend to read a mpda file that containts entity info. If you fail to tell the `ModelPart` which variables are available, you will se an error while trying to read or run your simulation.

As stated in the Basics section, the `Solver` is normaly in charge of reading the `ModelPart`, and the reason for that is in fact that the `Solver` is the class that knows which of those `Variables`` are used in the simulation.

## 1.3 Read a .mdpa file

We have initialized our model and modelpart and assigned the variables we will use. Its now time to actually read the modelpart.

The process of reading data into a `ModelPart` is handled by the `ModelPartIO` class in Kratos. This class is designed to streamline the loading of model data from `.mdpa` files. To use it, you need to specify the name of the file **excluding the .mdpa extension** and then call the `ReadModelPart` function, passing the target `ModelPart` as an argument. In this case we will go with `KratosWorkshop2019_high_rise_building_CSM` and `model_part_1`:

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

## 1.4 Adding Entities Manually

While you can read whole files, it is also possible to assign its entities manually. 

A `ModelPart` is generally constructed using `Nodes`, `Elements`, `Conditions` and `MasterSlaveConstraints`, and those entities can be created directly in python. 

Let's see how:

### 1.4.1. Nodes
Let's try to create a simple square using triangular elements. First we need to create the nodes, we will need 4 of them and we will add them to the `model_part_2` with the `CreateNewNode(ID, X, Y, Z)` function:

```python
node_1 = model_part_2.CreateNewNode(1, 0.00, 0.00, 0.00)
node_2 = model_part_2.CreateNewNode(2, 1.00, 0.00, 0.00)
node_3 = model_part_2.CreateNewNode(3, 1.00, 1.00, 0.00)
node_4 = model_part_2.CreateNewNode(4, 1.00, 0.00, 0.00)
```

### 1.4.2. Elements and Conditions

Once we have all nodes, we can add a couple of elements or conditions. Both entities have essentially the same signature, so we will only cover elements.

The signature for the function that creates elements is `CreateElement(Name, ID, [LIST OF NODES ID], PROPERTIES)`.

This function is a bit more complicated as Kratos has different concepts for the class `Element` which holds the phisical properties, and the `Geometry` which holds the spatical information. If you explore Kratos files you will se that thipically this is shown in the name, for example:

- `Element2D2N`: Means that we are using the generic `Element`, with two spatial coordinates `2D`, with two nodes `2N` (a line)
- `LaplacianElement3D4N`: Means that we are using a `LaplacianElement`, with three spatial coordinates `3D`, with 3 nodes `4N` (a teth)

Since we want a triangle and we don't care about the Z coordinates, we can use an `Element2D3N`

The properties are just another database that holds information but which is shared for all entitites of the `ModelPart`. For the properties field we will just use `model_part_2.GetProperties(1)`. This will create empty properties for us. 

```python
element_1 = model_part_2.CreateNewElement("Element2D3N", 1, [1,3,2], model_part_2.GetProperties(1))
element_2 = model_part_2.CreateNewElement("Element2D3N", 2, [1,4,3], model_part_2.GetProperties(1))
```

⚠️ Be mindful when asigning the list of nodes, as the order will determine the normal of the element if it is a surface.

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

### 1.4.3. MasterSlaveConstraints

This is the most advanced entity you can add to a modelpart and will represent a constrait that needs to hold while solving the problem. It requeires not only for existing entities but also for `Dofs` (degrees of freedom) which are normally managed by the solver.

You will see more detailed information in future sections, but a sneek peak of how would look like would be something like this:

```python
import KratosMultiphyscis

# Creating the model and a model_part
model = KratosMultiphysics.Model()
model_part = model.CreateModelPart("constraint_example")

# Setting up the database
model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

# Creating some nodes
node_1 = model_part.CreateNewNode(1, 0.00000, 0.00000, 0.00000)
node_2 = model_part.CreateNewNode(2, 0.00000, 1.00000, 0.00000)

# Define some dofs (you haven't seen this yet, don't worry)
KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)

# Add the MasterSlaveConstraint Entity
mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", id=1, 
    node_1, KratosMultiphysics.DISPLACEMENT_X, 
    node_2, KratosMultiphysics.DISPLACEMENT_X, 
    1.0, 0
)
```

As with the `Elements` and `Conditions`, the type (in this case `"LinearMasterSlaveConstraint"`) will determine its physcial properties.

## 1.5 Submodelparts

In addition to the capability of a `Model` to contain multiple `ModelParts`, Kratos provides the functionality to create subdivisions within a `ModelPart`, referred to as `SubModelParts`. This hierarchical organization allows for more efficient data management and analysis within complex models.

A `SubModelPart` represents a subset of the entities (`nodes`, `elements`, etc.) within a `ModelPart` that are grouped together for specific purposes. For instance, `elements` and `nodes` associated with the inlet of a wind tunnel can be organized into a `SubModelPart` to facilitate targeted simulations or analyses related to airflow through that inlet.

It is very important to notice that: 

1) While the different `ModelParts` in a `Model` a essentially different objects, `SubModelParts` are only sets of existing entities. New entitites will not be created when added to a `SubModelPart`. 
2) Entitites can belong to different `SubModelParts`
3) `SubModelPart` can have other `SubModelParts`, hence, providing you with a mechanism to heriarchicly divide your `ModelParts`. but an entity that belong to a given `SubModelPart` will also belong to all its parents.

# 2. Output
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
