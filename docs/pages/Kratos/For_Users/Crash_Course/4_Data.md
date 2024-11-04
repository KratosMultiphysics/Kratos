---
title: Kratos Data Management
keywords: 
tags: [Kratos-data-data-management]
sidebar: kratos_for_users
summary: 
---

## 1. Introduction

In the past section we have seen how to manipulate our geometry and perform several actions such as creating nodes, elements and writing or reading them.
If we go back to section [1.2 Creating the database]() there was one action that needed to be done before performing this operation and that was adding our variables.

In this section we will cover in detail this process and show how Kratos handles the databases.

## 2. Data typologies

When speaking about data in kratos we typically reffer to variables and we normally distigush three different types of data that one may use in a simulation.

## 2.1 Regular Variables

We understand by regular variables the ones offered by the programing language. In this case anything provided by python or C++ such an integer, a double or a string. 

These variables are used to represent generic properties that do not represent a physical magnitude or are common enough that are usefull to have in a generic way. 

For example, the coordinates of a node are a triplet of doubles that represents the physical position of that node in the space:

We can print the information about one of the nodes that we created in the previous section here:

```python
print(f"Node Info | ID: {node_3.Id}, X: {node_3.X}")
```

## 2.2. Kratos Variables

On the contraty, we understand by a Kratos variables a space that the Kratos frameworks provides to you in which you can store a value and are tied to a node, element, condition or similar. While you can create your own, there is a set of predefined variables that you can find in X and Y.

The main advantage of this variables is that can be mixed and shared among all components and applications in Kratos and its one of the pillars that allow the multiphyscis capabilites of Kratos.

In other words, this means that if your Strcutrual Solver modifies the value of the pressure in one node, everyone in Kratos will be able to account for that change, without explictly needing to communicate it.

Kratos variables are divided into two main categories:

### 2.2.1. Non Historical Variables: 

The non historical variables hold information for information that typically does not change over time, or that if it does, does it in a way that renders the old values obsolete.

These variables can be used at any point in time without needing to tell the Model which ones are needed before the simulation starts. typicaly in kratos `Nodes`, `Elements`. `Conditions`, `Properties` and `ProcessInfo` can have non-historical variables.

Non Historical Variables can be interacted with their own providded `GetValue` and `SetValue` functions. For example, setting and getting a value of PRESSURE of a node:

```python
import KratosMultiphyscis as KMP

# ...

node_3.SetValue(KMP.PRESSURE, 1.5)

print(node.GetValue(KMP.PRESSURE))
``` 

will return

```
> 1.5
```
        

### 2.2.2. Historical Variables: 
Historical variables on the other hand hold not only information about the their current value but also about the value they had in previous time steps of the simulation. The number of historical values that different historical variables have is the same and can be set with the `SetBufferSize` for every modelpart.

For performance reasons the list of historical variables needs to be known before the simulation starts (and for that matter, before reading a mdpa) and cannot change. The historical variables are currently only allowed in the `nodes`.

As with the non-historical variables, there are a set of functions provided to interact with them: 
- **`SetBufferSize`**: To specify the size of the buffer (by default 1). Note that this needs to be set inside the modelpart, and no for individual entities. 

    ```Python
    modelpart = Model.CreateModelPart("Test")
    modelpart.SetBufferSize(3) # Indicates that I want to keep the value for the last three time steps
    ```

- **`AddNodalSolutionStepVariable`**: To inform the `ModelPart` that a variable will need to have an historical database. Also note that this is done at the `ModelPart` level.

    ```Python
    modelpart = Model.CreateModelPart("Test")
    modelpart.AddNodalSolutionStepVariable(KMP.PRESSURE) # Indicates that I want to use PRESSURE as historical variable
    ```

- **`SetSolutionStepValue`**: Sets a value to the entity containing the variable. If not specified it will do it by adding it to the current time step, but that can be changed specifying how many times steps back you want to set the value:

    ```Python
    node_3.SetSolutionStepValue(KMP.PRESSURE, 2.0)    # Sets PRESSURE to 2 in the current time step
    node_3.SetSolutionStepValue(KMP.PRESSURE, 2.0, 0) # Sets PRESSURE to 3 in the time step 0 (current time step)
    node_3.SetSolutionStepValue(KMP.PRESSURE, 4.0, 1) # Sets PRESSURE to 4 in the time step 1 (previous time step)
    ```

- **`GetSolutionStepValue`**: Gets the variable in the current entity. If not specified it will do it by adding it to the current time step, but that can be changed specifying how many times steps back you want to set the value:

    ```Python
    print(node_3.GetSolutionStepValue(KMP.PRESSURE))    # Gets PRESSURE in the current time step
    print(node_3.GetSolutionStepValue(KMP.PRESSURE, 0)) # Gets PRESSURE in the time step 0 (current time step)
    print(node_3.GetSolutionStepValue(KMP.PRESSURE, 1)) # Gets PRESSURE in the time step 1 (previous time step)
    ```

    will return:

    ```Python
    > 2.0
    > 2.0
    > 4.0
    ```

## 2 The Node






























## 1. Introduction

In the first tutorial you have successfully installed Kratos and confirmed it works, but in order to run a real simulation you will need some files that define your problem. This sections aims to cover which are the essential files that you need to run a full fledged simulation and start to play with some of the Kratos base structures.

There are three different types of files that compose a Kratos simulation:
* `.py`: A python script to run Kratos
* `.json`: Contains settings for Kratos
* `.mdpa`: Contains the model part information

These files are usually created by a GUI (GiD, Salome, Flowgraph) and used to directly start the simulation. As the aim of the course is not to learn how to use these tools, you can download the input files for a simple structural mechanics case [here](https://github.com/KratosMultiphysics/Documentation/raw/refs/heads/master/Workshops_files/Kratos_Workshop_2019/Sources/2_Kratos_input_files_and_IO/2_Kratos_input_files_and_IO.zip).**

## 2. The Kratos python script
### 2.1. MainKratos.py
The most important file that you have downloade is called `MainKratos.py` and as its name suggest is a python script. This would be the equivalent to your well known `.exe` for a classic problem, with the advantage that you can customize it according to your needs. It is responsible to load the required Kratos applications and to call the main Kratos functionalities as desired by the user. 

Let's look at the content of this file for our structural analysis example:

```python
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Run()
```

In the first lines, Kratos and the structural analysis are imported. Then the settings are read from the `.json` and a Model is created. Finaly we use all that information to create a `StructuralMechanicsAnalysis` simulation. In the last line, the structural simulation executed. 

It's important to notice that, while this script works and its correct, it is not the only correct way to initialize a Kratos simluation. As we will see by the end of the Crash course, there are more advanced usages.

### 2.2. Run Kratos from the command line
Let's jump in and try the code:

Use the Kratos command prompt from your Kratos installation, navigate to the folder where your script is located and execute:

```
python MainKratos.py
```

If everything went well, you should have two folders with different output formats: The GiD post file ends with `.post.bin` and can be drag and droppen into GiD. Additionally VTK files are written to the VTK_Output folder. 

## 3. The project parameters file
While you have used the `MainKratos.py` script to indicate some basic options for the simulation to work, the settings for a Kratos are stored in a `.json` file. 

JSON is an open-standard format that uses human-readable text to transmit data objects consisting of attribute–value pairs. Kratos uses a thin wrapper arround this syntax, the `Parameters` object. 

This section is a short version of a more detailed [description about the JSON syntax](How-to-write-a-JSON-configuration-file) and a [tutorial on how to read and use it](https://github.com/KratosMultiphysics/Kratos/wiki/Python-Script-Tutorial:-Reading-ProjectParameters).

### 3.1 ProjectParameters.json
The project parameters file for Kratos is commonly named `ProjectParameters.json`. Let's look at the content of this file for our structural analysis example. It contains four main blocks:

* `problem_data`: General settings for the Kratos run
* `solver_settings`: Settings for the solvers, like analysis type, linear solver, etc.
* `processes`: Processes to e.g. apply boundary conditions.  
* `output_processes`: Settings for the output

Try to change the end time of the structural case from to `5.0` seconds and run the analysis again.

## 4. The Model and ModePart files
In KratosMultiphyscis, the information about your mesh is stored in a data structure named `Model`. The model is responsible of everything related with the geometrical part of Kratos. Only one can exist per simulation and will typically contain several `ModelParts`. The serialization of this modelparts. is what we call an `.mdpa`(**M**o**d**el**Pa**rt) and is your third input file.

It contains blocks for properties, nodes, elements, conditions and initial values. In addition the mesh entities can be grouped into sub model parts. A detailed description of the syntax is given [here](Input-data).

Don't worry to much about this right now as we will dip deeper into this file and how to read it. For now just asume that the `StructuralMechanicsAnalysis` is able to read it with the information in the `ProjectParameters.json`:

```json
"model_import_settings" : {
    "input_type"     : "mdpa",
    "input_filename" : "KratosWorkshop2019_high_rise_building_CSM"
},
```


### 4.1 Read a .mdpa file
The following exercise is a short version of a more detailed [tutorial](Python-Script-Tutorial:-Reading-ModelPart-From-Input-File).

A `ModelPart` has to be created via a `Model` object, which can contain several `ModelParts`. Right after creating the empty `ModelPart`, the variables needed for the following calculations have to be added. The empty `ModelPart` is then filled using the `ModelPartIO` that reads the information from an .mdpa file. In general, the analysis object takes care of these steps, especially because it knows which variables to add.

Generaly, you will not have to deal with the lecture of a ModelPart, but 

Here you will do it directly in your python script. 
If the .mdpa file contains application dependent elements, the corresponding Kratos application has to be imported. In our structural example, the elements are from the `StructuralMechanicsApplication`. 

Extend the small python script from the first part of the exercise with the following lines:

```python
import KratosMultiphysics.StructuralMechanicsApplication

this_model = KratosMultiphysics.Model()
this_model_part = this_model.CreateModelPart("MyModelPart")

# Adding variables BEFORE reading the .mdpa
this_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

model_part_io = KratosMultiphysics.ModelPartIO("KratosWorkshop2019_high_rise_building_CSM") #path to file without ".mdpa"
model_part_io.ReadModelPart(this_model_part)
```
You could also use the filename that you extracted from the ProjectParameters.json previously.

_Hint: You can_ ```>>> print(this_model_part)``` _to see its content._ 

## Output
An extensive example on writing GiD output can be found [here](Python-Script-Tutorial:-Writing-Output-File). In this part of the tutorial you will create a minimal configuration of a VTK output process. 

The `vtk_output` block in the ProjectParameters.json gives you an impression on the potential settings for the output. Here you will create just a minimal version of it.
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
