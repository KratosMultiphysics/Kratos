---
title: Kratos Input Files and IO
keywords: 
tags: [kratos_input_files_and_io.md]
sidebar: kratos_for_users
summary: 
---
# 1. Introduction

In the first tutorial you have successfully used Kratos to solve a few simple problems. The GiD preprocessor was used to create the input files for Kratos. There are three different types of files that compose a Kratos case:
* `.py`: A python script to run Kratos
* `.json`: Contains settings for Kratos
* `.mdpa`: Contains the model part information

These files where created by the GUI and used to directly start the simulation. In this tutorial we will have a closer look at the input files and their content. Also we will use Kratos without GUI but directly from the command line.
This is a first step in order to customize the input files for Kratos to special use cases as it will be done in the following tutorials and the more flexible usage of Kratos beyond the GUI.

**If you did not follow the [first tutorial](https://github.com/KratosMultiphysics/Kratos/wiki/Running-an-example-from-GiD), you can download the input files for a simple structural mechanics case [here](https://github.com/KratosMultiphysics/Documentation/tree/master/Workshops_files/Kratos_Workshop_2019/Sources/2_Kratos_input_files_and_IO).**

# 2. The Kratos python script
The main file of a Kratos simulation is a python script. It is responsible to load the required Kratos applications and to call the main Kratos functionalities as desired by the user.

## 2.1 MainKratos.py
The main python script for Kratos is commonly named `MainKratos.py`. Let's look at the content of this file for our structural analysis example:

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
In the first lines, Kratos and the structural analysis are imported. Then the settings are read from the `.json` file and used to create an object of the structural analysis. In the last line, the structural simulation executed.

## 2.2 Run Kratos from the command line
Use the Kratos command prompt from your Kratos installation, navigate to the folder where your script is located and execute:
```
kratos MainKratos.py
```
_Pro Tip: If you built Kratos yourself and set the paths properly as explained in the **Building Kratos** section of the Wiki, you can directly use python to execute your script._

The output of this analysis is written in two formats. The GiD post file ends with `.post.bin` and can be drag and droppen into GiD. Additionally VTK files are written to the VTK_Output folder.

## 2.3 Exercise
In order to show that the python script for Kratos indeed is just a simple python script, create a file named e.g. `my_python_script.py`. Import the Kratos core module and additionaly write some simple python commands.
```python
import KratosMultiphysics
# Custom python code
import math
a = 3
b = 4
c = math.sqrt(3**2 + 4**2)
print("-- Custom Code: c=", c, " --")
```

If you execute this script as described above you should see the following output in the terminal:
```
 |  /           |
 ' /   __| _` | __|  _ \   __|
 . \  |   (   | |   (   |\__ \
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 7.0.0
-- Custom Code: c= 5.0  --
KRATOS TERMINATED CORRECTLY
```
The python script is a powerful tool to customize a Kratos simulation, as you will see in the next tutorials.


# 3. The project parameters file
The settings for a Kratos simulation are stored in a `.json` file. JSON is an open-standard format that uses human-readable text to transmit data objects consisting of attribute–value pairs. Kratos uses a thin wrapper arround this syntax, the `Parameters` object. This section is a short version of a more detailed [description about the JSON syntax](https://github.com/KratosMultiphysics/Kratos/wiki/How-to-write-a-JSON-configuration-file) and a [tutorial on how to read and use it](https://github.com/KratosMultiphysics/Kratos/wiki/Python-Script-Tutorial:-Reading-ProjectParameters).

## 3.1 ProjectParameters.json
The project parameters file for Kratos is commonly named `ProjectParameters.json`. Let's look at the content of this file for our structural analysis example. It contains four main blocks:
* `problem_data`: General settings for the Kratos run
* `solver_settings`: Settings for the solvers, like analysis type, linear solver, etc.
* `processes`: Processes to e.g. apply boundary conditions.
* `output_processes`: Settings for the output

Try to change the end time of the structural case from to `5.0` seconds and run the analysis again.

## 3.2 Read the parameters from the .json file
Extend your custom python script by parsing the `.json` file into a `Parameters` object:
```python
with open("ProjectParameters.json",'r') as parameter_file:
    parameters = KratosMultiphysics.Parameters(parameter_file.read())
```
Extracting fields from the `Parameters` object works similar to a python dictionary. To get and set the value of a field, you have to use specific functions for the data type. In order to get the file name (a string) of the `.mdpa` file you need to type:

```python
model_part_file_name = parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
```

# 4. The Model part file
The `.mdpa`(**M**o**d**el**Pa**rt) file contains the model information in Kratos specific syntax. It contains blocks for properties, nodes, elements and conditions and initial values. In addition the mesh entities can be grouped into sub model parts. A detailed description of the syntax is given [here](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data).

## 4.1 Read a .mdpa file
The following exercise is a short version of a more detailed [tutorial](https://github.com/KratosMultiphysics/Kratos/wiki/Python-Script-Tutorial:-Reading-ModelPart-From-Input-File).

A `ModelPart` has to be created via a `Model` object, which can contain several `ModelParts`. Right after creating the empty `ModelPart`, the variables needed for the following calculations have to be added. The empty `ModelPart` is then filled using the `ModelPartIO` that reads the information from an .mdpa file. In general, the analysis object takes care of these steps, especially because it knows which variables to add.

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

# 5. Output
An extensive example on writing GiD output can be found [here](https://github.com/KratosMultiphysics/Kratos/wiki/Python-Script-Tutorial:-Writing-Output-File). In this part of the tutorial you will create a minimal configuration of a VTK output process.

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
