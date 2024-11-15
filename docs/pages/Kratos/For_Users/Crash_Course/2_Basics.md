---
title: 2 - Kraots Basics
keywords: 
tags: [Kratos Crash Course Basics]
sidebar: kratos_for_users
summary: 
---

## 1. Introduction

In the initial tutorial, you successfully installed Kratos and confirmed its basic functionality. To proceed with a complete simulation, you will need several essential files that define the parameters and settings for your specific problem. This section provides an overview of these required files, detailing their roles within the simulation workflow and introducing you to fundamental Kratos structures.

The Essential Files for a Kratos Simulation

- **Main Simulation Script `main.py`**:
This Python script serves as the entry point for your simulation. It typically initializes the Kratos environment, loads problem-specific parameters, and defines the execution workflow.

- **Project Parameters File `ProjectParameters.json`**:
This JSON file contains the key configuration details for the simulation, such as solver settings, output specifications, and general options for the model part (e.g., processes, input and output options and physical properties). It centralizes the settings required by Kratos solvers and other modules.

- **Materials File `Materials.json`**
This optional JSON file defines material properties associated with different parts of the model. It is necessary when simulating multiple materials or when complex property definitions are required, ensuring accurate behavior modeling across various physical domains.

- **Model Part Definition File `model_part.mdpa`**
The model part file defines the geometry, mesh, boundary conditions, and initial conditions for the simulation. It can be created or modified using Kratos pre-processing tools or external mesh generation software (Flowgraph, Salome, GiD,...). This file format (.mdpa) is specific to Kratos and provides the structural foundation for the simulation domain.

By organizing these files in a structured manner, you can define, configure, and execute a full Kratos simulation. Each of these files plays a critical role in building, solving, and post-processing your multiphysics problem. In the following sections, we’ll cover these files in detail, providing examples and explaining key configurations required for typical use cases.

These files are usually created by a GUI (GiD, Salome, Flowgraph) and used to directly start the simulation. As the aim of the course is not to learn how to use these tools, you can download the input files for a simple structural mechanics case [here](https://github.com/KratosMultiphysics/Documentation/raw/refs/heads/master/Crash_Course/2_Basics.zip).

### 1.2. Run Kratos from the command line
Let's jump in and try the code:

To execute your simulation script in Kratos, follow these steps:

- 1) Open the Kratos Command Prompt
- 2) Navigate to the Script Directory

```bash
cd path/to/your/script/directory
```
{: data-lang="Bash"}

- 3) Run the Simulation Script

```bash
python MainKratos.py
```
{: data-lang="Bash"}

If everything went well, you should have two folders with different output formats: The GiD post file ends with `.post.bin` and can be drag and droppen into GiD. Additionally `VTK` files are written to the `VTK_Output` folder. 

## 2. The project parameters file
While you have used the `MainKratos.py` script to indicate some basic options for the simulation to work, the settings for a Kratos are stored in a `.json` file. 

`JSON` is an open-standard format that uses human-readable text to transmit data objects consisting of attribute–value pairs. Kratos uses a thin wrapper arround this syntax, the `Parameters` object. 

This section is a short version of a more detailed [description about the JSON syntax](./Input_Output_and_Visualization/JSON_Configuration_File.md) and a [tutorial on how to read and use it](./Input_Output_and_Visualization/Project_Parameters.md).

### 2.1 ProjectParameters.json
The project parameters file for Kratos is commonly named `ProjectParameters.json`. Let's look at the content of this file for our structural analysis example. It contains four main blocks:

* `problem_data`: General settings for the Kratos run
* `solver_settings`: Settings for the solvers, like analysis type, linear solver, etc.
* `processes`: Processes to e.g. apply boundary conditions.  
* `output_processes`: Settings for the output

Try to change the end time of the structural case from to `5.0` seconds and run the analysis again.

## 3. The Model and ModePart files
In KratosMultiphyscis, the information about your mesh is stored in a data structure named `Model`. The model is responsible of everything related with the geometrical part of Kratos. Only one can exist per simulation and will typically contain several `ModelParts`. The serialization of this modelparts. is what we call an `.mdpa`(**M**o**d**el**Pa**rt) and is your third input file.

It contains blocks for properties, nodes, elements, conditions and initial values. In addition the mesh entities can be grouped into sub model parts. A detailed description of the syntax is given [here](./Input_Output_and_Visualization/Input_Data.md).

Don't worry to much about this right now as we will dip deeper into this file and how to read it. For now just asume that the `StructuralMechanicsAnalysis` is able to read it with the information in the `ProjectParameters.json`:

```json
"model_import_settings" : {
    "input_type"     : "mdpa",
    "input_filename" : "KratosWorkshop2019_high_rise_building_CSM"
},
```
{: data-lang="JSON"}

## 4. The Kratos python script
### 4.1. MainKratos.py
The last, but the most important file that you have downloade is called `MainKratos.py` and as its name suggest is a python script. This would be the equivalent to your well known `.exe` for a classic problem, with the advantage that you can customize it according to your needs. It is responsible to load the required Kratos applications and to call the main Kratos functionalities as desired by the user. 

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
{: data-lang="Python"}

In the first lines, Kratos and the structural analysis are imported. Then the settings are read from the `.json` and a Model is created. Finaly we use all that information to create a `StructuralMechanicsAnalysis` simulation. In the last line, the structural simulation executed. 

### 4.2. Analysis Stage
As you've observed, most of our work so far has involved reading configuration files and initializing a few objects to start the simulation. This simplicity is due to the encapsulation of simulation steps within the `AnalysisStage` class. For our purposes, we are using the `StructuralMechanicsAnalysis` class, which manages the main simulation loop for structural mechanics applications in Kratos.

Each Kratos application typically has a specialized `AnalysisStage` tailored to its requirements, ensuring that the necessary processes and solvers are set up correctly for the specific type of simulation (e.g., structural mechanics, fluid dynamics).

The foundational structure and logic of the AnalysisStage can be reviewed in the  [analysis_stage.py](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/analysis_stage.py) file in the Kratos repository. This file contains the base AnalysisStage class, providing the core logic and function definitions.

Again, don't worry to much for now, as we will cover this class in detailes later. For now Let's take a quick look and let's create a custom `AnalysisStage` that we will use during the coruse:

```python
import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class CourseAnalysisStage(StructuralMechanicsAnalysis):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        print("Custom Analysis Stage Created!")

    def Run(self):
        print("running Custom Analysis Stage!")
        super().Run()

if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = CourseAnalysisStage(model, parameters)
    simulation.Run()
```
{: data-lang="Python"}

As you can see, we have made a derived class from the `StructuralMechanicsAnalysis` that we were running which has two methods changed:
- `__init__`: This will be called every time a instance is created and we will greet us with the custom message we have added:
- `Run`: This is the method that we were calling in the `MainKratos.py` and we have modified it to also print an info message.

Be mindful that we are now creating an analysis stage of our custom `CourseAnalysisStage`:

```python
simulation = CourseAnalysisStage(model, parameters)
simulation.Run()
```
{: data-lang="Python"}

## 5. Wrap up
With this overview, you should now be familliar with the most Basic files (MainKratos.py, ProjectParameters.json, geometry.mdpa) and data strcutures (AnalysisStage, Model, Parameters) of Kratos.

In the following sections we will take a closer look at each one of those and use them to introduce more detailed concepts.