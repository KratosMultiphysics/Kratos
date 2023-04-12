---
title: ParOptBeam Guide
keywords: 
tags: [ParOptBeam_Guide.md]
sidebar: kratos_workshops
summary: 
---
# ParOptBeam Guide
ParOptBeam (Parametrical Optimizable Beam) is a useful program from the Chair of Structural Analysis, in which the user can run eigenvalue, static and dynamic analysis of a finite element (FE) beam model. After the simulation is finished, we receive time series of level forces, which we [convert to the ParOptBeam format](Postprocessing.html#3-converting-level-forces-to-paroptbeam){:target="_blank"}. For the project work, we will create an equivalent beam model of the building, and run a dynamic analysis using the level forces from the CFD simulation in Kratos.

This part of the guide will help you with setting up the model and running the dynamic analysis of the building as well as postprocessing the output .

___
## 1. ParOptBeam Master Structure
You can [download ParOptBeam under this link](https://github.com/mpentek/ParOptBeam){:target="_blank"}. Click on *Code &rarr; Download zip* to download the master file of ParOptBeam. The structure of the master folder you download is explained in this chapter.

### 1.1. Input
For the analysis to run, the program requires user input. The input required from the user are the **ProjectParameters** and the **Forces**.

- **ProjectParameters**: Under **input &rarr; parameters** you can find three different files with the project parameters. You will also have to create one file for your own building, which will be explained in [chapter 2](#2-setting-up-the-model-project-parameters) of this guide. The parameters consist of the Material, Geometry, Boundary Conditions, Optimization Parameters and types of Analysis (together with the respective output options of analysis).

- **Force**: Under **input &rarr; force &rarr; generic_building/generic_pylon** you can find three different files with with the dynamic forces in an .npy format. For your project you will also have to import your dynamic forces, which will be produced from the [convert_kratos_to_paroptbeam](Postprocessing.html#3-converting-level-forces-to-paroptbeam){:target="_blank"} script. Depending on how many level forces you have, it is advised that dynamic force should be named after the number (i.e. 15 nodes &rarr; dynamic_force_15 nodes). When modelling your beam model, **make sure that the beam has the same amount of nodes as the dynamic force!**

### 1.2. Source
The source folder is interesting for you to see and read, in order to understand the processes going on, however probably not necessary to work with during your project. It contains the different types of analysis definitions, beam elements, structure model and postprocessing as well as preprocessing settings.

### 1.3. Run File
To run ParOptBeam, you will need to run a python file, which calls the input (ProjectParameters and Dynamic Forces) and then runs the respective analysis based on these parameters. *run_generic_models* i.e. runs the simulation for all three ProjectParameters in the master ('TestParametersStraightBeam.json','ProjectParameters3DGenericBuilding.json', ProjectParameters3DGenericPylon.json'). In the project work, you will need to [create your own ProjectParameters file](#2-setting-up-the-model-project-parameters) and run it specifically. Another example of how a run file can look like is *run_generic_models_from_python*, where the project parameters are directly implemented in the run file instead of called from the input.

### 1.4. Output
If you run the *run_generic_models*, from the defined parameters, an eigenvalue, static, as well as dynamic analysis will run. You will notice that a new output folder will be created, displaying the results from the respective analysis (for the *GenericBuilding*, *GenericPylon* as well as *TestStraightBeam*). The results of interest are already defined by the user in the project parameters, and then printed out as .dat files and plotted in the *analyses_results_report* as a pdf file. The beam element has 6 degrees of freedom (dof-s) to choose from: 3 translations and 3 rotations. Each analysis has separate outputs:

- **Eigenvalue Analysis**: Eigenfrequencies, Mode Shapes and Types of Mode.
- **Static Analysis**: Deformations, External Forces and Reactions.
- **Dynamic Analysis**: Deformations, Velocity, Acceleration, External Forces and Reactions, however as time series.
- **Structure Model Properties**: Shows an overview of the beam model properties.

It is advised to have a look on the different examples from the ParOptBeam Master in order to understand the processes better.

___
## 2. Setting up the Model Project Parameters
This chapter will guide you through the project parameters setup for your own building. You will need to create your own seperate ProjectParameters and run file. You can start by copying the *ProjectParameters3DGenericBuilding.json* file, renaming it, and adapting to your own project.

### 2.1. Model Parameters
Under the model parameters you will create the geometry of the building and define its parameters. Keep in mind, that all the units in ParOptBeam are SI-Units (m, s, kg, N)

- **General Parameters**:
	- `"name"`: This will then be the name for your project which will be called from the run file.
	- `"domain_size"`: Should be 3D.
	- `"system_parameters"`:
		- `"type"`: Pick between "Timoshenko", "Bernoulli" and "CRBeam" (co-rotational beam), depending on how you want to model the elements of the beam. 
		- `"is_nonlinear"`: Pick false (linear). Nonlinear is still not established.

    ```json
    {
        "model_parameters": {
            "name": "GenericBuilding",
            "domain_size": "3D",
            "system_parameters": {
                "element_params": {
                    "type": "Timoshenko",
                    "is_nonlinear": false
                }
            }
        }
    }
    ```

- **Material Parameters**:
    - **is_nonlinear**: Pick false (linear). Nonlinear is still not established.
    - **Density [kg/m^3]**: Make an assumption upon the number of floors in your building, the area of each slab and also consider a certain factor for other structural elements, such as walls, columns, facade etc. 
    - **Young's Modulus [N/m^2]**: Assumption, based on the building material (concrete, steel, other?).
    - **Possion Ratio [-]**: Assumption, based on the building material (concrete, steel, other?).  
    - **Damping Ratio [-]**: Assumption, based on the building material (concrete, steel, other?) and structure. 

    ```json
    {
        "model_parameters": {
            "system_parameters": {
                "material": {
                    "is_nonlinear": false,
                    "density": 160.0,
                    "youngs_modulus": 2.861e8,
                    "poisson_ratio": 0.1,
                    "damping_ratio": 0.0
                }
            }
        }
    }
    ```

- **Geometrical Parameters**:
    - **Length_x [m]**: Length of the beam (=height of the building).
    - **Number of elements**: Pick number of elements based on (1) cross section of building along the height, (2) how detailed the beam model should be (i.e. 2 elements are too   little to achieve correct results) and (3) number of nodes of the beam model should be the same as the number of level forces. If you notice, that you need to model with  more elements than your current nodal forces, you can always go back to *convert_kratos_to_paroptbeam* and pick a larger number of sampling intervals.
    - **Defined on intervals**: Choose the number of intervals in which the following cross sectional parameters are defined. The number of intervals can be lower, however not larger, than the number of elements. 
    - **Interval bounds [m]**: Pick between which heights the beam segment is defined. Make sure that the first interval starts at 0.0 m and the last ends with "End", and that the end of one interval is the beginning of the interval after it.
    - **Length y & z [m]**: Lengths of cross section of the building. Make sure that it is clear, which is the y and which is the z- Direction of the building. 
    - **Area, Shear Area [m^2] & Moment of inertia y, z, torsional [m^4]**: Pick these based on the cross section parameters of the building.
    - **Outrigger**: Adds stiffness and mass effects of an outrigger. Input 0.0, in case your building has no outriggers.

    ```json
    {
        "model_parameters": {
            "system_parameters": {
                "geometry": {
                    "length_x": 180,
                    "number_of_elements": 3,
                    "defined_on_intervals": [{
                        "interval_bounds" : [0.0, 60.0],
                        "length_y": [55.0],
                        "length_z": [35.0],
                        "area"    : [1925.0],
                        "shear_area_y" : [1605.0],
                        "shear_area_z" : [1605.0],
                        "moment_of_inertia_y" : [196510.0],
                        "moment_of_inertia_z" : [458260.0],
                        "torsional_moment_of_inertia" : [691771.0],
                        "outrigger" :{
                            "mass": 947000.0,
                            "stiffness_ratio_y": 5,
                            "stiffness_ratio_z": 5
                        }
                    },{"next intervals"}
                    ]
                }
            }
        }
    }
    ```

- **Boundary conditions**: This will define how your beam model is supported. There are different options, such as:
    - fixed-fixed
    - pinned-pinned
    - fixed-pinned
    - pinned-fixed
    - fixed-free
    - free-fixed

    *In the case of a tall building, free-fixed is the correct condition (cantilever beam). Also free-fixed can be chosen, however be careful then with the choice of parameters, as they should be reversed.*

    ```json
    {
        "model_parameters": {
            "boundary_conditions": "fixed-free"
        }
    }
    ```

___
### 2.2. Optimization Parameters
Since it is difficult to adapt the beam geometrical and material parameters to a real building, an optimization process is helpful when knowing certain parameters of the building, such as the total mass and predicting its eigenfrequencies. By running the optimization, the geometrical and material parameters of the beam model are adapted to receive the target eigenfrequencies and total mass of the building, therefore receiving more realistic values for the cross section of the beam model.

The following parameters should be considered (non-mentioned parameters can be left as in the *ProjectParameters3DGenericBuilding*):
- **Density for total mass [kg]**: Input the total mass of the building. This will adapt the density of the beam model to its total mass.
- **Consider decomposed modes**: Choose which decomposed modes should be considered for the optimization. Sway z represents a rotation around the z axis (translation in y direction). For Sway y the same principle stands as for Sway z. Also, the torsional mode can be considered. 
- **Corresponding mode ids**: The corresponding mode ids for the considered decomposed modes.
- **Corresponding eigenfrequencies**: The corresponding eigenfrequencies for the considered decomposed modes.

After the optimization, the following parameters will be updated for the considered decomposed modes:

- **Sway y Optimization**: Update of Shear area z, Moment of Inertia y.
- **Sway z Optimization**: Update Shear area y, Moment of Inertia z.
- **Torsional Optimization**: Update of Torsional moment of inertia.

As well as the density for the total mass of the building. 

___
### 2.3. Analyses Parameters
As mentioned, you can run an eigenvalue, static and dynamic analysis in ParOptBeam. The following analysis parameters should be considered for the project (only most essential parameters will be explained):

#### 1.  General Parameters:
  - **Global Output Folder**: Defines the name of the output folder. It is suggested to have the same output folder name as project name, to avoid confusions.
  - **Model_Properties - Write/Plot**: This defines if the structural model properties should be printed as a .dat file and/or plotted in the .pdf report or not. 
  - **Report Options**: Settings for the .pdf output.
  - **Skin Model Parameters**: Settings for the visualization of the beam element as a skin model.

#### 2. Eigenvalue Analysis Parameters :
Here we select for which modes we want to print and/or plot the mode shapes and eigenfrequencies as well as animate, after running the eigenvalue analysis. Pick as many modes as you find reasonable and for which the information is important. For the animation, [installation of FFMPEG](Installation_Guides.html#14-ffmpeg){:target="_blank"} is required.

#### 3.  Dynamic Analysis Parameters:
 - **Time**: Pick the start, end and time step of simulation according to the CFD Simulation.
 - **File path**: Provide the load path in the required .npy format. 
 - **Plot/Write Step and/or Time**: Print or plot the results of a certain step or time.
 - **Animate**: True or false, depending if you are interested in the animation of the dynamic analysis (can be helpful to understand if the load input/direction is correct).
 - **Skin Model Animation Parameters**: Select start, end time as well as steps of the animation.
 - **Dof_List**: List of degrees of freedom considered, which are input as numbers. You can see the available dofs under "*source&rarr; global_definitions*". The following degrees of freedom are represented by these numbers (for a fixed-free boundary condition):
    - **Translation x (longitudinal)**: Dof 0
    - **Translation y**: Dof 1
    - **Translation z**: Dof 2
    - **Torsional Rotation a**: Dof 3
    - **Bending Rotation around y-axis**: Dof 4
    - **Bending Rotation around z-axis**: Dof 5

 - **Result type**: For each dof written in the dof list, you should add corresponding results of interest. It can be a "reaction", "ext_force", "displacement", "velocity" or "acceleration". You can also extract multiple results for the dof of interest, as long as they are listed inside of brackets, for example:

    ```json
    ["reaction", "displacement", "velocity"]
    ```
 - **Plot/Write result**: For the corresponding dofs and respective results of interest, you should input True/False if you want the results printed as .dat files or plotted in the report pdf.

- **Remark**:
The results of the dynamic analysis are written and plotted in the original ParOptBeam script for the base of the of the cantilever beam. This is correct when extracting reactions, however when extracting accelerations, velocities or displacements, other positions, possibly at the top of the building, are of interest.  In this case, go to *analysis &rarr; dynamic_analysis* to the definition of write_result_at_dof (same goes for plot_result_at_dof). Under selected_result = "displacement", "velocity" and "acceleration", dof represents the degree of freedom of the first node of the cantilever, so the base of the structure. Multiply the element number of interest with 6 (6 dofs per node) and add the dof, to extract the result of interest. For example, if you are interested in the acceleration at the top of the building and your model has 6 segments, write:


    result_data = self.solver.acceleration[6*6 + dof, :] 

    **instead of**

    result_data = self.solver.acceleration[dof, :] 

___
#### 4. Static Analysis Parameters:
In principle the same as the dynamic analysis parameters, however for a selected time step. In the project work, the dynamic analysis is probably of more interest, as it covers more information than the static analysis.


___
### 3. Postprocessing
After the parameters have been defined, update your run file, in order to call the ProjectParameters.json file of your building. After running the analysis, a new output folder for your project will be written, containing .dat files, animations, and an analysis pdf report.

If you chose "plot" for some of the analysis results, then the results will be plotted in the analysis report. However, it is suggested to still create some scripts (similar to [Postprocessing: Ascii output](Postprocessing.html#2-point-and-line-ascii-data){:target="_blank"}), to add additional plots/results which might be valuable for the project. Some interesting results are (pick relevant dofs!):

- Plot of time series of base reactions (with maximum value noted).
- Plot of time series of displacements at the top of the building (with maximum value noted).
- Plot of time series of accelerations at the top of the building (with maximum value noted).



