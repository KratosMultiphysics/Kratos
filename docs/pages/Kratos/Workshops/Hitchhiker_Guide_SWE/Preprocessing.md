---
title: Preprocessing
keywords: 
tags: [Preprocessing.md]
sidebar: kratos_workshops
summary: 
---
# Preprocessing
This section of the guide outlines the most important procedures necessary to prepare for running the simulation. The preprocessing stage involves the modelling of the  case study using GiD, followed by defining specific processes required to extract important outputs from the simulation. This guide serves to ensure proper configuration and preparation prior to the simulation.

An example file is available in Moodle, which can be downloaded under chapter *5-Computational Wind Engineering (CWE) → Example simulation of a Highrise with simulation results*. This download link contains the folder *CFD_HighRiseExampleFine*, which shows a prototype of a CFD simulation done in GiD + Kratos. It is very helpful to inspect this file and notice how the processes are implemented there.

___
## 1. Checklist

### 1.1. Preparatory Work
Before modelling in GiD, it is necessary to gather essential information regarding the case study, which will be discussed below.

#### 1. Wind profile, direction, etc.
The wind characteristics differ between each case study. Knowing the geographical position of the case study and its surroundings, you can manage to assess the main wind direction considered for the study and its respective velocity profile.

#### 2. Geometry of the structure 
It is essential to inspect the geometry of the building for the case study, which will then be modelled in GiD. Necessary information, such as height, width, depth, openings or changes of the building geometry along the height should be noted, in order to have an overview of how to model the building later.

#### 3. Sizing of the numerical wind tunnel (CFD domain size)
Make sure that your CFD domain size is correct. *CFDGuidelines* from chapter *5- Computational Wind Engineering* shows how to define the size of the domain depending on the case study. Also, have a plan for the refinement boxes with different mesh sizes. The parts close to the structure should have a finer mesh, in comparison to parts away from the structure, where the information received from the simulation is not of such importance.

___
### 1.2. GiD Model

#### 1. Prepare the geometrical model
Firstly, make sure to place your structure base on (0.0, 0.0, 0.0) and use x – streamwise, z – height, y –spanwise as directions. Afterwards, prepare necessary geometric entities, such as (in hierarchical order) nodes, lines, surfaces, volumes. Working on layers is very helpful. Prepare 2-3 bounding boxes for better mesh sizing control and refinement. 

Also, some useful troubleshooting commands are the as follows:
- **Menu bar &rarr; Utilities → Repair model**. This command checks and fixes modelling errors.
- **GiD → Geometry → Edit → Collapse → Model**. This deletes unnecessary (and perhaps hidden) lines and nodes. It is helpful especially when importing geometries from programms such as Autocad, Rhino and Sketchup, by exploding or ungrouping the geometry. 

#### 2. Prepare necessary groups 
The respective geometric entities should be put into groups and assigned the respective boundary conditions.

Make sure to have the following groups:

| Groups | Parts |
| ------ | ----- |
| Structure | All the surfaces of the structure. The bottom part of the structure should be hollow, meaning there is no surface for it, only the sides and the top. |
| Bottom Wall | Bottom surface of the CFD domain. Since the bottom of the structure is hollow, the bottom wall is a rectangular surface with an opening where the structure is positioned. |
| Top Wall | Top surface of the domain. |
| Side Walls | Both side surfaces of the domain. |
| Inlet | Surface, in which the wind is incoming (windward). |
| Outlet | Surface opposite to the inlet (leeward). |
| Fluid | All volumes of fluid. |

Boundary Conditions:

| Boundary Conditions | Groups |
| ------------------- | ------ |
| Slip boundary condition | Top Wall & Side Walls |
| NoSlip boundary condition | Bottom Wall & Structure |
| Prescribed velocity field | Inlet |
| Zero pressure field | Outlet |

#### 3. Simulation parameters
After the geometry is modelled and the groups with their respective boundary condition are defined, a proper choice of simulation parameters is important. The **CFD_HighRiseExampleFine**  file is very helpful for reference.

Make sure to have determined and assigned the following:
- Mesh sizing for structure and different boundary boxes. In the end there should be around 2.5 mio (+- 0.25 mio) elements when generating the mesh.
- Time step: Can be extracted by the CFL number.
- Simulation time: ~20T (T = maximum between estimated vortex shedding period or building fundamental period). 4T for the ramp-up phase and 16T for the simulation.
- Define the ramp-up phase function (to linearly get from 0.0 velocity to your target profile). This function will span between the time frame [0; 4T].
- Define the velocity profile function. This function will span between the time frame [4T; End(~20T)].
- Other simulation parameters can be extracted from the **CFD_HighRiseExampleFine** example.

#### 4. Run the preliminary (test) simulation 
Before running the full simulation, a good way to determine if everything has been defined correctly is to run the simulation from GiD for 3-5 time steps to make sure it runs. To see if it is running you can check GiD → Calculate → View process info. If the simulation is running for 3-5 time steps, it means that there are no errors in modelling and definition of the simulation of parameters.

___
## 2. Processes
After you are finished with the modelling of the structure in GiD and have tested that the simulation runs for 3-5 time steps, you can go to *Kratos&rarr;Write calculation files - No run*. This will create multiple file, including the *MainKratos.py* file, from which you run the simulation, and the **ProjectParameters.json** file, where you define new processes, which will be described here. The processes are necessary to define and extract information of interest from your simulation.

In order to have a copy of your original MainKratos and ProjectParameters file, we will create copies of them and name them MainKratosCustom.py and ProjectParametersCustom.json. **Make sure to update the line in MainKratosCustom referring to the project parameters from ProjectParameters &rarr; ProjectParametersCustom.**
The ProjectParametersCustom.json file will undergo a number of changes in the next steps, therefore it is advised to always refer to how the processes have been applied in the CFD_HighRiseExampleFine example. Also, discuss with your supervisors, as different projects require other process input.

The process files can be found in the the example. Do not forget to copy the files to your .gid folder.

**Note that you are not restricted only to the following processes. You are free to add more outputs of your choice.*

___
### 2.1. Point output process

The `"point_output_process"` writes results from a chosen geometrical position (under "position") in the model to a file. It first searches the entity containing the requested output location and then interpolates the requested variable(s). The output can be requested for elements, conditions and nodes (under "entity_type"). If you choose "nodes", no geometrical interpolation is performed and the exact coordinates have to be specified, therefore it is advised to specify it for "elements". Output ranges from pressure to velocities in different directions (under "output_variables").

The `"multiple_points_output_process"` serves the same purpose as the `"point_output_process"`, however, here you can define several points of interest (under `"positions"`), instead of one.

For the project work, with the structure base positioned at (0,0,0) and structural height H, we are interested in the pressure and velocities around the point (-2H,0,H). Therefore, we should add a `"point_output_process"` with these coordinates in the ProjectParametersCustom file, similar to how it is done in the example. This point is important for the Pressure Coefficient calculation.

Below is an example of the json parameters for the `"point_output_process"`:
```json
{
    "processes"        : {
        "auxiliar_process_list"            : [
            {
                "python_module" : "point_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "PointOutputProcess",
                "Parameters"    : {
                    "model_part_name"   : "FluidModelPart.fluid_computational_model_part",
                    "entity_type"       : "element",
                    "position"          : [-850.0, 0.0, 425.5],
                    "output_variables"  : ["PRESSURE", "VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z"],
                    "output_file_settings"  : {
                        "file_name"  : "reference_point_output",
                        "output_path": "results/ascii_output/ref_point_minus2H_H"
                    }
                }
            }
        ]
    }
}
```
___
### 2.2. Line output process

The `"line_output_process"` extracts output for several points along a line to a file. Internally it holds an object of type "MultiplePointsOutputProcess". By defining the "start_point" and "end_point" coordinates as well as the number of sampling points, we can receive multiple point outputs, without having to define each point separately.

For the project work, we are interested in different line outputs, such as:
- Around the contour of the building at h= 0.6H (midline) for pressure output P. Define the lines slightly outside the contour and not exactly on the surface to avoid numerical problems with location of the output points on the line. In case of a i.e. rectangular building cross sections, you will need to define 4 line output processes, with small distances ~0.1 m from the building contour.
- Pressure output in the centerline (front-top-back) in streamwise direction. In case the building is i.e. shaped as a cuboid, you will need to define 3 line output processes, with small distances ~0.1 m from the building contour.
- 3 lines, positioned respectively at -2H, -H and 0.5H, where we want to extract the velocities from the base to the height of 1.5H, in order to analyze the flow field

Below is an example of the json parameters for the `"line_output_process"`:
```json
{
    "processes"        : {
        "auxiliar_process_list"            : [
            {
                "python_module"   : "line_output_process",
                "kratos_module"   : "KratosMultiphysics",
                "process_name"    : "LineOutputProcess",
                "Parameters" : {
                    "model_part_name"   : "FluidModelPart.fluid_computational_model_part",
                    "start_point"       : [-1270.0, 0.0, 0.0],
                    "end_point"         : [-1270.0, 0.0, 650.0],
                    "sampling_points"   : 30,
                    "output_variables"  : ["PRESSURE", "VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z"],
                    "output_file_settings": {            
                        "file_name"  : "line_output_at_minus3H",
                        "output_path": "results/ascii_output/line/line_output_at_minus3H/"}
                }
            }
        ]
    }
}   
```
___
### 2.3. Force output process
After analyzing the flow field and the pressure around the building, we need to analyze the forces created from the wind loading on the building. There are 2 processes to analyze the forces on the building:

- `compute_global_force_process` is used to get the global forces (and moments) on the building in the flow- and body attached axis system. Here you need to define the reference point at the base of the building (0,0,0), in order to get the time series of the base forces and moments.
- `compute_level_force_process` is used to to get the level forces (and moment) in the flow- and body-attached axis system. We define the "start_point" and "end_point" coordinates from the building base (0,0,0) to the top of the building (0,0,H), as well as the number of intervals. The intervals are used as sampling points, meaning the number of intervals will represent the number of level forces you receive after the simulation. The forces at each level are then saved as a time history load, which will later be exported to ParOptBeam as a time history load for the computational structural dynamics (CSD), one way coupling (OWC).


Below is an example of the json parameters for the `"compute_global_force_process"` and `"compute_level_force_process"`:
```json
{
    "processes"        : {
        "auxiliar_process_list"            : [
            {
                "python_module" : "compute_global_force_process",
                "process_name"  : "ComputeGlobalForceProcess",
                "Parameters"    : {
                    "model_part_name"       : "FluidModelPart.Drag_Structure",
                    "write_output_file"     : true,
                    "print_to_screen"       : false,
                    "reference_point"       : [0.0, 0.0, 0.0], 
                    "z_rotation_angle"      : 15.0,
                    "interval"              : [0.0, "End"],
                    "output_file_settings"  : {
                        "output_path"   : "results/ascii_output/forces"
                    }
                }
            },
            {
                "python_module" : "compute_level_force_process",
                "process_name"  : "ComputeLevelForceProcess",
                "Parameters"    : {
                    "model_part_name"       : "FluidModelPart.Drag_Structure",
                    "write_output_file"     : true,
                    "print_to_screen"       : false,
                    "start_point"           : [0.0, 0.0, 0.0],
                    "end_point"             : [0.0, 0.0, 425.5],
                    "z_rotation_angle"      : 15.0,
                    "intervals"             : 20,
                    "interval"              : [0.0, "End"],
                    "output_file_settings"  : {
                        "output_path"   : "results/ascii_output/forces/level_forces"
                    }
                }
            }
        ]
    }
}   
```

___
### 2.4. CFL output process
The CFL number depends on the flow field, and differ from each case of simulation. Kratos has a module to print out the CFL number throughout the simulation. This process will give these following output for each time step:
- **max_value**: The maximum CFL number of the simulation.
- **max_id**: The ID of the element with the maximum CFL number.
- **mean**: Average CFL number in the simulation.
- **%>1.0**: Percentage of elements with CFL number above 1.0.
- **%>cfl_output_limit**: Percentage of elements with CFL number above the `"cfl_output_limit"`.

Below is an example of the json parameters for the `"cfl_output_process"`:
```json
{
    "processes"        : {
        "auxiliar_process_list"            : [
            {
                "python_module" : "cfl_output_process",
                "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
                "process_name"  : "CFLOutputProcess",
                "Parameters"    : {
                    "model_part_name"       : "FluidModelPart.fluid_computational_model_part",
                    "write_output_file"     : true,
                    "print_to_screen"       : true,
                    "cfl_output_limit"      : 2.5,
                    "interval"              : [0.0, "End"],
                    "output_step"           : 1,
                    "output_file_settings"  : {
                        "file_name"         : "cfl_results",
                        "output_path"       : "results/ascii_output/",
                        "write_buffer_size" : 1
                    }
                }  
            }
        ]
    }
} 
```

___
### 2.5. HDF5 output process
This process output the .h5 files which are necessary for visualizing the flow field in Paraview. The files contains the pressure and velocity field in the fluid domain, and the pressure distribution in the surface of the structure. 
You can change the file naming convention and its directory in the `"file_name"` field. The parameters with `<>` inside `"file_name"` indicate a placeholder value and will change depending on what is written.
- `<model_part_name>` will change depending on the `"model_part_name"` chosen in the process.
- `<time>` will change depending on the simulation time that is currently printed by the process. 
- `<step>` will change depending on the time step that is currently printed by the process. You can use this as an alternative, but using `<time>` is of course more natural.

To use the output of this process, follow the Hitchhiker guide in [postprocessing with paraview](Postprocessing.html#1-postprocessing-in-paraview){:target="_blank"}.

```json
{
    "processes"        : {
        "auxiliar_process_list"            : [
            {
                "python_module": "single_mesh_temporal_output_process",
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "Parameters": {
                    "model_part_name": "FluidModelPart.NoSlip3D_Structure",
                    "file_settings": {
                        "file_access_mode": "truncate",
                        "echo_level": 1,
                        "file_name": "results/hdf5_output/structure/<model_part_name>_T-<time>.h5",
                        "time_format": "0.2f"
                    },
                    "nodal_solution_step_data_settings": {
                        "list_of_variables": ["PRESSURE"]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": ["SCALAR_MEAN", "SCALAR_VARIANCE"]
                    },
                    "output_time_settings": {
                        "time_frequency": 0.02,
                        "step_frequency": 1
                    }
                }
            },
            {
                "python_module": "single_mesh_temporal_output_process",
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "Parameters": {
                    "model_part_name": "FluidModelPart",
                    "file_settings": {
                        "file_access_mode": "truncate",
                        "echo_level": 1,
                        "file_name": "results/hdf5_output/domain/<model_part_name>-<time>.h5",
                        "time_format": "0.2f"
                    },
                    "nodal_solution_step_data_settings": {
                        "list_of_variables": ["PRESSURE", "VELOCITY"]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": ["SCALAR_MEAN", "SCALAR_VARIANCE",
                                              "VECTOR_3D_MEAN", "VECTOR_3D_VARIANCE"]
                    },
                    "output_time_settings": {
                        "time_frequency": 0.02,
                        "step_frequency": 1
                    }
                } 
            }
        ]
    }
}   
```