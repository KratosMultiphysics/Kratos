---
title: Overview
keywords: process mpm output
tags: [output mpm process]
sidebar: mpm_application
summary: 
---

# Output Processes

In `KratosMultiphysics`, simulation results are typically written to a file (or multiple files) using Python classes that extend the base `OutputProcess` class, which is implemented in the `KratosCore`.
More details about the `OutputProcess` class can be found [here](../Output_Processes/Overview). The parameters required by these output processes must be defined in the `"output_processes": {}` section of the `ProjectParameters.json` file (more details about this input file can be found [here](../Input_Files/json#projectparametersjson)).


In the `MPMApplication`, two kinds of output processes can be used, i.e., the **grid** and the **material point** output processes.

## Material Point Output Processes
The **material point** output processes, which are implemented within the `MPMApplication`, provide information contained in the material point elements or conditions.
This type of output processes includes:
- [MPM GiD Output Process](./mpm_gid_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_gid_output_process.py))
- [MPM Vtk Output Process](./mpm_vtk_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_vtk_output_process.py))
- [MPM Json Output Process](./mpm_json_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_json_output_process.py))
- [MPM Write Energy Output Process](./mpm_write_energy_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_write_energy_output_process.py))

## Grid Output Processes
The **grid** output processes, which are implemented in the Kratos Core (see [here](../../../Kratos/Processes/Output_Process/Output_Process)), provide information at the background grid level.
This type of output processes includes:
- [GiD Output Process](../../../Kratos/Processes/Output_Process/GiD_Output_Process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/gid_output_process.py))
- [Vtk Output Process](../../../Kratos/Processes/Output_Process/VTK_Output_Process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/vtk_output_process.py))
- Vtu Output Process ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/vtu_output_process.py))
- [Point Output Process](../../../Kratos/Processes/Output_Process/Point_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/point_output_process.py))
- Multiple Points Output Process ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/multiple_points_output_process.py))

## Restart Output Process

Within the `MPMApplication`, the [`SaveRestartProcess`](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/save_restart_process.py) implemented in the Kratos Core can be used for the periodic saving of all data contained in the computing (MPM) model part. This functionality acts as a checkpoint system, allowing users to resume a simulation from the last saved state using the [`RestartUtility`](https://github.com/KratosMultiphysics/Kratos/blob/master/python_scripts/restart_utility.py) and without having to restart from the beginning.

- [Save & Restart a Simulation](./restart)
