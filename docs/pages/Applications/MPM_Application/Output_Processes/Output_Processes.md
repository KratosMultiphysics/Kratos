---
title: Output Processes
keywords: process core Output
tags: [Output process]
sidebar: mpm_application
summary: 
---

# Output Processes

Within the `MPMApplication`, two kinds of output processes are available:
* **grid** output processes, which are implemented in the Kratos Core (see [here](../../../Kratos/Processes/Output_Process/Output_Process)) and provide information at the background grid level;
* **material point** output processes, which are implemented within the `MPMApplication` and provide information contained in the material point elements or conditions.


## List of Grid Output Processes

- [GiD Output Process](../../../Kratos/Processes/Output_Process/GiD_Output_Process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/gid_output_process.py))
- [Vtk Output Process](../../../Kratos/Processes/Output_Process/VTK_Output_Process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/vtk_output_process.py))
- Vtu Output Process ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/vtu_output_process.py))
- [Point Output Process](../../../Kratos/Processes/Output_Process/Point_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/point_output_process.py))
- Multiple Points Output Process ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python_scripts/multiple_points_output_process.py))


## List of Material Point Output Processes

- [MPM GiD Output Process](./mpm_gid_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_gid_output_process.py))
- [MPM Vtk Output Process](./mpm_vtk_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_vtk_output_process.py))
- [MPM Point Output Process](./mpm_point_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_point_output_process.py))
- [MPM Multiple Points Output Process](./mpm_multiple_points_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_multiple_points_output_process.py))
- [MPM Json Output Process](./mpm_json_output_process) ([<i class="fa fa-github"></i> code](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MPMApplication/python_scripts/mpm_json_output_process.py))
