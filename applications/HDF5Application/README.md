# HDF5Application

The *HDF5Application* enables the serialization of a model part with or without MPI using the [HDF5 library](https://support.hdfgroup.org/HDF5/). The model part is stored in and HDF5 file, which can be used for:

* Viewing a model part with [HDFVIEW](https://support.hdfgroup.org/products/java/hdfview/)
* Scientific visualization with tools supporting [XDMF](http://www.xdmf.org/index.php/Main_Page). Tools, which are known to work, include [ParaView 5.4](https://www.paraview.org/) and [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/).
* Checkpointing (under development).
* Re-partitioning (not supported yet).

- [HDF5Application](#hdf5application)
  - [Installing HDF5 (minimum version 1.8)](#installing-hdf5-minimum-version-18)
    - [Serial HDF5 Library](#serial-hdf5-library)
      - [Ubuntu](#ubuntu)
    - [Parallel HDF5 Library](#parallel-hdf5-library)
      - [Ubuntu](#ubuntu-1)
  - [Build Instructions](#build-instructions)
  - [Installing h5py](#installing-h5py)
  - [Note](#note)
  - [Kratos processes](#kratos-processes)
    - [Initialization from hdf5 process](#initialization-from-hdf5-process)
    - [Single mesh temporal input process](#single-mesh-temporal-input-process)
    - [Single mesh temporal output process](#single-mesh-temporal-output-process)
    - [Multiple mesh temporal output process](#multiple-mesh-temporal-output-process)
    - [Single mesh xdmf output process for Paraview](#single-mesh-xdmf-output-process-for-paraview)
    - [User defined I/O process](#user-defined-io-process)

## Installing HDF5 (minimum version 1.8)
The HDF5 C libraries are used with the *HDF5Application*. If Kratos is configured with MPI then the parallel HDF5 library must be installed. Otherwise, the serial HDF5 library is used.

### Serial HDF5 Library
#### Ubuntu
```
    sudo apt-get install libhdf5-dev
```
### Parallel HDF5 Library
#### Ubuntu
```
    sudo apt-get install libhdf5-openmpi-dev
```

## Build Instructions

1. Install serial or parallel HDF5 library

2. Configure Kratos to build the *HDF5Application* (following the standard [instructions](https://github.com/KratosMultiphysics/Kratos/blob/master/INSTALL.md))

    - For GNU/Linux:
    ```
    add_app ${KRATOS_APP_DIR}/HDF5Application
    ```
    - For Windows:
    ```
    CALL :add_app %KRATOS_APP_DIR%\HDF5Application;
    ```

3. Build Kratos

## Installing h5py
This package is needed to use some of the python-files, e.g. `xdmf_utils.py`
```
    sudo apt-get install python-h5py / python3-h5py
```

## Note
The minimum version for the GCC compiler is **4.9**. This is because earlier version don't fully support *regular expressions*.

## Kratos processes
There are few available HDF5 processes which can be integrated into work flow via ProjectParameters.json file.

### Initialization from hdf5 process
This process can be used to initialize a given model part using existing HDF5 files. Illustrated example reads in VELOCITY and PRESSURE variables from
"hdf5_output/example_initialization_file.h5" file to nodal historical data value container, and nodal, elemental and condition non-historical data value containers. This process also can read FLAGS from HDF5 file and populate model part accordingly.

```json
            {
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "python_module": "initialization_from_hdf5_process",
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "file_settings": {
                        "file_name": "hdf5_output/example_initialization_file.h5",
                        "echo_level": 1
                    },
                    "nodal_solution_step_data_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "element_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "condition_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_flag_value_settings": {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "element_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "condition_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    }
                }
            }
```

### Single mesh temporal input process

This process is used to initialize model part variables at each time step from hdf5 files. You can specify the list of hdf5 files to read from by changing the "file_name" with tags ```<time>``` as shown in following example. If no ```<time>``` tag is present in "file_name", then same file will be read at each time step and the variables specified will be initialized. In this process, mesh is only read from the initial time step file. Results for each container will be read from time series h5 files.

```json

            {
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "python_module": "single_mesh_temporal_input_process",
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "file_settings": {
                        "file_name": "hdf5_output/<model_part_name>-<time>.h5",
                        "time_format": "0.4f",
                        "echo_level": 1
                    },
                    "nodal_solution_step_data_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "element_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "condition_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_flag_value_settings": {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "element_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "condition_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    }
                }
            }
```

### Single mesh temporal output process
This process can be used to output data from model parts to HDF5. This will write mesh only in the first time step h5 file. Rest of the time series h5 files will contain only results for each nodal/element/condition containers.

```json
            {
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "python_module": "single_mesh_temporal_output_process",
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "file_settings": {
                        "file_name": "hdf5_output/<model_part_name>-<time>.h5",
                        "time_format": "0.4f",
                        "max_files_to_keep": "unlimited",
                        "echo_level": 1
                    },
                    "output_time_settings": {
                        "step_frequency": 1,
                        "time_frequency": 1.0
                    },
                    "nodal_solution_step_data_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "element_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "condition_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_flag_value_settings": {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "element_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "condition_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    }
                }
            }
```

### Multiple mesh temporal output process
This process is used to output model part variable data to HDF5 with mesh written for each time step. This is useful in the case if required to write down deformed mesh in each time step.

```json
            {
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "python_module": "multiple_mesh_temporal_output_process",
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "file_settings": {
                        "file_name": "hdf5_output/<model_part_name>-<time>.h5",
                        "time_format": "0.4f",
                        "max_files_to_keep": "unlimited",
                        "echo_level": 1
                    },
                    "output_time_settings": {
                        "step_frequency": 1,
                        "time_frequency": 1.0
                    },
                    "model_part_output_settings": {
                        "prefix": "/ModelData"
                    },
                    "nodal_solution_step_data_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "element_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "condition_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_flag_value_settings": {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "element_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "condition_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    }
                }
            }
```

### Single mesh xdmf output process for Paraview
This process outputs model part variable data for each time step, additionally it writes down the XDMF file after each time step which is required to visualize HDF5 data in paraview. This process requires ```h5py``` to be installed in the ```python``` version (which is compiled with ```Kratos Multiphysics```)

```json
            {
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "python_module": "single_mesh_xdmf_output_process",
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "file_settings": {
                        "file_name": "hdf5_output/<model_part_name>-<time>.h5",
                        "time_format": "0.4f",
                        "max_files_to_keep": "unlimited",
                        "echo_level": 1
                    },
                    "output_time_settings": {
                        "step_frequency": 1,
                        "time_frequency": 1.0
                    },
                    "model_part_output_settings": {
                        "prefix": "/ModelData"
                    },
                    "nodal_solution_step_data_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_data_value_settings": {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "element_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "condition_data_value_settings" : {
                        "list_of_variables": [
                            "VELOCITY",
                            "PRESSURE"
                        ]
                    },
                    "nodal_flag_value_settings": {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "element_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    },
                    "condition_flag_value_settings" : {
                        "list_of_variables": [
                            "SLIP"
                        ]
                    }
                }
            }
```

### User defined I/O process
Users can build their own I/O processes by using components of HDF5 application. Following is such example. This example writes model part, and ```DISPLACEMENT``` values for each time step at ```finalize_solution_step```

```json
            {
                "kratos_module": "KratosMultiphysics.HDF5Application",
                "python_module": "user_defined_io_process",
                "Parameters": {
                    "model_part_name": "MainModelPart",
                    "process_step": "finalize_solution_step",
                    "controller_settings": {
                        "controller_type": "temporal_controller",
                        "time_frequency": 0.5
                    },
                    "io_settings": {
                        "file_name": "results/<model_part_name>-<time>.h5"
                    },
                    "list_of_operations": [
                        {
                            "operation_type": "model_part_output"
                        },
                        {
                            "operation_type": "nodal_solution_step_data_output",
                            "list_of_variables": ["DISPLACEMENT"]
                        }
                    ]
                }
            }
```
