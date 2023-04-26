---
title: Temporal Process Example
keywords: statistics, temporal, process, example
tags: [temporal_example.md]
sidebar: statistics_application
summary: 
---
Following is an example on how to use `temporal_statistics_process` in json. The output variables data type should be matched with the input variables order and the data type for `"norm_type" = "none"`. If `"norm_type" != "none"`, then output variables should be `scalars`. If `"container" = "nodal_historical_historical"` is used as the container type, then output variables should be added to `NodalSolutionStepVariables` list in Kratos since this container type outputs temporal statistics variables to nodal historical container. This json settings also can be added to `auxiliary_processes` list.

For details about all the available statistical methods, norm_types, etc, please refer to rest of the `README.md` file.

```json
        {
            "kratos_module": "KratosMultiphysics.StatisticsApplication",
            "python_module": "temporal_statistics_process",
            "Parameters": {
                "model_part_name": "FluidModelPart.fluid_computational_model_part",
                "input_variable_settings": [
                    {
                        "method_name": "variance",
                        "norm_type": "none",
                        "container": "nodal_historical_non_historical",
                        "echo_level": 1,
                        "method_settings": {
                            "input_variables": [
                                "VELOCITY",
                                "PRESSURE"
                            ],
                            "output_mean_variables": [
                                "VECTOR_3D_MEAN",
                                "SCALAR_MEAN"
                            ],
                            "output_variance_variables": [
                                "VECTOR_3D_VARIANCE",
                                "SCALAR_VARIANCE"
                            ]
                        }
                    }
                ],
                "statistics_start_point_control_variable_name": "TIME",
                "statistics_start_point_control_value": 2.5
            }
        }
```