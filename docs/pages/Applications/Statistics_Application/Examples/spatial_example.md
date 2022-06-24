---
title: Spatial Process Example
keywords: statistics, spatial, process, example
tags: [spatial_example.md]
sidebar: statistics_application
summary: 
---
Following example illustrates different methods used in different containers with different norms. `input_variable_settings` holds an array of methods for specified containers, specified norm and specified variables. They can be customized for your requirement. `output_settings` holds information about how the output should be handled.

```json
{
                "kratos_module" : "KratosMultiphysics.StatisticsApplication",
                "python_module" : "spatial_statistics_process",
                "Parameters" : {
                    "model_part_name" : "test_model_part",
                    "input_variable_settings" : [
                        {
                            "method_name"    : "sum",
                            "norm_type"      : "none",
                            "container"      : "nodal_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "mean",
                            "norm_type"      : "none",
                            "container"      : "nodal_non_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "variance",
                            "norm_type"      : "none",
                            "container"      : "element_non_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "rootmeansquare",
                            "norm_type"      : "none",
                            "container"      : "condition_non_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "sum",
                            "norm_type"      : "magnitude",
                            "container"      : "nodal_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "mean",
                            "norm_type"      : "pnorm_2.5",
                            "container"      : "nodal_non_historical",
                            "variable_names" : ["VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "variance",
                            "norm_type"      : "component_x",
                            "container"      : "condition_non_historical",
                            "variable_names" : ["VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "rootmeansquare",
                            "norm_type"      : "index_3",
                            "container"      : "nodal_non_historical",
                            "variable_names" : ["LOAD_MESHES"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "min",
                            "norm_type"      : "frobenius",
                            "container"      : "nodal_non_historical",
                            "variable_names" : ["GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        }
                    ],
                    "output_settings" : {
                        "output_control_variable": "STEP",
                        "output_time_interval"   : 1,
                        "write_kratos_version"   : false,
                        "write_time_stamp"       : false,
                        "output_file_settings"   : {
                            "file_name"  : "<model_part_name>_<container>_<norm_type>_<method_name>.dat",
                            "output_path": "spatial_statistics_process",
                            "write_buffer_size" : -1
                        }
                    }
                }
            }
```