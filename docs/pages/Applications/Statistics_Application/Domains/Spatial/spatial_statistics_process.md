---
title: Spatial Statistics Process
keywords: statistics, spatial, process
tags: [spatial_statistics_process.md]
sidebar: statistics_application
summary: 
---
This is a seperate process, which can be included in the json file as an auxiliary process. Followings are the defaults used in this process.

```json
{
    "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
    "input_variable_settings" : [
        {
            "method_name"    : "sum",
            "norm_type"      : "none",
            "container"      : "nodal_historical",
            "variable_names" : [],
            "method_settings": {}
        }
    ],
    "output_settings" : {
        "output_control_variable": "STEP",
        "output_time_interval"   : 1,
        "write_kratos_version"   : true,
        "write_time_stamp"       : true,
        "output_file_settings"   : {
            "file_name"  : "<model_part_name>_<container>_<norm_type>_<method_name>.dat",
            "output_path": "spatial_statistics_output",
            "write_buffer_size" : -1
        }
    }
}
```

`model_part_name` indicates which model part to be operated on for statistics calculation. `input_variable_settings` contains list of statistics and/or norm methods along with their acting containers of the model part. `method_name` is the statistics method name (always everything is in lower case). `norm_type` is the applied norm, `none` means, value methods are used, otherwise specified norm is used. `container` is the [container](#spatial-method-containers) to be specified for statistical value calculation. `method_settings` is used to pass additional information required by statistics method such as in [Distribution](#distribution) method. Output is written to a file in the format given at `file_name`, under the folder `folder_name`. This process can be used, if someone prefers not to modify their python scripts to calculate statistics as mentioned in the previous sections