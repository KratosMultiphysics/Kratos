---
title: Temporal Statistics Process
keywords: statistics, temporal, process
tags: [temporal_statistics_process.md]
sidebar: statistics_application
summary: 
---
This is a seperate process, which can be included in the json file as an auxiliary process. Followings are the defaults used in this process.

```json
    "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
    "input_variable_settings" : [
        {
            "method_name"     : "sum",
            "norm_type"       : "none",
            "container"       : "nodal_historical_non_historical",
            "echo_level"      : 0,
            "method_settings" : {}
        }
    ],
    "echo_level" : 0,
    "statistics_start_point_control_variable_name" : "TIME",
    "statistics_start_point_control_value"         : 0.0
```

`model_part_name` is the model part which all of the mentioned statistics will be calculated. `input_variable_settings` contains list of statistical methods, their norms and container combinations where statistics will be calculated. under each `input_variable_settings_block`, `method_name` is the statistics method, `norm_type` is the norm type (`"none"` for value methods, other wise norms are used), `container` to tell from which container to read/write, `method_settings` to specifiy input/output variables for the method. Apart from that, there are few global settings. `statistics_start_point_control_variable_name` is the temporal statistics calculation start control variable name. This variable should be present in the model_part specified in order to calculate temporal statistics. Statistics calculation will start when this variable's value passes user specified value in `statistics_start_point_control_value`.