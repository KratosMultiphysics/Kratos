---
title: Compare Two Files Check
keywords: process core
tags: [Compare Two Files Check]
sidebar: kratos_core_processes
summary: 
---

# Compare Two Files Check

## Description

This process compares files that are written during a simulation against reference files.

## Parameters & Defaults

```json
{
    "reference_file_name"   : "",
    "output_file_name"      : "",
    "remove_output_file"    : true,
    "comparison_type"       : "deterministic",
    "tolerance"             : 1e-6,
    "relative_tolerance"    : 1e-9,
    "dimension"             : 3
}
```

##### `reference_file_name` 
Reference file.

##### `output_file_name` 
Output file.

##### `remove_output_file` 
if set to `true` the output file will be removed once the process finalizes.

##### `comparison_type` 
Tye of comparion. The follwing values are accepted:

###### `deterministic`:
Checks that two files have exactly the same content.

###### `mesh_file`:
Checks that two mesh files have the same numerical values with the provided tolerance.

###### `sol_file`:
Checks that two solution files have the same numerical values with the provided tolerance.

###### `post_res_file`:
Checks that two post res files have the same numerical values.

###### `dat_file`:
Checks that two `*.dat` files have exactly the same content.

###### `csv_file`:
Checks that two `*.cvs` files have exactly the same content.

###### `mesh_file`:
Checks that two files have exactly the same content.

##### `tolerance` 
Absolute numerical tolerance used when comparing numerical values.

##### `relative_tolerance` 
Absolute numerical tolerance used when comparing numerical values.
