---
title: Calculate Nodal Area
keywords: process core
tags: [calculate nodal area process]
sidebar: kratos_core_processes
summary: The Calculate Nodal Area process in Kratos calculates the area associated with each node in a mesh, supporting both historical and non-historical data. This documentation outlines its usage, including descriptions of its functionality, parameters, and defaults.
---

# Calculate Nodal Area

## Description

The Calculate Nodal Area process is designed to compute the nodal area for each node within a specified model part in Kratos. This process is crucial for various simulations where nodal areas impact the calculation of forces, material properties, or other phenomena. It supports calculations for both historical (time-dependent) and non-historical data, making it versatile for different simulation needs. The process iterates through all elements in the model part, calculates the area contributed by each element to its nodes, and sums up these contributions to find the total area associated with each node.

This process is an essential tool for calculating nodal areas, facilitating various types of analyses where such information is critical.

## `CalculateNodalAreaProcess` constructors

The process first checks the provided parameters against the default parameters and validates them to ensure they are correctly specified.

It is important to ensure that the `DOMAIN_SIZE` is correctly defined in the `ProcessInfo` of the model part or explicitly provided in the parameters to avoid runtime errors.

## Execution

### `Execute`

When executed proceeds to calculate the nodal areas by integrating over the elements' geometries using their shape functions and Jacobian determinants. The calculated nodal areas are stored as nodal variables, either as historical or non-historical data depending on the template parameter `THistorical`.

This function first checks the availability of required variables and initializes nodal area values to zero. It then iterates over all elements in the model part, calculates the area contribution from each element to its nodes, and aggregates these contributions to compute the total area associated with each node.

Finally, it synchronizes the nodal area data across different processors if running in parallel.

Additionally, the process requires that the `NODAL_AREA` variable is available for the nodes in the model part; otherwise, it will terminate with an error.

## Parameters & Defaults

The process is initialized with a set of parameters specified in JSON format. The available parameters and their default values are as follows:

```json
{
    "model_part_name" : "PLEASE_DEFINE_A_MODEL_PART_NAME",
    "domain_size"     : 0
}
```

##### `model_part_name`
A string specifying the name of the model part for which nodal areas will be calculated. Default is `PLEASE_DEFINE_A_MODEL_PART_NAME`, which is a placeholder indicating that the user must specify the model part name.

##### `domain_size`
An integer defining the spatial dimension of the domain (2 for two-dimensional domains, 3 for three-dimensional domains). Default is `0`, which indicates that the domain size should be inferred from the model part's properties.