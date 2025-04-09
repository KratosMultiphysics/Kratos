---
title: Overview
keywords: mpm
tags: [mpm]
sidebar: mpm_application
summary: 
---

Each class of material point elements and conditions defines a set of variables required for solving the given problem.
The name and type of these variables depend on what the specific material point element or condition simulates.

* [List of material point element variables](./element_variables)
* [List of material point condition variables](./condition_variables)

The values of these variables can be accessed and modified in both C++ and Python. Below, we show how to interact with them using Python.

## Get Values

Material point element and condition variables can be retrieved using the class method `CalculateOnIntegrationPoints`.
The first argument is the variable to be accessed, whereas the second is the model part `process_info`.

```python
elem_coords = mpm_element.CalculateOnIntegrationPoints(
    KratosMultiphysics.MPMApplication.MP_COORD, mpm_model_part.process_info)[0]
cond_coords = mpm_condition.CalculateOnIntegrationPoints(
    KratosMultiphysics.MPMApplication.MPC_COORD, mpm_model_part.process_info)[0]
```

{% include important.html content='Note that after calling the method we use `[0]` because the returned value is always a list' %}

## Set Values

Similarly, the values of material point element and condition variables can be set by using the class method
`SetValuesOnIntegrationPoints`.
The first argument is the variable whose value is to be assigned, the second is the actual value *always
contained in a list* and the third is the model part `process_info`.

```python
mpm_element.SetValuesOnIntegrationPoints(
    KratosMultiphysics.MPMApplication.MP_DISPLACEMENT, [[0.1, 2.2, 3.0]], mpm_model_part.process_info)
mpm_condition.SetValuesOnIntegrationPoints(
    KratosMultiphysics.MPMApplication.MPC_DISPLACEMENT, [[0.1, 2.2, 3.0]], mpm_model_part.process_info)
```
