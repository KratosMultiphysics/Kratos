---
title: Run optimization
keywords: 
tags: [Run_optimization.md]
sidebar: shape_optimization_application
summary: 
---

This section selects how to run the defined optimization problem.

## Run optimization dialog

<p align="center">
    <img src="images/run_optimization.png" alt="Run optimization dialog box"/>
</p>
<p align="center">Figure 1: Run optimization dialog box</p>

|Option|Description|
|------|-----------|
|Submit| This descides the optimization problem submission type. **Write & Optimize**: This will write a new optimization files to be read by KratosMultiphysics, and then starts the optimization procedure.**Write Input File**: This will only write the input files which are used by KratosMultiphysics. **Optimize existing Input File**: This will not write a new input files for KratosMultiphysics, instead it will use the existing input files to start the otpimization procedure. |