---
title: Getting Kratos Binaries for Windows
keywords: 
tags: [Getting-Kratos-Binaries-for-Windows.md]
sidebar: kratos_for_users
summary: 
---

This page covers the process of downloading a compiled version of Kratos for Windows and how to use it to run an example. 

## Minimum requirements 

* Windows 7 x64 or greater
* Advanced text editor is recommended

## Download and Installation

You can find the Kratos binaries for Windows in the [Release Section](https://github.com/KratosMultiphysics/Kratos/releases/tag/7.0).

1 - Download the zip file for windows

2 - Extract the contents of the zip in a location of your choice

This finalizes the installation process.

## Usage

1 - Download any example of your interest from the [Examples Repo](https://github.com/KratosMultiphysics/Examples)

2 - Open the `Kratos` command line.

3 - Navigate to the folder and execute the script using kratos
```cmd
kratos MainKratos.py
```

### Example Script

```python
import KratosMultiphysics as KM
import KratosMultiphysics.ExternalSolversApplication as ESA
```