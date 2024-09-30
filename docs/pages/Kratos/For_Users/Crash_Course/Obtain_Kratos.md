---
title: Creating an Example
keywords: 
tags: [Crash Course Tutorial Creating example]
sidebar: kratos_for_users
summary: 
---

# Obtain Kratos

## 1. Introduction

Before starting the course, you need to acquire a working version of Kratos. Kratos can either be compiled from source or downloaded as a Python package for both Microsoft Windows and Linux, supporting Python versions `3.8` to `3.12`. 

For this course, you will use the Python package, and we will guide you through the installation process.

## 2. Install Python

If Python is not installed on your system, you can download it from the official [Python Website](https://www.python.org/) website.

## 2.1 Create a python environment:

To avoid conflicts with existing installations or packages, we strongly recommend creating a dedicated environment for Kratos.

To create a virtual environment, follow these steps:

- Create a folder for your environment.
- Run the following command:

```
python -m venv kratos_venv
```

This will create an isolated environment for Kratos. To activate or deactivate it, use the following commands:

Windows:
```ps1
# In cmd.exe
kratos_venv\bin\activate.bat

# In PowerShell
kratos_venv\bin\Activate.ps1
```

Linux:
```shell
source kratos_venv/Scripts/activate
```

If everything is set up correctly, the name of your virtual environment will appear before the console prompt.

## 3. Install Kratos

Once your virtual environment is activated, install the `KratosMultiphysics-all` package by running:

```
python -m pip install KratosMultiphysics-all
```

To verify the installation, run the following test in a Python session:

```python
import KratosMultiphysics as KMP
```

Expected output:
```console
 |  /           |
 ' /   __| _` | __|  _ \   __|
 . \  |   (   | |   (   |\__ \
_|\_\_|  \__,_|\__|\___/ ____/
           Multi-Physics 9.4.2
           Compiled for GNU/Linux and Python3.12 with GCC-15.2
Compiled with threading and MPI support.
Maximum number of threads: 24.
Running without MPI.
Process Id: 517572
```