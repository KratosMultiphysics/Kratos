---
title: Creating an Example
keywords: 
tags: [Crash Course Tutorial Creating example]
sidebar: kratos_for_users
summary: 
---

# Obtain Kratos

This is the crash course for Kratos Multiphyscis. During the course you will learn the fundamentals for out multiphyscis framework from scrach.

## 1. Preparation

KratosMultiphyscis, like may other scientific software is a set of libraries rather than a progam. You can think of it like numpy, pytorch, etc...

This means that you cannot double-click Kratos and start, but rather use its components. In order to do so, there will be a set of programs that you will need to install:

1) A Python interpreter
2) A Code Editor
3) Pre and Post processor (optional)

So the first part of this course will cover where to find these software, how to install it, and how to prepare it to run Kratos.

⚠️ If you are using MacOS, you are unable to install software in your computer, or you have any other problem, you can follow this course in a GoogleColab Notebook that you can copy from [here]().

## 2. Install Python

If Python is not installed on your system, you can download it from the official [Python Website](https://www.python.org/) website. 

Kratos is availiable for **Windows** and **Linux** from versions **3.8** to **3.12** so please select one of those.

## 2.1 Create a python virtual environment:

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

## 4. Install a Text editor

While not mandatory, is extremly recomended that you install a code editor. This is nothing more that a text editor with some tweaks to ease the production of code, like sintaxt highlight, code suggestions, etc...

We recommend you to install VSCode as is one of the most stright forward and friendly editor, but feel free to use any other editor of your choice