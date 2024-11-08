---
title: 1 - Obtain Kratos & First Steps
keywords: 
tags: [Kratos Crash Course Download First Editor Pre Post]
sidebar: kratos_for_users
summary: 
---

# Obtain Kratos

Welcome to the Kratos Multiphysics Crash Course! In this course, you’ll gain a comprehensive understanding of the fundamentals of Kratos, a powerful and flexible multiphysics framework.

This course is designed to guide you through Kratos from the ground up, covering essential concepts, components, and workflows. By the end, you will have a solid foundation for setting up, running, and extending Kratos simulations across multiple physical domains.

## 1. Preparation

KratosMultiphysics, like many scientific software tools, is a collection of libraries rather than a standalone program. This means that Kratos doesn’t operate as a traditional application you can double-click to launch. Instead, similar to libraries like NumPy or PyTorch, Kratos provides components and functions that you integrate into your own scripts and workflows.

This means that you cannot double-click Kratos and start, but rather use its components. To effectively work with Kratos, you'll need to set up an environment that includes specific programs and tools to interface with its components and manage simulations. Below is a list of essential programs and tools for using Kratos:

- A Python interpreter
- A Code Editor
- Pre and Post processor (optional)

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
```bash
# In cmd.exe
kratos_venv\bin\activate.bat

# In PowerShell
kratos_venv\bin\Activate.ps1
```
{: data-lang="Ps1"}

Linux:
```bash
source kratos_venv/Scripts/activate
```
{: data-lang="Bash"}

If everything is set up correctly, the name of your virtual environment will appear before the console prompt.

## 3. Install Kratos

Once your virtual environment is activated, install the `KratosMultiphysics-all` package by running:

```
python -m pip install KratosMultiphysics-all
```
{: data-lang="Bash"}

⚠️ If you are a windows user, please also install [Windows Redistributable](https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170)

To verify the installation, run the following test in a Python session:

```python
import KratosMultiphysics
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
{: data-lang="Output"}

## 4. Install a Text editor

While not mandatory, is extremly recomended that you install a code editor. This is nothing more that a text editor with some tweaks to ease the production of code, like sintaxt highlight, code suggestions, etc...

We recommend you to install [VSCode](https://code.visualstudio.com/) as is one of the most stright forward and friendly editor, but feel free to use any other editor of your choice

## 5. Pre and Post processor

To fully utilize KratosMultiphysics, it is highly recommended that you install pre- and post-processing tools. These tools are essential for setting up simulation models, visualizing results, and performing detailed analysis. Kratos is compatible with several popular pre- and post-processors, each designed to address different needs and user segments like GiD, Salome and Paraview.

Since all the geometries for this course are already done, we recommend you to just install [Paraview](https://www.paraview.org/) for visualization.