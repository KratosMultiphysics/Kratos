---
title: Wheel Building
keywords: Wheels
tags: [Wheels, Python]
sidebar: kratos_for_developers
summary: Introduction to python wheel building
---

## Overview

This page guides you step by step on how the wheel building can be carried out with specific python environment.

## Pre-requisites

First, you need to install the following packages for the wheel building in the python environment of which you would like to build the wheels.

```python
python -m pip install toml shutil fnmatch logging platform subprocess
```

## Wheel generation

You may need to change the following block of file `scripts/wheels/build_wheel.py`

```json
PLATFORM_CONFIG = {
    'Linux': {
        'PYTHONS': ["38", "39", "310", "311", "312", "313", "314"],
        'PYTHON_BUILD_VER': "314",
        'BUILD_SCRIPT': "scripts/wheels/linux/configure.sh",
        'BUILD_CORES': 24,
        'KRATOS_ROOT': "/workspace/kratos/Kratos",
        'WHEEL_ROOT': "/workspace/wheel",
        'WHEEL_OUT': "/data_swap_guest"
    },
    'Windows': {
        'PYTHONS': ["38", "39", "310", "311", "312", "313", "314"],
        'PYTHON_BUILD_VER': "39",
        'BUILD_SCRIPT': "scripts/wheels/windows/configure.bat",
        'BUILD_CORES': 24,
        'BASE_LD_LIBRARY_PATH': "",
        'KRATOS_ROOT': "C:/kratos/Kratos",
        'WHEEL_ROOT': "C:/dist/wheel",
        'WHEEL_OUT': "C:/data_swap_guest/"
    },
    'Darwin': {
        'PYTHONS': ["39"],
        'PYTHON_BUILD_VER': "39",
        'BUILD_SCRIPT': "scripts/wheels/darwin/configure.sh",
        'BUILD_CORES': 8,
        'KRATOS_ROOT': "/workspace/kratos/Kratos",
        'WHEEL_ROOT': "/workspace/wheel",
        'WHEEL_OUT': "/workspace/dist",
    }
}
```

. In this file, you will find the following block. If you want to make the wheels for a linux os, then make the `'PYTHONS' = [""]` (a list with one empty string item) and `'PYTHON_BUILD_VER'=""`. You can set the the `'KRATOS_ROOT'= "YOUR KRATOS PATH"`. This is the root path of the kratos, from which you would like to make the wheels. The `'WHEEL_OUT'="OUTPUT_PATH"` defines where the built wheels are stored. The path `'WHEEL_ROOT'="WHEEL_TEMP_PATH"` is used to store the temporary files generated in the wheel generation process. This folder will be automatically deleted at the end of a successful wheel generation. The `'BUILD_SCRIPT'="YOUR_BUILD_SCRIPT` defines the Kratos - Multiphysics building script for your specific wheel generation. An example for each os can be found in :

| OS      | Example build script                 |
| --------| ------------------------------------ |
| Linux   | scripts/wheels/linux/configure.sh    |
| Windows | scripts/wheels/windows/configure.bat |
| Mac OS  | scripts/wheels/darwin/configure.sh   |



An example setup for a linux system is as follows:
```python
PLATFORM_CONFIG = {
    'Linux': {
        'PYTHONS': [""],
        'PYTHON_BUILD_VER': "",
        'BUILD_SCRIPT': "scripts/wheels/linux/configure.sh",
        'BUILD_CORES': 24,
        'KRATOS_ROOT': "/software/kratos/dev",
        'WHEEL_ROOT': "/software/kratos/dev/wheels",
        'WHEEL_OUT': "/software/kratos/dev/temp"
    },
    'Windows': {
        'PYTHONS': ["38", "39", "310", "311", "312", "313", "314"],
        'PYTHON_BUILD_VER': "39",
        'BUILD_SCRIPT': "scripts/wheels/windows/configure.bat",
        'BUILD_CORES': 24,
        'BASE_LD_LIBRARY_PATH': "",
        'KRATOS_ROOT': "C:/kratos/Kratos",
        'WHEEL_ROOT': "C:/dist/wheel",
        'WHEEL_OUT': "C:/data_swap_guest/"
    },
    'Darwin': {
        'PYTHONS': ["39"],
        'PYTHON_BUILD_VER': "39",
        'BUILD_SCRIPT': "scripts/wheels/darwin/configure.sh",
        'BUILD_CORES': 8,
        'KRATOS_ROOT': "/workspace/kratos/Kratos",
        'WHEEL_ROOT': "/workspace/wheel",
        'WHEEL_OUT': "/workspace/dist",
    }
}
```