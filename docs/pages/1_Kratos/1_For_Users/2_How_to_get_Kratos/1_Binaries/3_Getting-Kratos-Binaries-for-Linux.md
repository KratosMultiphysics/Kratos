---
title: Getting Kratos Binaries for Linux
keywords: 
tags: [Getting-Kratos-Binaries-for-Linux.md]
sidebar: kratos_for_users
summary: 
---

This page covers the process of downloading a compiled version of Kratos for Linux and how to use it to run an example. 

## Minimum requirements

* Ubuntu 14.04 or greater (other flavors of Linux are supported but not tested. Ie. Suse, CentOS, RedHat...)

## Download and Installation

You can find the Kratos binaries for Linux in the [Release Section](https://github.com/KratosMultiphysics/Kratos/releases/tag/7.0)

1 - Download the tgz file for linux

2 - Extract the contents of the tgz in a location of your choice

3 - kratos.sh is the script you will use to launch problems.

This finalizes the installation process.

## Usage

1 - Download any example of your interest from the [Examples Repo](https://github.com/KratosMultiphysics/Examples) and execute:

```
kratos.sh MainKratos.py
```

## Advanced

If you prefer to skip the installation part, or have multiple versions of kratos installed simultaneously:

1 - Copy this `runner.sh` script and edit the contents of the `KRATOS_ROOT` variable:

```Bash
KRATOS_ROOT=/path/to/kratos   # Please change this variable

export PYTHONHOME=${KRATOS_ROOT}
export PYTHONPATH=${KRATOS_ROOT}
export LD_LIBRARY_PATH=${KRATOS_ROOT}/libs:${KRATOS_ROOT}/OpenMPI/lib:/home/roigcarlo/KratosInstall/libs
${KRATOS_ROOT}/runkratos $*
```

Do not forget to give the script execution rights:
```bash
chmod +x runner.sh
```

2 - Make an alias for the script

```bash
echo "alias kratos=runner.sh" >> $HOME/.bashrc
```

3 - Open the terminal and execute the problem script using kratos:

```bash
kratos MainKratos.py
```