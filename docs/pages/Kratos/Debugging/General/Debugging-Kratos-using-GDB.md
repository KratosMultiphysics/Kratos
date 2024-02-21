---
title: Debugging Kratos using GDB
keywords: 
tags: [Debugging-Kratos-using-GDB.md]
sidebar: kratos_debugging
summary: 
---

# Debugging Kratos using GDB
This page details the process of debugging Kratos with GDB. The guide is intended for Linux systems using a terminal, if you prefer to debug on Windows or through gdb using externals tools like visual studio code we have some pages detailing those processes as well (see https://github.com/KratosMultiphysics/Kratos/wiki/How-to-cross-debug-Kratos-under-Windows or https://github.com/KratosMultiphysics/Kratos/wiki/Debugging-Kratos-using-Visual-Studio)

## Required software
* GDB
* Python / Python-dbg

To correctly compile a version of Kratos for debug you must first chose one of the following compilation modes, which will result in different degrees of reported information and execution speed:

* **RelWithDebInfo**: Optimization flags will be still active when possible and debug symbols are enabled. This is by far the fastest option but since optimizations are enabled the code being reported may or may not be an accurate representation of the one you have in the source files. This also could (and will) lead to potentially "optimized values" when inspecting variables. 
* **Debug**: This is the standard Debug mode. Optimizations are disabled and debug symbols enabled. Code execution will be slow but no variables or pieces of code should change, meaning that you can inspect all parts of the code.
* **FullDebug**: This option is similar as the option above but also introduces some flags which enable debug only checks which will make the code significantly slower and introduce an extra layer of security. This option is general good if you are not sure where your bug occurs.

As python is used as the main for the Kratos executions, some times you will want to know exactly which are the entry points that python used to execute C++ code. While this information is usually not available from GDB itself, you can obtain it if you execute te code through a python binary compiled with debug symbols. In order to do so, you must first install the python-dbg package:

```console
sudo apt-get install python3-dbg
```

And then modify you `configure.sh` script so you target the debug Python version:

```console
# Set basic configuration
export KRATOS_BUILD_TYPE="Debug"
export PYTHON_EXECUTABLE="/usr/bin/python3-dbg"
```

Then it is required to run the `configure.sh` script again in order to compile the Kratos Python export with the debug version of Python. We note that this might take a while.

## Debugging the code

Once the code is compiled, and the env paths set as in any other compilation you may start debugging your code with GDB by calling it directly with the python command and passing your script either as argument, or from GDB itself by doing

```
gdb python3-dbg
run MainKratos.py
```

Is not the scope of this guide to detail how to use GDB (you can find a complete manual here: https://www.gnu.org/software/gdb/documentation/) but some basic commands that will get you going are:

* **b [file]**: Adds a breakpoint in the specified line of a file. Is it common that you receive the message "code not found, do you want to set the breakpoints when the code is loaded?". This is expected when you try to set a breakpoint before Kratos libs are loaded.

* **d [num]**: This deletes a given breakpoint

* **s**: Once the code stops, lets to step into the next instruction. This is specially useful for breakpoints in functions as it lets you dive in the function that is being called.

* **n**: Similar to `s` but it will execute the next line in the current frame, executing full functions calls without stooping (unless there is another breakpoint inside)

* **c**: Continue the execution until the next breakpoint.

* **l**: Shows you the current line of code being executed and the context around it. 

* **bt**: Shows the current stack call. This information will be missing at python level unless compiled with python-dbg

* **up**: Allows you to move one frame up in the stack (going to the function above you)

* **down**: Allows you to move one frame down in the stack (going to the function below you)

* **p [symbol]**: prints the value of any symbol or expression in the current frame.

## MPI

_**Note**: This section is very introductory. We are aware that there are more sophisticated methods in order to perform some of the actions described here, but this should serve as start point to begin to work with._

It is also possible to debug MPI code using GDB. Although there are different approaches here (e.g. attaching the debugger to individual processes), we usually obtain bests results by launching different GDB instances and executing each one of the processes in those instances of the debugger.
This processes can be laborious if you want to debug even a small number of processes. Hence, we provide the next launcher in order to automatize such process:

*mpirun_debug.sh*:
```console
#!/bin/bash

BREAKPOINTS=`echo "${GDB_BREAKPOINTS}" | sed -r -e "s/([^;]++);?/-ex 'break \1' /g"`
echo "Using breakpoints: ${BREAKPOINTS}"

if ! [ -n "${MPI_RUNRUN_DBG+1}" ]; then
  MPI_RUNRUN_DBG="gdb -ex 'set breakpoint pending on' ${BREAKPOINTS} -ex run --args"
fi

if ! [ -n "${MPI_RUNRUN_TERM+1}" ]; then
  MPI_RUNRUN_TERM="xterm -hold -e"
fi

eval "mpirun -np $1 ${MPI_RUNRUN_TERM} ${MPI_RUNRUN_DBG} ${@:2}"
```

This script requires execution permissions. This can be done in the command line by doing
```
chmod +x mpirun_debug.sh
```

Then, this script will automatically set up N instances of gdb, launch the command specified and run the code with a single instruction:

```console
alias mpirund=~/.scripts/mpirun_debug.sh
```

```console
mpirund NUM_PROCS python3-dbg MainKratos.py
```

Be aware that in order to correctly debug the code, most of the times you will want to set the same breakpoints in all processes. You can do this by using `GDB_BREAKPOINTS` environment variable as you would within GDB and the breakpoints there will be automatically loaded to all processes open:

```console
export GDB_BREAKPOINTS="my_file1.cpp:20;my_file2.h:13140"
```