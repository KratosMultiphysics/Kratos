---
title: Profiling Kratos with cProfile and VTune
keywords: 
tags: [Profiling Kratos with cProfile and VTune.md]
sidebar: kratos_debugging
summary: 
---

# Profiling Kratos with cProfile and VTune

## Profiling Python code with cProfile
This page details the steps to follow in order to profile Kratos. Profiling the python part can be done with cProfile. For this no modification of the code is necessary. To visualize the profilling results **SnakeViz** is recommended. SnakeViz can be installed with the following command using pip:

```console
pip install snakeviz
```
or using anaconda:
```console
conda install -c anaconda snakeviz
```
More details can be found here: https://jiffyclub.github.io/snakeviz/

In order to run the profiler type the following command:
```console
python -m cProfile -o outputFile.prof MainKratos.py
```
To view the results graphically in the browser:
```console
snakeviz outputFile.prof
```
The profiling results contain the number of calls for each function, the time per call, the total time and the cummulative time.

## Profiling with VTune

VTune can be used as a stand-alone program or as part of Intel oneAPI.
It can be downloaded from:
https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler-download.html

If further details are required a look at the documentation is recommended: https://www.intel.com/content/www/us/en/develop/documentation/vtune-help/top.html

After installation VTune can be run with the gui or the command line. For viewing the profiling results the gui is recommended. Inside the gui the following things are needed:

* Application = python executable
* Application parameters = MainKratos.py
* Working directory = directory of MainKratos.py
* Analysis type: e.g. Hotspots

After starting kratos:

```console
startkratos
```
vtune can be run with the gui or the command line.

The hotspot analysis gives the run time for each **function**. To see the runtime for each line of code kratos needs to be compiled with debug symbols:
```console
compilekratosrelwdbg
```