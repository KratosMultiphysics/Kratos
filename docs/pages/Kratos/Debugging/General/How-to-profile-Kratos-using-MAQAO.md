---
title: How to profile Kratos using MAQAO
keywords: 
tags: [How-to-profile-Kratos-using-MAQAO.md]
sidebar: kratos_debugging
summary: 
---

# Overview
MAQAO is a lightweight tool which is useful for analyzing the quality of binary code and identifying hotspots. It is available at [www.maqao.org](http://www.maqao.org/). MAQAO can be used for any build type, but it provides more detailed information when debugging symbols are available.

# Identifying Hotspots

## Profiling a Kratos Simulation
Set up a problem to analyze in Kratos and run

`maqao lprof -ldi=on xp=lprof_output -- /usr/bin/python3 MainKratos.py`

to profile the Kratos simulation and save the results in the local directory `lprof_output`.

For MPI run

`mpirun -np 4 maqao lprof -ldi=on xp=lprof_output -- /usr/bin/python3 MainKratos.py`

## Analyzing Results
To print a summary of results in the terminal run `maqao lprof -df xp=lprof_output`.

To generate html files run `maqao lprof -df xp=lprof_output of=html`.

# Code Quality Analysis

Code quality analysis can be applied to binary code in the Kratos libraries. Function bodies and function loops can be analyzed.

## Analyzing Function Bodies

Run `maqao cqa Kratos.so fct-body="factorize"` to analyze all functions with names containing _factorize_ in Kratos.so.

## Analyzing Function Loops

Run `maqao cqa Kratos.so fct-loops="factorize"`.
