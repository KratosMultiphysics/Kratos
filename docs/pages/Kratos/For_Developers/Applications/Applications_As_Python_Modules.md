---
title: Applications as Python Modules
keywords: 
tags: [Applications Python Modules Legacy]
sidebar: kratos_for_developers
summary: 
---

https://github.com/KratosMultiphysics/Kratos/pull/3217/ enables the applications can be used as python-modules.

This means that the python-scripts in the Applications can now be used e.g. as:
```python
# old:
import KratosMultiphysics
import analysis_stage # this works because the folder "python_scripts" is added to the python_path

# new, pythonic way
from KratosMultiphysics import analysis_stage
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic
```

The new way is the standard already for new applications created with the Application-Generator. For already existing applications currently both ways are supported, but it is recommended to the developers to use the new way (pythonic way) of importing python-scripts because at some point the old import-mechanism (which adds the python_scripts to the python-path) will be removed eventually, and then all python-files should be using the new way already.

# The following three steps are recommended for developers:
## 1. Modify the `CMakeLists.txt` in your Application (replacing `DummyApplication` with the real name of your Application):

Change the location where the `DummyApplication.py` is being installed to: 
- from: `install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DummyApplication.py" DESTINATION KratosMultiphysics )` 
- to `install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DummyApplication.py" DESTINATION "KratosMultiphysics/DummyApplication" RENAME "__init__.py")`

## 2. Update `DummyApplication.py` to use the new import mechanism

The `DummyApplication.py` should look like this:

```python
# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosDummyApplication import *
application = KratosDummyApplication()
application_name = "KratosDummyApplication"
application_folder = "DummyApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
```

## 3. Update all python-files to use the new/pythonic import-mechanism
As shown above

## The changes result in the following structure in the **KratosMultiphysics**-Folder:

previously:

```
KratosMultiphysics
 |-- __init__.py
 |-- application_importer.py    
 |-- kratos_globals.py
 |-- ... (other kratos files)
 |-- DummyApplication.py
```

new:

```
KratosMultiphysics
 |-- __init__.py
 |-- application_importer.py    
 |-- kratos_globals.py
 |-- ... (other kratos files)
 |-- DummyApplication
      |-- __init__.py
```