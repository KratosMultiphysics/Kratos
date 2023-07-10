---
title: How to cross debug Kratos under Windows
keywords: 
tags: [How-to-cross-debug-Kratos-under-Windows.md]
sidebar: kratos_debugging
summary: 
---

# Cross debugging (Python and C++) of Kratos under Windows
The presented development environment allows to cross debug Kratos (Phyton and C++). In addition, it is possible to easily switch between the _Debug_ and _Release_ configuration. Moreover some hints how to bypass common errors in setting up Kratos under Windows (see also [Windows Install](Windows-Install)) are provided.

## Required software
* Visual Studio 2015 (see [Windows Install](Windows-Install)) for debugging C++
* Visual Code (https://code.visualstudio.com/) for debugging Python

### Visual Code
For debugging Kratos, Visual Code requires a Pyhton extension. The following screenshot shows some useful extensions.
 
![Useful Visual Code extensions](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Home/Extensions.PNG)

## Setting up Kratos
To allow for fast switching between the _Debug_ and _Release_ configuration it is strongly recommended to use two different build folders (e.g. build_release and build_debug) with two different configuration files (see -DCMAKE_BUILD_TYPE=Debug or Release in configure.bat).
### Hints
* Make sure you have compiled the boost library in the _Release_ and _Debug_ modus
* Make sure the following environment variables are set properly PYTHONPATH=[e.g. C:\Kratos_install], LD_LIBRARY_PATH=[e.g. C:\Kratos_install\libs] and add e.g. C:\Kratos_install and C:\Kratos_install\libs to the variable PATH
* Make sure you have copied the necessary _Release_ AND _Debug_ .dll files into the Kratos installation folder (see Post Compilation in [Windows Install](Windows-Install))

## Cross Debugging (Python and C++)
The first step is to perform one debugging step (Python) in Visual Code (Hint: use "stopOnEntry": true in the launch.json file). This creates a Python process.

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Home/Kratos_VisualCode_Python.PNG)

Once the first debugging step in Visual Code is done, one can attach the Visual Studio Debugger (C++) to the before created Python process.

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Home/Attach_to_process.png)

Done! Now you can cross debug Kratos (Python and C++) under Windows.

### Hints
* For switching between _Debug_ and _Release_ you need to open the two different KratosMultiphysics.sln files created by CMake in the corresponding build folders e.g. build_release and build_debug

Have fun by cross debugging Kratos under WINDOWS ;-)

Michael Breitenberger