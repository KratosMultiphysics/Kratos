---
title: Debugging Kratos using Visual Studio
keywords: 
tags: [Debugging-Kratos-using-Visual-Studio.md]
sidebar: kratos_debugging
summary: 
---

# Debugging Kratos using Visual Studio
As Kratos runs using Python scripts the debugging process using Visual Studio may be a bit different of what a normal user may be used to. You have two different known ways to be able to debug your code: Cross Debugging Python or debugging directly the C++ part. This page will focus in the second option. You can find a dedicated page of the wiki for cross debugging [here](How-to-cross-debug-Kratos-under-Windows).

## Required software
* Visual Studio 2015/2017 (see [Windows Install](Windows-Install)) for debugging C++
* Boost libraries in debug mode.

## Preparing the environment 

First make sure the code is configured and compiled in `Debug` or `Fulldebug` modes. We recommend you to have a separate build directories for `Release` and `Debug`.

![](https://user-images.githubusercontent.com/1935791/35916965-859b3fe8-0c0c-11e8-867a-46a0f8f62f69.png)


Go to Kratos project, on the solution explorer, right click and select __Properties__. From here go to __Debugging__ menu. you should see a window like this:

![](https://user-images.githubusercontent.com/1935791/35917027-d5f36934-0c0c-11e8-8730-943ac6aeb213.png)

In this menu, you will need to change some options

- __Command__: You should select python.exe you used to compile Kratos. Ex: _C:\python36\python.exe_
- __Command arguments__: Name of the script of the case you want to run. Ex: _MainKratos.py_
- __Working directory__: Path to the folder where the script above is located. Ex: _C:\MyExamples\MyProblem.gid_
- __Environment__: You have to a add
  - __PATH__: so it also points to the `Libs` directory of Kratos where you have compiled the code. _Ex: C:\Kratos_install\Libs_
  - __PYTHONPATH__: it has to point to your kratos folder where the KratosMultiphysics directory is located. typically the same as above but without libs. Ex: _C:\Kratos_install\Libs_

## Debugging
After that just right click in the __Kratos__ -> __Debug__ -> __Start new debug instance__:

![](https://user-images.githubusercontent.com/1935791/35917255-bec136c8-0c0d-11e8-90ce-389c384cad57.png)

Notice that you can put breakpoints and stop executions in all applications and components, not only in _Kratos_. Running from _Kratos_ project just serves as an entry point in Visual Studio.

If you have any problem please feel free to open an issue in our project's page: https://github.com/KratosMultiphysics/Kratos/issues
