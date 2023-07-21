---
title: How to Compile Kratos in Debug mode
keywords: 
tags: [How-to-Compile-Kratos-in-Debug-mode.md]
sidebar: kratos_for_users
summary: 
---

## Build modes for debugging Kratos Multiphysics

Kratos Multiphysics supports different build modes that can be used to debug your code, at the cost of increased run time. The build mode is defined by the `CMAKE_BUILD_TYPE` configure variable. The following options are supported, ordered from faster (less detailed debug information) to slower (all info):

**Release** (default)

Fastest executable, no debug information included.

**RelWithDebInfo**

Compile Kratos in release mode, but keep extra information in the executable (which allows the debugger to keep track of functions and variables). Note that some parts of the code might be optimized out, and therefore not debuggable.

**Debug**

This compiles Kratos in debug mode and defines the `KRATOS_DEBUG` macro, which is used to perform additional checks during runtime.

**FullDebug** 

Compiles Kratos in debug mode and enables all available runtime checks within the code. The additional checks (compared to regular Debug mode) can be _extremely_ slow. In particular, this mode enables runtime checks for ublas matrix and vector operations, which might be helpful to detect bugs related to inconsistent matrix sizes.

Should you need more information, in Linux you can also add the compilation flag `-D_GLIBCXX_DEBUG`. This flag will define the preprocessor variable `_GLIBCXX_DEBUG`, which activates many verifications in the STL libraries. Run the code normally and it will crash when it finds some error. This also makes your code very slow. To find out the call stack, it is recommended to use a debugger, because the code crashes and tells you the type of error, but not the location in the code. 

## Selecting a build mode

To select a build type in Linux, add the following line to your configure script:

```sh
-DCMAKE_BUILD_TYPE=<build_mode>     \
```

replacing `<build_mode>` with one of `Release`, `RelWithDebInfo`, `Debug` or `FullDebug` (capitalization matters!).

In Visual Studio, you can select the build mode from the `Solution configurations` menu, as shown in the image:

![Build mode selection in Visual Studio](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/How_to_Compile_Kratos_in_Debug_mode/build_modes_win.png)

If you compile Kratos using a mode other than `Release`, you should be able to see the build mode next to the Kratos version when importing Kratos in Python:

![Running in debug mode](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/How_to_Compile_Kratos_in_Debug_mode/fulldebug_mode.png)

## Keeping both a Release and a Debug version of Kratos at the same time

You will notice that switching from one build mode to another takes long, since the entire code base has to be re-complied in the new mode. A solution to this problem is to use different build folders to compile the code in each mode. When you download kratos, it contains a `cmake_build` folder for compilation. You might want to create a different folder, for instance `cmake_build_debug`, to compile in debug mode. You can do so from the main Kratos folder: 

```sh
mkdir cmake_build_debug
cp cmake_build/configure.sh cmake_build_debug/configure_debug.sh
```

Then, **edit** 'configure_debug.sh' by setting the **build mode** with: 

```sh
-DCMAKE_BUILD_TYPE=Debug     \
```

Now you can compile the debug version from the new folder, while keeping the intermediate files of the release build intact: 

```sh
cd cmake_build_debug
sh configure_debug.sh
```

and both versions will coexist.

When you run a Python script that calls the Kratos library, it will call **the last one you compiled** or, more specifically, the one you last ran `make install` for. 
