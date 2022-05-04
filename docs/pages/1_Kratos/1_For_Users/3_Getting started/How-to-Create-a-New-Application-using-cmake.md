---
title: How to Create a New Application using cmake
keywords: 
tags: [How-to-Create-a-New-Application-using-cmake.md]
sidebar: kratos_for_users
summary: 
---

# Overview

This method is outdated. Please go to https://github.com/KratosMultiphysics/Kratos/wiki/Creating-a-base-application for the newer version of creating a new application.

Using Cmake the creation of new application is a straightforward task. As an example we will create here an app that is empty. 
# Content
* [Downloading the 'TESTAPPLICATION.zip files to create your standard directories][link1]
* [Modifying and renaming the test_files into the desired name][link2]
    * [Understanding what is automatically performed][link3]
* [Editing other files to include the test application][link4]
    * [applications/CMakeLists.txt][link5]
    * [cmake_build/configure.sh (configure.bat in Windows)][link6]
* [Importing the new application from the scripts ( .py files)][link7]

[link1]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Create-a-New-Application-using-cmake#downloading-the-testapplicationzip-files-to-create-your-standard-directories
[link2]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Create-a-New-Application-using-cmake#modifying-and-renaming-the-test_files-into-the-desired-name
[link3]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Create-a-New-Application-using-cmake#understanding-what-is-automatically-performed
[link4]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Create-a-New-Application-using-cmake#editing-other-files-to-include-the-test-application
[link5]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Create-a-New-Application-using-cmake#applicationscmakeliststxt
[link6]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Create-a-New-Application-using-cmake#cmake_buildconfiguresh-configurebat-in-windows
[link7]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Create-a-New-Application-using-cmake#importing-the-new-application-from-the-scripts--py-files

## Downloading the 'TESTAPPLICATION.zip files to create your standard directories 

You can find a model for the structure of your application inside the [Test_application.zip](http://kratos-wiki.cimne.upc.edu/images/7/71/Test_application.zip)  file. You just have to take care to change the name of your application in the files it appears in and to start inserting your own elements, utilities etc. like in a standard, already existing application.

To begin with uncompress this folder into the applications folder. Inside of it you will find the following main parts:

* **TestApplication.py** : This file will be copied later to the KratosMultyphisics folder so that the application can be run when called in the script of your problem.
* **TestApplication.cpp** and **TestApplication.h** plus several folders including the utilities, elements, etc.
* **cMakeLists.txt** : Here you must include the source files that you want to be compiled. 

The next step is telling Kratos that we want to include and compile this application . 

## Modifying and renaming the test_files into the desired name 

We created a simple linux script to make automatic the name change (see lines below).

In order to use the script follow the following steps

1. Download and extract the TestApplication.zip file into a temporary directory (say /home/username/temp)
2. Download and extract the script [prepare_cmake.sh.zip](http://kratos-wiki.cimne.upc.edu/images/d/d5/Prepare_cmake.sh.zip)  into the same directory
3. Extract the zip files
4. Edit the file "prepare_cmake.sh" and select the name you want to provide. You need to provide lowercase,uppercase and mixed versions of the name
5. Launch the ".sh" file by typing "sh prepare_cmake.sh" at the command line (under Windows you can use the [Subsystem for Linux](https://www.microsoft.com/store/productId/9NBLGGH4MSV6))
6. Copy the modified application directory in "kratos/applications"
7. Compile and go! 

### Understanding what is automatically performed 

In order to understand what is done by the automatic script, we should consider that the following operations are needed to change the name:

If your new application name is newname you have firstly to change the name of the following files: 

```sh
test_application.h    ----> newname_application.h
test_application.cpp    ----> newname_application.cpp
test_python_application.cpp    ----> newname_python_application.cpp
```

Then you have to look inside the following files: test_application.h test_application.cpp test_python_application.cpp and perform the following substitutions: 

```sh
KRATOS_TEST_APPLICATION_H_INCLUDED    ----> KRATOS_NEWNAME_APPLICATION_H_INCLUDED
KratosTestApplication    ----> KratosNewNameApplication
include test_application.h    ----> #include newname_application.h
```

and in the CMakeLists.txt 

```sh
KRATOS_TEST_APPLICATION_SOURCES   ----> #KRATOS_NEWNAME_APPLICATION_SOURCES
test_application.cpp    ----> newname_application.cpp
test_python_application.cpp    ----> newname_python_application.cpp
TestApplication.py    ----> NewNameApplication.py
```

## Editing other files to include the test application 

### applications/CMakeLists.txt

To begin with me must include the folder we have just created so that it can be compiled: in the first line, next to the other messages we add: 

```sh
message("TEST_APPLICATION.....................${TEST_APPLICATION}")
```

and scrolling a bit below we add the new directory: 

```sh
if(${TEST_APPLICATION} MATCHES ON)
  add_subdirectory(test_application)
endif(${TEST_APPLICATION} MATCHES ON)
```

### cmake_build/configure.sh (configure.bat in Windows)

Finally you have to edit your .sh file (the one that was example_configure.sh.do_not_touch) so that your application is compiled too.

to do so, simply add the line:

```sh
-DTEST_APPLICATION=ON		\
```

among the other applications and you are ready. The new application should be compiled.

**NOTE**: It is important that the last thing on the line must be the \ . If you add spaces or other characters cmake will not be able to read the file correctly. 

## Importing the new application from the scripts ( .py files) 

To do so, if the Kratos path was set correctly, you only need to import KratosMultyphisics and then the test app: 

```sh
# Including kratos path
from KratosMultiphysics import *
from KratosMultiphysics.TestApplication import *
# And any other apps that you need to help solving your problem
```