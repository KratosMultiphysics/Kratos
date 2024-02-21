---
title: Editing the main files
keywords: 
tags: [Tutorial-Editing-the-main-files.md]
sidebar: kratos_for_developers
summary: 
---

## Overview

In the main files of our application: my_laplacian_application.h and my_laplacian_application.cpp we will have to include the code describing our element and condition and tell _Kratos_ how to use it. The original files created in the previous section require some minor modifications, which are explained below: 

## Header file: my_laplacian_application.h

The first lines are just the typical legal issues related to the Kratos license. Just add your name as the main author of this application.

```cpp
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//
```

In order to avoid class redefinition, you have to define a label for your application definition, as follows: 

```cpp
#if !defined(KRATOS_MY_LAPLACIAN_APPLICATION_H_INCLUDED )
#define KRATOS_MY_LAPLACIAN_APPLICATION_H_INCLUDED
```

the `#endif` is at the end of the code (the rest of the lines must be included between these ones): 

```cpp
#endif // KRATOS_MY_LAPLACIAN_APPLICATION_H_INCLUDED defined
```

We will need to add the "include" files, such as the element and condition:

```cpp
#include "custom_elements/my_laplacian_element.h"
#include "custom_conditions/point_source_condition.h"
```

The class `MyLaplacianApplication` definition comes inside the _Kratos_ namespace. Then, you will need to declare the element and conditions as private member variables:

```cpp
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const MyLaplacianElement mMyLaplacianElement;
    const PointSourceCondition mPointSourceCondition;

    ///@}
```

As specified in the [Style Guide](Style-Guide#class-members), all member variables inside a class are identified by a lower case _m_ at the beginning of the name.

## Source file: my_laplacian_application.cpp 

First of all, remember to customize the header and add the include files. The elements initialization will need the geometries headers:

```cpp
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    HERE_YOUR_NAME
//

// System includes

// External includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/line_2d_2.h"
#include "geometries/point_2d.h"
#include "my_laplacian_application.h"
#include "my_laplacian_application_variables.h"
```

The implementation file includes the initialization of the element and conditions inside the constructor, and then, registering them:

```cpp
namespace Kratos {

KratosMyLaplacianApplication::KratosMyLaplacianApplication():
    KratosApplication("MyLaplacianApplication"),
    mMyLaplacianElement( 0, Element::GeometryType::Pointer( new Triangle2D3<Node>( Element::GeometryType::PointsArrayType (3) ) ) ),
    mPointSourceCondition( 0, Element::GeometryType::Pointer( new Point2D  <Node>( Element::GeometryType::PointsArrayType (1) ) ) )
    {}

void KratosMyLaplacianApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosMyLaplacianApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE( POINT_HEAT_SOURCE )

    KRATOS_REGISTER_ELEMENT( "MyLaplacianElement", mMyLaplacianElement )

    KRATOS_REGISTER_CONDITION( "PointSourceCondition", mPointSourceCondition )

}
}  // namespace Kratos.
```

## Variables files
We only need to create `POINT_HEAT_SOURCE` variable since the other variables needed in this app are already included in the kernel  (`CONDUCTIVITY` and `TEMPERATURE`). The application generator has carried out the next tasks:

* **Define** the variables in the `application_variables.h` (just remember to customize the header)
* **Create** the variables in the `application_variables.cpp`
* **Register** the variables in **Kratos** in the `application.h`, as we have seen above
* **Register** the variables in **Python** in the `custom_python/python_application.cpp` file.


## CMakeLists.txt

Finally, we will need to add the new sources to the cmake file (only the .cpp are needed).


This is currently done automatically with the commands. Check that you have these lines in your CMakeLists.txt. 

```txt
## MyLaplacian Core sources
file(GLOB_RECURSE KRATOS_MY_LAPLACIAN_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/my_laplacian_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/my_laplacian_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)
```
