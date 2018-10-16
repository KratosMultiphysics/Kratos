//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_DEFINE_PYTHON_H_INCLUDED )
#define  KRATOS_DEFINE_PYTHON_H_INCLUDED

/* System includes */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "includes/define.h"


#ifdef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#endif
#define KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,variable) \
	module.attr(#variable) = &variable;

#ifdef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_X) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_Y) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_Z)

#ifdef KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XX) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_YY) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XY)

#ifdef KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XX) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_YY) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_ZZ) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XY) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_YZ) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XZ)

#ifdef KRATOS_REGISTER_IN_PYTHON_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_2D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_2D_TENSOR_VARIABLE_WITH_COMPONENTS(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XX) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XY) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_YX) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_YY)

#ifdef KRATOS_REGISTER_IN_PYTHON_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_3D_TENSOR_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_3D_TENSOR_VARIABLE_WITH_COMPONENTS(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XX) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XY) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_XZ) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_YX) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_YY) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_YZ) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_ZX) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_ZY) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(module,name##_ZZ)

#ifdef KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION
#undef KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION
#endif
#define KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(module,flag) \
    module.attr(#flag) = &flag;


#ifdef KRATOS_REGISTER_IN_PYTHON_FLAG
#undef KRATOS_REGISTER_IN_PYTHON_FLAG
#endif
#define KRATOS_REGISTER_IN_PYTHON_FLAG(module,flag) \
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(module,flag);   \
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(module,NOT_##flag)

// This macro is used to print the ofstream-operator
// i.e. printing an object will give the same result in Python as in C++
// To be defined as the "__str__" function
// It replicates the function "self_ns::str(self))" of boost-python
#define KRATOS_DEF_PYTHON_STR(ObjectType) \
    [](const ObjectType& self) -> \
    const std::string { std::stringstream ss; \
    ss << self; return ss.str(); }

#endif /* KRATOS_DEFINE_H_INCLUDED  defined */
