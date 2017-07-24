//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/mapping_variables.h"
#include "python/add_mapping_variables_to_python.h"

#ifdef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE
#endif
#define KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE(variable) \
    scope().attr(#variable) = boost::ref(variable);

#ifdef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE(name) \
    KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE(name##_X) \
    KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE(name##_Y) \
    KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE(name##_Z)

namespace Kratos
{

namespace Python
{
    using namespace boost::python;

    void  AddMappingVariablesToPython()
    {

    }
}  // namespace Python.
} // Namespace Kratos

#undef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_FLAG
#undef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
