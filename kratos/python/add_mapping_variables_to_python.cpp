//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/mapping_variables.h"
#include "python/add_mapping_variables_to_python.h"

namespace Kratos
{

namespace Python
{
    using namespace boost::python;

    void  AddMappingVariablesToPython()
    {
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(TANGENT_FACTOR)
    }
}  // namespace Python.
} // Namespace Kratos

