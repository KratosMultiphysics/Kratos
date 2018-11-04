//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    msandre
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/mesh_moving_variables.h"
#include "python/add_mesh_moving_variables_to_python.h"

namespace Kratos
{

namespace Python
{
    namespace py = pybind11;

    void  AddALEVariablesToPython(pybind11::module& m)
    {
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MESH_DISPLACEMENT)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MESH_ACCELERATION)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MESH_REACTION)
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MESH_RHS)

        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,LAPLACIAN_DIRECTION)

    }
}  // namespace Python.
} // Namespace Kratos
