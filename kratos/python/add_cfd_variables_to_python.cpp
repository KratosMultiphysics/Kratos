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
//                   Riccardo Rossi
//



// System includes

// External includes


// Project includes
#include "includes/define_python.h"
#include "includes/cfd_variables.h"
#include "python/add_cfd_variables_to_python.h"

namespace Kratos
{
//KRATOS_CREATE_FLAG(STRUCTURE,   63);

namespace Python
{
    namespace py = pybind11;

    void  AddCFDVariablesToPython(pybind11::module& m)
    {

        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ADVPROJ );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONV_PROJ );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PRESS_PROJ );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MATERIAL_ACCELERATION );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VORTICITY );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RELAXED_ACCELERATION );

        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DIVPROJ );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_OLD_IT );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, C_SMAGORINSKY );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CFL_NUMBER );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MOLECULAR_VISCOSITY );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENT_VISCOSITY );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, Y_WALL);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_COEFFICIENT);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, OSS_SWITCH );

        // Legacy variables
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DYNAMIC_TAU );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DYNAMIC_VISCOSITY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EFFECTIVE_VISCOSITY );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, KINEMATIC_VISCOSITY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THAWONE );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THAWTWO );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, M );

        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CROSS_WIND_STABILIZATION_FACTOR );

    }
}  // namespace Python.
} // Namespace Kratos

