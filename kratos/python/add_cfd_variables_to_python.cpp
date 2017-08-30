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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/data_value_container.h"
#include "containers/flags.h"
#include "includes/cfd_variables.h"
#include "python/add_cfd_variables_to_python.h"

namespace Kratos
{
//KRATOS_CREATE_FLAG(STRUCTURE,   63);

namespace Python
{
    using namespace boost::python;

    void  AddCFDVariablesToPython()
    {

        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( CONV_PROJ );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PRESS_PROJ );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ACCELERATION );
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VORTICITY );

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DIVPROJ );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE_OLD_IT );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( C_SMAGORINSKY );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CFL_NUMBER );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MOLECULAR_VISCOSITY );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( TURBULENT_VISCOSITY );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( Y_WALL);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE_COEFFICIENT);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( OSS_SWITCH );

        // Legacy variables
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DYNAMIC_TAU );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DYNAMIC_VISCOSITY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( EFFECTIVE_VISCOSITY );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( KINEMATIC_VISCOSITY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( THAWONE );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( THAWTWO );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( M );

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CROSS_WIND_STABILIZATION_FACTOR );

    }
}  // namespace Python.
} // Namespace Kratos

