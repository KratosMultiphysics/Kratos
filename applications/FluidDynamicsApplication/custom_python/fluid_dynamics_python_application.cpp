//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>
#include "boost/python/detail/wrap_python.hpp"

// Project includes
#include "includes/define.h"
#include "fluid_dynamics_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosFluidDynamicsApplication)
{

    class_<KratosFluidDynamicsApplication,
           KratosFluidDynamicsApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosFluidDynamicsApplication")
           ;

    AddCustomConstitutiveLawsToPython();
    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();
    AddCustomProcessesToPython();

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PATCH_INDEX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TAUONE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(TAUTWO);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PRESSURE_MASSMATRIX_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(Y_WALL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SUBSCALE_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(C_DES);
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE(C_SMAGORINSKY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(CHARACTERISTIC_VELOCITY);

    // For Non-Newtonian constitutive relations
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(REGULARIZATION_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(BINGHAM_SMOOTHER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(GEL_STRENGTH);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(SUBSCALE_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(COARSE_VELOCITY);

    // Embedded fluid variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EMBEDDED_IS_ACTIVE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EMBEDDED_WET_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(EMBEDDED_WET_VELOCITY);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(Q_VALUE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(VORTICITY_MAGNITUDE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(RECOVERED_PRESSURE_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_WEIGHTS);
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
