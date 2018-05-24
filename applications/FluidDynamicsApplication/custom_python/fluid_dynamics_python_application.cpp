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
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "fluid_dynamics_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;


PYBIND11_MODULE(KratosFluidDynamicsApplication,m)
{

    class_<KratosFluidDynamicsApplication,
           KratosFluidDynamicsApplication::Pointer,
           KratosApplication >(m,"KratosFluidDynamicsApplication")
           .def(init<>())
           ;

    AddCustomConstitutiveLawsToPython(m);
    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PATCH_INDEX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TAUONE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TAUTWO);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PRESSURE_MASSMATRIX_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,Y_WALL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SUBSCALE_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,C_DES);
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,C_SMAGORINSKY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CHARACTERISTIC_VELOCITY);

    // For Non-Newtonian constitutive relations
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,REGULARIZATION_COEFFICIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,BINGHAM_SMOOTHER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,GEL_STRENGTH);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,SUBSCALE_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,COARSE_VELOCITY);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FIC_BETA);

    // Embedded fluid variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,EMBEDDED_IS_ACTIVE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,EMBEDDED_WET_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,EMBEDDED_WET_VELOCITY);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,Q_VALUE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,VORTICITY_MAGNITUDE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,RECOVERED_PRESSURE_GRADIENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,NODAL_WEIGHTS);

    // Compressible fluid variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,HEAT_CAPACITY_RATIO);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,REACTION_DENSITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,REACTION_ENERGY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MACH);
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
