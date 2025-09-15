//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Hofer, https://github.com/PhiHo-eng
//                   Erich Wehrle, https://github.com/e-dub
//  based on original file from
//                   Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
#if defined(KRATOS_PYTHON)
#include "includes/define_python.h"
#include "topology_optimization_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{

namespace Python
{
    PYBIND11_MODULE(KratosTopologyOptimizationApplication,m)
    {
    namespace py = pybind11;

        py::class_<KratosTopologyOptimizationApplication,
                            KratosTopologyOptimizationApplication::Pointer,
                            KratosApplication >(m, "KratosTopologyOptimizationApplication")
                            .def(py::init<>())
                            ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YOUNGS_MODULUS_MIN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YOUNGS_MODULUS_0 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PENAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAT_INTERP )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, X_PHYS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, X_PHYS_OLD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DCDX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DVDX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SOLID_VOID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_STRAIN_ENERGY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_ELEMENT_SIZE )
    }

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
