// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "meshing_application.h"
#include "meshing_application_variables.h"
#include "custom_python/add_meshers_to_python.h"
#include "custom_python/add_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_io_to_python.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosMeshingApplication,m)
{
    py::class_<KratosMeshingApplication,
        KratosMeshingApplication::Pointer,
        KratosApplication >(m, "KratosMeshingApplication")
        .def(py::init<>())
        ;
    AddMeshersToPython(m);
    AddProcessesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomIOToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AVERAGE_NODAL_ERROR)                                  // The average nodal error
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ANISOTROPIC_RATIO)                                    // The anisotropic aspect ratio
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUXILIAR_GRADIENT)                                    // An auxiliar gradient needed to compute the metric
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUXILIAR_HESSIAN)                                     // An auxiliar hessian needed to compute the metric
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, METRIC_SCALAR);                                       // A single scalar metric
    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(m, METRIC_TENSOR_2D) // A 2D metric vector
    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(m, METRIC_TENSOR_3D) // A 3D metric vector

    //for ULF (surface_tension) application:
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,TRIPLE_POINT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CONTACT_ANGLE)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SUBSCALE_INDEX)
}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
