//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "metis_application.h"
#include "custom_python/add_processes_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;



PYBIND11_MODULE(KratosMetisApplication,m)
{

    class_<KratosMetisApplication,
           KratosMetisApplication::Pointer,
           KratosApplication >(m,"KratosMetisApplication")
           .def(init<>())
           ;
    AddProcessesToPython(m);

    //registering variables in python
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(VAUX)

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
