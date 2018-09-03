//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes

#if defined(KRATOS_PYTHON)
// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

#include "fluid_transport_application.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosFluidTransportApplication, m)
{
	class_<KratosFluidTransportApplication,
	KratosFluidTransportApplication::Pointer,
    KratosApplication>(m, "KratosFluidTransportApplication")
	.def(init<>());

	//registering variables in python
}


}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
