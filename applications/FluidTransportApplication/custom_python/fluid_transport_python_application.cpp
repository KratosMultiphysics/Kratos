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
// External includes
#include <boost/python.hpp>

// Project includes
#include "fluid_transport_application.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

BOOST_PYTHON_MODULE(KratosFluidTransportApplication)
{
	class_<KratosFluidTransportApplication,
	KratosFluidTransportApplication::Pointer,
	bases<KratosApplication>, boost::noncopyable >("KratosFluidTransportApplication");

	//registering variables in python
}


}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
