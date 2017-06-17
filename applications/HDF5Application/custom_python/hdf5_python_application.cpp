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


// Project includes
#include "includes/define.h"
#include "hdf5_application.h"
#include "hdf5_application_variables.h"
#include "custom_python/add_custom_io_to_python.h"

namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosHDF5Application)
  {

	  class_<KratosHDF5Application,
			  KratosHDF5Application::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosHDF5Application")
			;

	AddCustomIOToPython();

	//registering variables in python


  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
