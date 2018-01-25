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
#include "includes/define.h"
#include "co_simulation_application.h"
#include "co_simulation_application_variables.h"
#include "custom_python/add_base_classes_to_python.h"
#include "custom_python/add_custom_application_interfaces_to_python.h"
#include "custom_python/add_custom_coupling_strategies_to_python.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_relaxation_schemes_to_python.h"
//#include "custom_python/"


namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosCoSimulationApplication)
  {

	  class_<KratosCoSimulationApplication,
			  KratosCoSimulationApplication::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosCoSimulationApplication")
			;

    AddCustomBaseClassesToPython();
    AddCustomApplicationIterfacesToPython();            
	AddCustomCouplingStrategiesToPython();
    AddCustomIoToPython();
    AddCustomRelaxationSchemesToPython();
    
	//registering variables in python

  }

  // We can add a second module for stand alone users.

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
