//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/kutta_condition_process.h"



namespace Kratos
{

namespace Python
{


  void  AddCustomProcessesToPython()
  {
	using namespace boost::python;

        class_<KuttaConditionProcess, bases<Process>, boost::noncopyable >("KuttaConditionProcess",init<ModelPart&>())
            .def("Execute",&KuttaConditionProcess::Execute)
            ;


  }





}  // namespace Python.

} // Namespace Kratos
