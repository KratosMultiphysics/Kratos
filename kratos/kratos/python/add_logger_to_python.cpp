//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

// External includes
#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>

// Project includes
#include "input_output/logger.h"
using namespace boost::python;

namespace Kratos
{

namespace Python
{

/**
 * Prints the arguments from the python script using the Kratos Logger class
 * @args: tuple boost::python::object representing the arguments of the function
 * @kwargs: dictionary of boost::python::objects resenting key-value pairs for name arguments
 **/
object print(tuple args, dict kwargs) {
  std::stringstream buffer;
  LoggerOutput output(buffer);

  LoggerMessage message;

  list keys = kwargs.keys();

  // Extract the tuple part
  for(int i = 0; i < len(args); ++i) {
    object curArg = args[i];
    if(curArg) {
      message << extract<const char *>(boost::python::str(args[i])) << std::endl;
    }
  }

  // Extract the dict part
  for(int i = 0; i < len(keys); ++i) {
    object curArg = kwargs[keys[i]];
    if(curArg) {
      message << extract<const char *>(keys[i]) << ": " << extract<const char *>(boost::python::str(kwargs[keys[i]])) << std::endl;
    }
  }

  output.WriteMessage(message);

  return object();
}

void  AddLoggerToPython() {
	using namespace boost::python;

  class_<Logger, boost::shared_ptr<Logger>, boost::noncopyable>("Logger", no_init)
  .def("Print", raw_function(print,1))
  .staticmethod("Print");
}

}  // namespace Python.

} // Namespace Kratos
