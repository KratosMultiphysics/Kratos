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

object printImpl(tuple args, dict kwargs, Logger::Severity severity) {
    std::cout << "?" << std::endl;
    if(len(args) == 0)
        std::cout << "ERROR" << std::endl;

    Logger& loggerRef = boost::python::extract<Logger&>(args[0]);
    
    std::stringstream buffer;

    // Extract the tuple part
    std::cout << len(args)-1 << std::endl;
    for(int i = 1; i < len(args); ++i) {
        object curArg = args[i];
        std::cout << "Printing arg " << i << std::endl;
        std::cout << extract<const char *>(boost::python::str(args[i])) << std::endl;
        if(curArg) {
            buffer << extract<const char *>(boost::python::str(args[i])) << std::endl;
        }
    }

    loggerRef << KRATOS_CODE_LOCATION << buffer.str() << severity;
    std::cout << "--?" << std::endl;

    return object();
}

/**
 * Prints the arguments from the python script using the Kratos Logger class
 * @args: tuple boost::python::object representing the arguments of the function The first argument is the label
 * @kwargs: dictionary of boost::python::objects resenting key-value pairs for
 * name arguments
 **/
object printInfo(tuple args, dict kwargs) {
    return printImpl(args, kwargs, Logger::Severity::INFO);
}

/**
 * Prints the arguments from the python script using the Kratos Logger class
 * @args: tuple boost::python::object representing the arguments of the function The first argument is the label
 * @kwargs: dictionary of boost::python::objects resenting key-value pairs for
 * name arguments
 **/
object printWarning(tuple args, dict kwargs) {
    return printImpl(args, kwargs, Logger::Severity::WARNING);
}

void  AddLoggerToPython() {
	using namespace boost::python;

    class_<Logger, boost::shared_ptr<Logger>, boost::noncopyable>("Logger", init<std::string const &>())
    .def("PrintInfo", raw_function(printInfo,1))
    .def("PrintWarning", raw_function(printWarning,1))
    ;
}

}  // namespace Python.

} // Namespace Kratos
