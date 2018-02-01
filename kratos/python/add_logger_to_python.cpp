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
    if(len(args) == 0)
        std::cout << "ERROR" << std::endl;
    
    std::stringstream buffer;

    // Extract the tuple part
    for(int i = 1; i < len(args); ++i) {
        object curArg = args[i];
        if(curArg) {
            buffer << extract<const char *>(boost::python::str(args[i]));
        }
    }

    const char* label = extract<const char *>(boost::python::str(args[0]));
    Logger(label) << KRATOS_CODE_LOCATION << "NEWLOG: " << buffer.str() << severity << std::endl;

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
    .staticmethod("PrintInfo")
    .staticmethod("PrintWarning")
    ;

    // Enums
    enum_<Logger::Category>("Category")
    .value("STATUS", Logger::Category::STATUS)
    .value("CRITICAL", Logger::Category::CRITICAL)
    .value("STATISTICS", Logger::Category::STATISTICS)
    .value("PROFILING", Logger::Category::PROFILING)
    .value("CHECKING", Logger::Category::CHECKING);
}

}  // namespace Python.

} // Namespace Kratos
