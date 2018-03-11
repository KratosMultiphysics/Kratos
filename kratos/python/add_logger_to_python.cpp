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
#include "includes/define.h"
#include "input_output/logger.h"

using namespace boost::python;

namespace Kratos {
namespace Python {

/**
 * Prints the arguments from the python script using the Kratos Logger class. Implementation
 * @args tuple boost::python::object representing the arguments of the function The first argument is the label
 * @kwargs dictionary of boost::python::objects resenting key-value pairs for
 * @severity Logger::Severity The message level of severity @see Logger::Severity
 * @useKwargLabel bool Indicates if the label must be gather from kwargs (true) or is the first argument of the call (false)
 * name arguments
 **/
object printImpl(tuple args, dict kwargs, Logger::Severity severity, bool useKwargLabel) {
    if(len(args) == 0)
        std::cout << "ERROR" << std::endl;
    
    std::stringstream buffer;
    Logger::Severity severityOption = severity;
    Logger::Category categoryOption = Logger::Category::STATUS;

    const char* label;

    // Get the label
    if(useKwargLabel) {
        if(kwargs.contains("label")) {
            label = extract<const char *>(kwargs["label"]);
        } else {
            label = "undefined";
        }
    } else {
        label = extract<const char *>(boost::python::str(args[0]));
    }

    // Extract the tuple part
    for(int i = (useKwargLabel ? 0 : 1); i < len(args); ++i) {
        object curArg = args[i];
        if(curArg) {
            buffer << extract<const char *>(boost::python::str(args[i])) << ((i != len(args)) ? " " : "");
        }
    }

    // Extract the options
    if(kwargs.contains("severity")) {
        severityOption = extract<Logger::Severity>(kwargs["severity"]);
    }

    if(kwargs.contains("category")) {
        categoryOption = extract<Logger::Category>(kwargs["category"]);
    }

    // Send the message and options to the logger
    Logger(label) << buffer.str() << severityOption << categoryOption << std::endl;

    return object();
}

/**
 * Prints the arguments from the python script using the Kratos Logger class. Default function uses INFO severity.
 * @args tuple boost::python::object representing the arguments of the function The first argument is the label
 * @kwargs dictionary of boost::python::objects resenting key-value pairs for
 * name arguments
 **/
object printDefault(tuple args, dict kwargs) {
    return printImpl(args, kwargs, Logger::Severity::INFO, true);
}

/**
 * Prints the arguments from the python script using the Kratos Logger class using INFO severity.
 * @args tuple boost::python::object representing the arguments of the function The first argument is the label
 * @kwargs dictionary of boost::python::objects resenting key-value pairs for
 * name arguments
 **/
object printInfo(tuple args, dict kwargs) {
    return printImpl(args, kwargs, Logger::Severity::INFO, false);
}

/**
 * Prints the arguments from the python script using the Kratos Logger class using WARNING severity.
 * @args tuple boost::python::object representing the arguments of the function The first argument is the label
 * @kwargs dictionary of boost::python::objects resenting key-value pairs for
 * name arguments
 **/
object printWarning(tuple args, dict kwargs) {
    return printImpl(args, kwargs, Logger::Severity::WARNING, false);
}

void  AddLoggerToPython() {

    class_<LoggerOutput, boost::shared_ptr<LoggerOutput>, boost::noncopyable>("LoggerOutput", no_init)
    .def("SetMaxLevel", &LoggerOutput::SetMaxLevel)
    .def("GetMaxLevel", &LoggerOutput::GetMaxLevel)
    .def("SetSeverity", &LoggerOutput::SetSeverity)
    .def("GetSeverity", &LoggerOutput::GetSeverity)
    .def("SetCategory", &LoggerOutput::SetCategory)
    .def("GetCategory", &LoggerOutput::GetCategory)
    ;

    scope logger_scope = class_<Logger, boost::shared_ptr<Logger>, boost::noncopyable>("Logger", init<std::string const &>())
    .def("Print", raw_function(printDefault,1))
    .def("PrintInfo", raw_function(printInfo,1))
    .def("PrintWarning", raw_function(printWarning,1))
    .def("GetDefaultOutput", &Logger::GetDefaultOutputInstance, return_value_policy<reference_existing_object>())
    .staticmethod("PrintInfo")
    .staticmethod("PrintWarning")
    .staticmethod("GetDefaultOutput")
    ;

    // Enums for Severity
    enum_<Logger::Severity>("Severity")
    .value("WARNING", Logger::Severity::WARNING)
    .value("INFO", Logger::Severity::INFO)
    .value("DETAIL", Logger::Severity::DETAIL)
    .value("DEBUG", Logger::Severity::DEBUG)
    .value("TRACE", Logger::Severity::TRACE);

    // Enums for Category
    enum_<Logger::Category>("Category")
    .value("STATUS", Logger::Category::STATUS)
    .value("CRITICAL", Logger::Category::CRITICAL)
    .value("STATISTICS", Logger::Category::STATISTICS)
    .value("PROFILING", Logger::Category::PROFILING)
    .value("CHECKING", Logger::Category::CHECKING);
}

}  // namespace Python.
} // Namespace Kratos
