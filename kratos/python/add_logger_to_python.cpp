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

// Project includes
#include "includes/define_python.h"
#include "includes/data_communicator.h"
#include "input_output/logger.h"
#include "input_output/file_logger_output.h"


namespace Kratos {
namespace Python {

namespace py = pybind11;

const DataCommunicator& getDataCommunicator(pybind11::kwargs kwargs) {
    if (kwargs.contains("data_communicator")) {
        const DataCommunicator& r_data_communicator = py::cast<DataCommunicator&>(kwargs["data_communicator"]);
        return r_data_communicator;
    }
    else {
        return DataCommunicator::GetDefault();
    }
}

namespace Internals{
    std::string GetMessage(pybind11::args Args, bool useKwargLabel){
        if(len(Args) == 0) {
            return ""; 
        }

        std::stringstream buffer;
        // Get the label
        unsigned int to_skip = useKwargLabel ? 0 : 1; //if the kwargs label is false, consider the first entry of the args as the label

        unsigned int counter = 0;
        for(auto item : Args)
        {
            if(counter >= to_skip)
            {
                buffer << item;
                if(counter < len(Args))
                    buffer << " ";
            }
            counter++;
        }
        return buffer.str();
    }

    std::string GetLabel(pybind11::args Args, pybind11::kwargs KWargs, bool useKwargLabel){
        if(useKwargLabel) {
            if(KWargs.contains("label")) {
                return py::str(KWargs["label"]);
            } else {
                return "";
            }
        } 

        if(len(Args) == 0) {
            return ""; 
        }

        return py::str(Args[0]); // Consider the first entry of the args as the label
   }

    Logger::Severity GetSeverity(pybind11::kwargs KWargs, Logger::Severity DefaultSeverity){
        if(KWargs.contains("severity")) {
            DefaultSeverity = py::cast<Logger::Severity>(KWargs["severity"]);
        }

        return DefaultSeverity;
    }

    Flags GetCategory(pybind11::kwargs KWargs, Flags DefaultCategory){
        if(KWargs.contains("category")) {
            DefaultCategory = py::cast<Flags>(KWargs["category"]);
        }

        return DefaultCategory;
    }
}

/**
 * Prints the arguments from the python script using the Kratos Logger class. Implementation
 * @args tuple  representing the arguments of the function The first argument is the label
 * @kwargs dictionary  resenting key-value pairs for
 * @severity Logger::Severity The message level of severity @see Logger::Severity
 * @useKwargLabel bool Indicates if the label must be gather from kwargs (true) or is the first argument of the call (false)
 * name arguments
 * @printRank bool record the MPI rank in the output message.
 **/
void printImpl(pybind11::args args, pybind11::kwargs kwargs, Logger::Severity severity, bool useKwargLabel, LoggerMessage::DistributedFilter filterOption) {
    // Send the message and options to the logger
    Logger logger(Internals::GetLabel(args, kwargs, useKwargLabel));
    logger << Internals::GetMessage(args, useKwargLabel) << Internals::GetSeverity(kwargs, severity) << Internals::GetCategory(kwargs, LoggerMessage::STATUS) << std::endl;
    logger << filterOption << getDataCommunicator(kwargs);
}

bool isPrintingRank(pybind11::kwargs kwargs) {
    const DataCommunicator& r_data_communicator = getDataCommunicator(kwargs);
    return r_data_communicator.Rank() == 0;
}

/**
 * Prints the arguments from the python script using the Kratos Logger class. Default function uses INFO severity.
 * @args pybind11::args pybind11::object representing the arguments of the function The first argument is the label
 * @kwargs pybind11::dictionary of pybind11::objects resenting key-value pairs for
 * name arguments
 **/
void printDefault(pybind11::args args, pybind11::kwargs kwargs) {
    if (isPrintingRank(kwargs)) {
        printImpl(args, kwargs, Logger::Severity::INFO, true, LoggerMessage::DistributedFilter::FromRoot());
    }
}

/**
 * Prints the arguments from the python script using the Kratos Logger class using INFO severity.
 * @args pybind11::args pybind11::object representing the arguments of the function The first argument is the label
 * @kwargs pybind11::dictionary of pybind11::objects resenting key-value pairs for
 * name arguments
 **/
void printInfo(pybind11::args args, pybind11::kwargs kwargs) {
    if (isPrintingRank(kwargs)) {
        printImpl(args, kwargs, Logger::Severity::INFO, false, LoggerMessage::DistributedFilter::FromRoot());
    }
}

/**
 * Prints the arguments from the python script using the Kratos Logger class using WARNING severity.
 * @args pybind11::args pybind11::object representing the arguments of the function The first argument is the label
 * @kwargs pybind11::dictionary of pybind11::objects resenting key-value pairs for
 * name arguments
 **/
void printWarning(pybind11::args args, pybind11::kwargs kwargs) {
    if (isPrintingRank(kwargs)) {
        printImpl(args, kwargs, Logger::Severity::WARNING, false, LoggerMessage::DistributedFilter::FromRoot());
    }
}

void printDefaultOnAllRanks(pybind11::args args, pybind11::kwargs kwargs) {
    printImpl(args, kwargs, Logger::Severity::INFO, true, LoggerMessage::DistributedFilter::FromAllRanks());
}

void printInfoOnAllRanks(pybind11::args args, pybind11::kwargs kwargs) {
    printImpl(args, kwargs, Logger::Severity::INFO, false, LoggerMessage::DistributedFilter::FromAllRanks());
}

void printWarningOnAllRanks(pybind11::args args, pybind11::kwargs kwargs) {
    printImpl(args, kwargs, Logger::Severity::WARNING, false, LoggerMessage::DistributedFilter::FromAllRanks());
}

void LoggerStart(pybind11::args args, pybind11::kwargs kwargs){
     if (isPrintingRank(kwargs)) {
        Logger logger(Internals::GetLabel(args, kwargs, false));
        logger.Start() << Internals::GetMessage(args, false) << Internals::GetSeverity(kwargs, Logger::Severity::INFO) << Internals::GetCategory(kwargs, LoggerMessage::STATUS) << std::endl
         << LoggerMessage::DistributedFilter::FromRoot() << getDataCommunicator(kwargs);
    }
}

void LoggerStop(pybind11::args args, pybind11::kwargs kwargs){
     if (isPrintingRank(kwargs)) {
        Logger logger(Internals::GetLabel(args, kwargs, false));
        logger.Stop() << Internals::GetMessage(args, false) << Internals::GetSeverity(kwargs, Logger::Severity::INFO) << Internals::GetCategory(kwargs, LoggerMessage::STATUS) << std::endl
         << LoggerMessage::DistributedFilter::FromRoot() << getDataCommunicator(kwargs);
    }
}

void  AddLoggerToPython(pybind11::module& m) {

    auto logger_output = py::class_<LoggerOutput, Kratos::shared_ptr<LoggerOutput>>(m,"LoggerOutput")
    .def("SetMaxLevel", &LoggerOutput::SetMaxLevel)
    .def("GetMaxLevel", &LoggerOutput::GetMaxLevel)
    .def("SetSeverity", &LoggerOutput::SetSeverity)
    .def("GetSeverity", &LoggerOutput::GetSeverity)
    .def("SetCategory", &LoggerOutput::SetFlags)
    .def("GetCategory", &LoggerOutput::GetFlags)
    .def("SetOption", &LoggerOutput::SetOption)
    .def("GetOption", &LoggerOutput::GetOption)
    ;
    logger_output.attr("WARNING_PREFIX") = LoggerOutput::WARNING_PREFIX;
    logger_output.attr("INFO_PREFIX") = LoggerOutput::INFO_PREFIX;
    logger_output.attr("DETAIL_PREFIX") = LoggerOutput::DETAIL_PREFIX;
    logger_output.attr("DEBUG_PREFIX") = LoggerOutput::DEBUG_PREFIX;
    logger_output.attr("TRACE_PREFIX") = LoggerOutput::TRACE_PREFIX;

    py::class_<FileLoggerOutput, Kratos::shared_ptr<FileLoggerOutput>, LoggerOutput>(m,"FileLoggerOutput")
    .def(py::init<std::string>())
    ;

    py::class_<Logger, Kratos::shared_ptr<Logger>> logger_scope(m,"Logger");
    logger_scope.def(py::init<std::string const &>());
    logger_scope.def_static("Print", printDefault); // raw_function(printDefault,1))
    logger_scope.def_static("PrintInfo",printInfo); // raw_function(printInfo,1))
    logger_scope.def_static("PrintWarning", printWarning); //raw_function(printWarning,1))
    logger_scope.def_static("PrintOnAllRanks", printDefaultOnAllRanks);
    logger_scope.def_static("PrintInfoOnAllRanks",printInfoOnAllRanks);
    logger_scope.def_static("PrintWarningOnAllRanks", printWarningOnAllRanks);
    logger_scope.def_static("Flush", Logger::Flush);
    logger_scope.def_static("GetDefaultOutput", &Logger::GetDefaultOutputInstance, py::return_value_policy::reference); //_internal )
    logger_scope.def_static("AddOutput", &Logger::AddOutput);
    logger_scope.def_static("Start", LoggerStart); // raw_function(printDefault,1))
    logger_scope.def_static("Stop", LoggerStop); // raw_function(printDefault,1))
    ;

    // Enums for Severity
    py::enum_<Logger::Severity>(logger_scope,"Severity")
    .value("WARNING", Logger::Severity::WARNING)
    .value("INFO", Logger::Severity::INFO)
    .value("DETAIL", Logger::Severity::DETAIL)
    .value("DEBUG", Logger::Severity::DEBUG)
    .value("TRACE", Logger::Severity::TRACE);

    logger_output.attr("STATUS") = LoggerMessage::STATUS;
    logger_output.attr("CRITICAL") = LoggerMessage::CRITICAL;
    logger_output.attr("STATISTICS") = LoggerMessage::STATISTICS;
    logger_output.attr("PROFILING") = LoggerMessage::PROFILING;
    logger_output.attr("CHECKING") = LoggerMessage::CHECKING;


}

}  // namespace Python.
} // Namespace Kratos
