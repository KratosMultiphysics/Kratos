// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes 
#include <boost/python.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "containers/flags.h"

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/master_slave_process.h"
#include "custom_processes/alm_fast_init_process.h"
#include "custom_processes/alm_variables_calculation_process.h"

namespace Kratos
{
    namespace Python
    {
        void  AddCustomProcessesToPython()
        {
            using namespace boost::python;
            typedef Process  ProcessBaseType;

            class_<ALMFastInit, bases<ProcessBaseType>, boost::noncopyable >
            (
                "ALMFastInit", init<ModelPart&>()
            )
            .def("Execute", &ALMFastInit::Execute)
            ;
            
            class_<MasterSlaveProcess, bases<ProcessBaseType>, boost::noncopyable >
            (
                "MasterSlaveProcess", init<ModelPart&>()
            )
            .def("Execute", &MasterSlaveProcess::Execute)
            ;
            
            class_<ALMVariablesCalculationProcess, bases<ProcessBaseType>, boost::noncopyable >
            (
                "ALMVariablesCalculationProcess", init<ModelPart&, Variable<double>&, Parameters>()
            )
            .def(init<ModelPart&, Variable<double>&>()) // Considering default variables
            .def(init<ModelPart&>()) 
            .def("Execute", &ALMVariablesCalculationProcess::Execute)
            ;
        }
    }  // namespace Python.
} // Namespace Kratos

