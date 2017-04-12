// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
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
#include "custom_processes/alm_variables_calculation_process.h"

namespace Kratos
{
    namespace Python
    {
        void  AddCustomProcessesToPython()
        {
            using namespace boost::python;
            typedef Process  ProcessBaseType;

            class_<ALMVariablesCalculationProcess, bases<ProcessBaseType>, boost::noncopyable >
            (
                "ALMVariablesCalculationProcess", init<ModelPart&, Variable<double>&>()
            )
            .def(init<ModelPart&, ModelPart&>()) // Considering default variables
            .def("Execute", &ALMVariablesCalculationProcess::Execute)
            ;
        }
    }  // namespace Python.
} // Namespace Kratos

