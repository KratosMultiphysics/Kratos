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

// External includes

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "containers/flags.h"

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/master_slave_process.h"
#include "custom_processes/alm_fast_init_process.h"
#include "custom_processes/alm_variables_calculation_process.h"
#include "custom_processes/contact_spr_error_process.h"

namespace Kratos
{
namespace Python
{
using namespace pybind11;

void  AddCustomProcessesToPython(pybind11::module& m)
{
    typedef Process  ProcessBaseType;

    class_<ALMFastInit, ALMFastInit::Pointer, ProcessBaseType >
    (m, "ALMFastInit")
    .def(init<ModelPart&>())
    ;

    class_<MasterSlaveProcess, MasterSlaveProcess::Pointer, ProcessBaseType >
    (m, "MasterSlaveProcess")
    .def(init<ModelPart&>())
    ;

    class_<ALMVariablesCalculationProcess, ALMVariablesCalculationProcess::Pointer, ProcessBaseType >
    (m, "ALMVariablesCalculationProcess")
    .def(init<ModelPart&, Variable<double>&, Parameters>())
    .def(init<ModelPart&, Variable<double>&>()) // Considering default variables
    .def(init<ModelPart&>())
    ;

    //SPR_ERROR
    class_<ContactSPRErrorProcess<2>, ContactSPRErrorProcess<2>::Pointer, Process >(m, "ContactSPRErrorProcess2D")
    .def(init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    ;

    class_<ContactSPRErrorProcess<3>, ContactSPRErrorProcess<3>::Pointer, Process >(m, "ContactSPRErrorProcess3D")
    .def(init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    ;
}
}  // namespace Python.
} // Namespace Kratos

