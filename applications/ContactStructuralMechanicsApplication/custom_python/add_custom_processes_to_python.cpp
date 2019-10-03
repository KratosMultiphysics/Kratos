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

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

//Processes
#include "custom_processes/master_slave_process.h"
#include "custom_processes/alm_fast_init_process.h"
#include "custom_processes/alm_variables_calculation_process.h"
#include "custom_processes/contact_spr_error_process.h"
#include "custom_processes/compute_dynamic_factor_process.h"
#include "custom_processes/contact_search_wrapper_process.h"
#include "custom_processes/simple_contact_search_process.h"
#include "custom_processes/advanced_contact_search_process.h"
#include "custom_processes/find_intersected_geometrical_objects_with_obb_for_contact_search_process.h"
#include "custom_processes/normal_gap_process.h"
#include "custom_processes/normal_check_process.h"
#include "custom_processes/mpc_contact_search_process.h"
#include "custom_processes/mpc_contact_search_wrapper_process.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void  AddCustomProcessesToPython(pybind11::module& m)
{
    typedef Process  ProcessBaseType;

    py::class_<ALMFastInit, ALMFastInit::Pointer, ProcessBaseType >
    (m, "ALMFastInit")
    .def(py::init<ModelPart&>())
    ;

    py::class_<MasterSlaveProcess, MasterSlaveProcess::Pointer, ProcessBaseType >
    (m, "MasterSlaveProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<ComputeDynamicFactorProcess, ComputeDynamicFactorProcess::Pointer, ProcessBaseType >
    (m, "ComputeDynamicFactorProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<ALMVariablesCalculationProcess, ALMVariablesCalculationProcess::Pointer, ProcessBaseType >
    (m, "ALMVariablesCalculationProcess")
    .def(py::init<ModelPart&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, Variable<double>&>()) // Considering default variables
    .def(py::init<ModelPart&>())
    ;

    //SPR_ERROR
    py::class_<ContactSPRErrorProcess<2>, ContactSPRErrorProcess<2>::Pointer, Process >(m, "ContactSPRErrorProcess2D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<ContactSPRErrorProcess<3>, ContactSPRErrorProcess<3>::Pointer, Process >(m, "ContactSPRErrorProcess3D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    // Wrapper contact search
    py::class_<ContactSearchWrapperProcess, typename ContactSearchWrapperProcess::Pointer, Process>(m, "ContactSearchProcess")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    // Simple contact search
    py::class_<SimpleContactSearchProcess<2, 2>, typename SimpleContactSearchProcess<2, 2>::Pointer, Process>(m, "SimpleContactSearchProcess2D2N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearchProcess<2, 2>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearchProcess<2, 2>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearchProcess<2, 2>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearchProcess<2, 2>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearchProcess<2, 2>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearchProcess<2, 2>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearchProcess<2, 2>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearchProcess<2, 2>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearchProcess<2, 2>::InvertSearch)
    ;
    py::class_<SimpleContactSearchProcess<3, 3>, typename SimpleContactSearchProcess<3, 3>::Pointer, Process>(m, "SimpleContactSearchProcess3D3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearchProcess<3, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearchProcess<3, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearchProcess<3, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearchProcess<3, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearchProcess<3, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearchProcess<3, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearchProcess<3, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearchProcess<3, 3>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearchProcess<3, 3>::InvertSearch)
    ;
    py::class_<SimpleContactSearchProcess<3, 4>, typename SimpleContactSearchProcess<3, 4>::Pointer, Process>(m, "SimpleContactSearchProcess3D4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearchProcess<3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearchProcess<3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearchProcess<3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearchProcess<3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearchProcess<3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearchProcess<3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearchProcess<3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearchProcess<3, 4>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearchProcess<3, 4>::InvertSearch)
    ;
    py::class_<SimpleContactSearchProcess<3, 3, 4>, typename SimpleContactSearchProcess<3, 3, 4>::Pointer, Process>(m, "SimpleContactSearchProcess3D3N4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearchProcess<3, 3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearchProcess<3, 3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearchProcess<3, 3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearchProcess<3, 3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearchProcess<3, 3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearchProcess<3, 3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearchProcess<3, 3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearchProcess<3, 3, 4>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearchProcess<3, 3, 4>::InvertSearch)
    ;
    py::class_<SimpleContactSearchProcess<3, 4, 3>, typename SimpleContactSearchProcess<3, 4, 3>::Pointer, Process>(m, "SimpleContactSearchProcess3D4N3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&SimpleContactSearchProcess<3, 4, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&SimpleContactSearchProcess<3, 4, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&SimpleContactSearchProcess<3, 4, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&SimpleContactSearchProcess<3, 4, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&SimpleContactSearchProcess<3, 4, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&SimpleContactSearchProcess<3, 4, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&SimpleContactSearchProcess<3, 4, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&SimpleContactSearchProcess<3, 4, 3>::CheckMortarConditions)
    .def("InvertSearch",&SimpleContactSearchProcess<3, 4, 3>::InvertSearch)
    ;

    // Advanced contact search
    py::class_<AdvancedContactSearchProcess<2, 2>, typename AdvancedContactSearchProcess<2, 2>::Pointer, Process>(m, "AdvancedContactSearchProcess2D2N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearchProcess<2, 2>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearchProcess<2, 2>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearchProcess<2, 2>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearchProcess<2, 2>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearchProcess<2, 2>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearchProcess<2, 2>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearchProcess<2, 2>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearchProcess<2, 2>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearchProcess<2, 2>::InvertSearch)
    ;
    py::class_<AdvancedContactSearchProcess<3, 3>, typename AdvancedContactSearchProcess<3, 3>::Pointer, Process>(m, "AdvancedContactSearchProcess3D3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearchProcess<3, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearchProcess<3, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearchProcess<3, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearchProcess<3, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearchProcess<3, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearchProcess<3, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearchProcess<3, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearchProcess<3, 3>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearchProcess<3, 3>::InvertSearch)
    ;
    py::class_<AdvancedContactSearchProcess<3, 4>, typename AdvancedContactSearchProcess<3, 4>::Pointer, Process>(m, "AdvancedContactSearchProcess3D4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearchProcess<3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearchProcess<3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearchProcess<3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearchProcess<3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearchProcess<3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearchProcess<3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearchProcess<3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearchProcess<3, 4>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearchProcess<3, 4>::InvertSearch)
    ;
    py::class_<AdvancedContactSearchProcess<3, 3, 4>, typename AdvancedContactSearchProcess<3, 3, 4>::Pointer, Process>(m, "AdvancedContactSearchProcess3D3N4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearchProcess<3, 3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearchProcess<3, 3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearchProcess<3, 3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearchProcess<3, 3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearchProcess<3, 3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearchProcess<3, 3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearchProcess<3, 3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearchProcess<3, 3, 4>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearchProcess<3, 3, 4>::InvertSearch)
    ;
    py::class_<AdvancedContactSearchProcess<3, 4, 3>, typename AdvancedContactSearchProcess<3, 4, 3>::Pointer, Process>(m, "AdvancedContactSearchProcess3D4N3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&AdvancedContactSearchProcess<3, 4, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&AdvancedContactSearchProcess<3, 4, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&AdvancedContactSearchProcess<3, 4, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&AdvancedContactSearchProcess<3, 4, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&AdvancedContactSearchProcess<3, 4, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&AdvancedContactSearchProcess<3, 4, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&AdvancedContactSearchProcess<3, 4, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&AdvancedContactSearchProcess<3, 4, 3>::CheckMortarConditions)
    .def("InvertSearch",&AdvancedContactSearchProcess<3, 4, 3>::InvertSearch)
    ;

    // Wrapper contact search
    py::class_<MPCContactSearchWrapperProcess, typename MPCContactSearchWrapperProcess::Pointer, Process>(m, "MPCContactSearchProcess")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    // MPC contact search
    py::class_<MPCContactSearchProcess<2, 2>, typename MPCContactSearchProcess<2, 2>::Pointer, Process>(m, "MPCContactSearchProcess2D2N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&MPCContactSearchProcess<2, 2>::InitializeMortarConditions)
    .def("ClearMortarConditions",&MPCContactSearchProcess<2, 2>::ClearMortarConditions)
    .def("CheckContactModelParts",&MPCContactSearchProcess<2, 2>::CheckContactModelParts)
    .def("CreatePointListMortar",&MPCContactSearchProcess<2, 2>::CreatePointListMortar)
    .def("UpdatePointListMortar",&MPCContactSearchProcess<2, 2>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&MPCContactSearchProcess<2, 2>::UpdateMortarConditions)
    .def("ResetContactOperators",&MPCContactSearchProcess<2, 2>::ResetContactOperators)
    .def("CheckMortarConditions",&MPCContactSearchProcess<2, 2>::CheckMortarConditions)
    .def("InvertSearch",&MPCContactSearchProcess<2, 2>::InvertSearch)
    ;
    py::class_<MPCContactSearchProcess<3, 3>, typename MPCContactSearchProcess<3, 3>::Pointer, Process>(m, "MPCContactSearchProcess3D3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&MPCContactSearchProcess<3, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&MPCContactSearchProcess<3, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&MPCContactSearchProcess<3, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&MPCContactSearchProcess<3, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&MPCContactSearchProcess<3, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&MPCContactSearchProcess<3, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&MPCContactSearchProcess<3, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&MPCContactSearchProcess<3, 3>::CheckMortarConditions)
    .def("InvertSearch",&MPCContactSearchProcess<3, 3>::InvertSearch)
    ;
    py::class_<MPCContactSearchProcess<3, 4>, typename MPCContactSearchProcess<3, 4>::Pointer, Process>(m, "MPCContactSearchProcess3D4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&MPCContactSearchProcess<3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&MPCContactSearchProcess<3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&MPCContactSearchProcess<3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&MPCContactSearchProcess<3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&MPCContactSearchProcess<3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&MPCContactSearchProcess<3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&MPCContactSearchProcess<3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&MPCContactSearchProcess<3, 4>::CheckMortarConditions)
    .def("InvertSearch",&MPCContactSearchProcess<3, 4>::InvertSearch)
    ;
    py::class_<MPCContactSearchProcess<3, 3, 4>, typename MPCContactSearchProcess<3, 3, 4>::Pointer, Process>(m, "MPCContactSearchProcess3D3N4N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&MPCContactSearchProcess<3, 3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&MPCContactSearchProcess<3, 3, 4>::ClearMortarConditions)
    .def("CheckContactModelParts",&MPCContactSearchProcess<3, 3, 4>::CheckContactModelParts)
    .def("CreatePointListMortar",&MPCContactSearchProcess<3, 3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&MPCContactSearchProcess<3, 3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&MPCContactSearchProcess<3, 3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&MPCContactSearchProcess<3, 3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&MPCContactSearchProcess<3, 3, 4>::CheckMortarConditions)
    .def("InvertSearch",&MPCContactSearchProcess<3, 3, 4>::InvertSearch)
    ;
    py::class_<MPCContactSearchProcess<3, 4, 3>, typename MPCContactSearchProcess<3, 4, 3>::Pointer, Process>(m, "MPCContactSearchProcess3D4N3N")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&MPCContactSearchProcess<3, 4, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&MPCContactSearchProcess<3, 4, 3>::ClearMortarConditions)
    .def("CheckContactModelParts",&MPCContactSearchProcess<3, 4, 3>::CheckContactModelParts)
    .def("CreatePointListMortar",&MPCContactSearchProcess<3, 4, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&MPCContactSearchProcess<3, 4, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&MPCContactSearchProcess<3, 4, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&MPCContactSearchProcess<3, 4, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&MPCContactSearchProcess<3, 4, 3>::CheckMortarConditions)
    .def("InvertSearch",&MPCContactSearchProcess<3, 4, 3>::InvertSearch)
    ;

    // Normal gap process
    py::class_<NormalGapProcess<2, 2>, typename NormalGapProcess<2, 2>::Pointer, Process>(m, "NormalGapProcess2D2N")
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, const bool>())
    ;
    py::class_<NormalGapProcess<3, 3>, typename NormalGapProcess<3, 3>::Pointer, Process>(m, "NormalGapProcess3D3N")
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, const bool>())
    ;
    py::class_<NormalGapProcess<3, 4>, typename NormalGapProcess<3, 4>::Pointer, Process>(m, "NormalGapProcess3D4N")
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, const bool>())
    ;
    py::class_<NormalGapProcess<3, 3, 4>, typename NormalGapProcess<3, 3, 4>::Pointer, Process>(m, "NormalGapProcess3D3N4N")
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, const bool>())
    ;
    py::class_<NormalGapProcess<3, 4, 3>, typename NormalGapProcess<3, 4, 3>::Pointer, Process>(m, "NormalGapProcess3D4N3N")
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, const bool>())
    ;

    // Normal check process
    py::class_<NormalCheckProcess, NormalCheckProcess::Pointer, Process>(m, "NormalCheckProcess")
    .def(py::init<ModelPart&>())
    ;

    // FindIntersectedGeometricalObjectsWithOBBContactSearchProcess
    py::class_<FindIntersectedGeometricalObjectsWithOBBContactSearchProcess, FindIntersectedGeometricalObjectsWithOBBContactSearchProcess::Pointer, Process>(m,"FindIntersectedGeometricalObjectsWithOBBContactSearchProcess")
    .def(py::init<ModelPart&,ModelPart&>())
    .def(py::init<ModelPart&,ModelPart&, const double>())
    .def(py::init<Model&, Parameters>())
    ;
}
}  // namespace Python.
} // Namespace Kratos

