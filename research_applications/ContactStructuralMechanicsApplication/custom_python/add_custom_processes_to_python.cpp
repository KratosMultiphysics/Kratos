// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
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
#include "custom_python/add_custom_processes_to_python.h"

/* Processes */
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
#include "custom_processes/assign_parent_element_conditions_process.h"

namespace Kratos::Python
{
namespace py = pybind11;

/**
 * @brief RegisterContactSearchProcess is a function that registers a contact search process
 * in the given py::module object.
 * @param m the py::module object where the contact search process will be registered
 * @param rName the name of the contact search process
 */
template<class TClass>
void RegisterContactSearchProcess(py::module& m, const std::string& rName)
{
    py::class_<TClass, typename TClass::Pointer, Process>(m, rName.c_str())
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def(py::init<ModelPart&, Parameters, Properties::Pointer>())
    .def("InitializeMortarConditions", &TClass::InitializeMortarConditions)
    .def("ClearMortarConditions", &TClass::ClearMortarConditions)
    .def("CheckContactModelParts", &TClass::CheckContactModelParts)
    .def("CreatePointListMortar", &TClass::CreatePointListMortar)
    .def("UpdatePointListMortar", &TClass::UpdatePointListMortar)
    .def("UpdateMortarConditions", &TClass::UpdateMortarConditions)
    .def("ResetContactOperators", &TClass::ResetContactOperators)
    .def("CheckMortarConditions", &TClass::CheckMortarConditions)
    .def("InvertSearch", &TClass::InvertSearch)
    ;
}

/**
 * @brief Registers simple contact search process.
 * @param m the module to register the laws in
 * @param rEndName the name to append to the law names
 */
template<std::size_t TDim, std::size_t TNumNodes, std::size_t TNumNodesMaster>
void RegisterSimpleContactSearchProcess(py::module& m, const std::string& rEndName)
{
    const std::string name = "SimpleContactSearchProcess" + rEndName;
    RegisterContactSearchProcess<SimpleContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>>(m, name);
}

/**
 * @brief Registers advanced contact search process.
 * @param m the module to register the laws in
 * @param rEndName the name to append to the law names
 */
template<std::size_t TDim, std::size_t TNumNodes, std::size_t TNumNodesMaster>
void RegisterAdvancedContactSearchProcess(py::module& m, const std::string& rEndName)
{
    const std::string name = "AdvancedContactSearchProcess" + rEndName;
    RegisterContactSearchProcess<AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>>(m, name);
}

/**
 * @brief Registers the normal gap process
 * @param m the module to register the laws in
 * @param rEndName the name to append to the law names
 */
template<std::size_t TDim, std::size_t TNumNodes, std::size_t TNumNodesMaster>
void RegisterNormalGapProcess(py::module& m, const std::string& rEndName)
{
    const std::string name = "NormalGapProcess" + rEndName;
    py::class_<NormalGapProcess<TDim, TNumNodes, TNumNodesMaster>, typename NormalGapProcess<TDim, TNumNodes, TNumNodesMaster>::Pointer, Process>(m, name.c_str())
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, const bool>())
    ;
}

/**
 * @brief Registers MPC contact search process.
 * @param m the module to register the laws in
 * @param rEndName the name to append to the law names
 */
template<std::size_t TDim, std::size_t TNumNodes, std::size_t TNumNodesMaster>
void RegisterMPCContactSearchProcess(py::module& m, const std::string& rEndName)
{
    const std::string name = "MPCContactSearchProcess" + rEndName;
    RegisterContactSearchProcess<MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>>(m, name);
}

void  AddCustomProcessesToPython(pybind11::module& m)
{
    py::class_<ALMFastInit, ALMFastInit::Pointer, Process >
    (m, "ALMFastInit")
    .def(py::init<ModelPart&>())
    ;

    py::class_<MasterSlaveProcess, MasterSlaveProcess::Pointer, Process >
    (m, "MasterSlaveProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<ComputeDynamicFactorProcess, ComputeDynamicFactorProcess::Pointer, Process >
    (m, "ComputeDynamicFactorProcess")
    .def(py::init<ModelPart&>())
    ;

    py::class_<ALMVariablesCalculationProcess, ALMVariablesCalculationProcess::Pointer, Process >
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
    .def(py::init<ModelPart&, Parameters, Properties::Pointer>())
    ;

    // Simple contact search
    RegisterSimpleContactSearchProcess<2, 2, 2>(m, "2D2N");
    RegisterSimpleContactSearchProcess<3, 3, 3>(m, "3D3N");
    RegisterSimpleContactSearchProcess<3, 4, 4>(m, "3D4N");
    RegisterSimpleContactSearchProcess<3, 3, 4>(m, "3D3N4N");
    RegisterSimpleContactSearchProcess<3, 4, 3>(m, "3D4N3N");

    // Advanced contact search
    RegisterAdvancedContactSearchProcess<2, 2, 2>(m, "2D2N");
    RegisterAdvancedContactSearchProcess<3, 3, 3>(m, "3D3N");
    RegisterAdvancedContactSearchProcess<3, 4, 4>(m, "3D4N");
    RegisterAdvancedContactSearchProcess<3, 3, 4>(m, "3D3N4N");
    RegisterAdvancedContactSearchProcess<3, 4, 3>(m, "3D4N3N");

    // Wrapper contact search
    py::class_<MPCContactSearchWrapperProcess, typename MPCContactSearchWrapperProcess::Pointer, Process>(m, "MPCContactSearchProcess")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def(py::init<ModelPart&, Parameters, Properties::Pointer>())
    ;

    // MPC contact search
    RegisterMPCContactSearchProcess<2, 2, 2>(m, "2D2N");
    RegisterMPCContactSearchProcess<3, 3, 3>(m, "3D3N");
    RegisterMPCContactSearchProcess<3, 4, 4>(m, "3D4N");
    RegisterMPCContactSearchProcess<3, 3, 4>(m, "3D3N4N");
    RegisterMPCContactSearchProcess<3, 4, 3>(m, "3D4N3N");

    // Normal gap process
    RegisterNormalGapProcess<2, 2, 2>(m, "2D2N");
    RegisterNormalGapProcess<3, 3, 3>(m, "3D3N");
    RegisterNormalGapProcess<3, 4, 4>(m, "3D4N");
    RegisterNormalGapProcess<3, 3, 4>(m, "3D3N4N");
    RegisterNormalGapProcess<3, 4, 3>(m, "3D4N3N");

    // Normal check process
    py::class_<NormalCheckProcess, NormalCheckProcess::Pointer, Process>(m, "NormalCheckProcess")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    // FindIntersectedGeometricalObjectsWithOBBContactSearchProcess
    py::class_<FindIntersectedGeometricalObjectsWithOBBContactSearchProcess, FindIntersectedGeometricalObjectsWithOBBContactSearchProcess::Pointer, Process>(m,"FindIntersectedGeometricalObjectsWithOBBContactSearchProcess")
    .def(py::init<ModelPart&,ModelPart&>())
    .def(py::init<ModelPart&,ModelPart&, const double>())
    .def(py::init<Model&, Parameters>())
    ;

    // AssignParentElementConditionsProcess
    py::class_<AssignParentElementConditionsProcess, AssignParentElementConditionsProcess::Pointer, Process>(m,"AssignParentElementConditionsProcess")
    .def(py::init<ModelPart&,ModelPart&>())
    .def(py::init<Model&, Parameters>())
    ;
}
}  // namespace Kratos::Python.

