//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2018 $
//   Revision:            $Revision:                     0.0 $
//
//

// External includes
#include <pybind11/stl.h>

//Application includes
#include "custom_python/add_custom_processes_to_python.h"

// Project includes
#include "includes/node.h"
#include "processes/process.h"

//PreMeshing processes
#include "includes/model_part.h"

//MiddleMeshing processes
// #include "custom_processes/refine_mesh_elements_on_size_process.hpp"
#include "custom_processes/remove_mesh_nodes_for_fluids_process.hpp"

//PostMeshing processes
#include "custom_processes/recover_volume_losses_process.hpp"
#include "custom_processes/select_mesh_elements_for_fluids_process.hpp"
#include "custom_processes/generate_new_nodes_before_meshing_process.hpp"
#include "custom_processes/inlet_management_process.hpp"
#include "custom_processes/set_inlet_process.hpp"
#include "custom_processes/model_start_end_meshing_for_fluids_process.hpp"
#include "custom_processes/split_elements_process.hpp"
#include "custom_processes/set_active_flag_process.hpp"
#include "custom_processes/set_active_flag_mesher_process.hpp"
#include "custom_processes/set_material_properties_to_fluid_nodes_process.hpp"
#include "custom_processes/set_material_properties_from_fluid_to_rigid_nodes_process.hpp"
#include "custom_processes/set_material_properties_to_solid_nodes_process.hpp"
#include "custom_processes/adaptive_time_interval_process.hpp"
#include "custom_processes/transfer_model_part_elements_process.hpp"
#include "custom_processes/build_mesh_boundary_for_fluids_process.hpp"
#include "custom_processes/build_model_part_boundary_for_fluids_process.hpp"
#include "custom_processes/generate_new_conditions_mesher_for_fluids_process.hpp"
#include "custom_processes/lagrangian_rotation_process.hpp"
#include "custom_processes/compute_average_pfem_mesh_parameters_process.hpp"
#include "custom_processes/fix_scalar_pfem_dof_process.hpp"
#include "custom_processes/free_scalar_pfem_dof_process.hpp"
#include "custom_processes/set_dummy_property_for_rigid_boundaries_process.hpp"

#include "custom_processes/assign_scalar_variable_to_pfem_entities_process.hpp"
#include "custom_processes/assign_vector_variable_to_pfem_conditions_process.hpp"
#include "custom_processes/assign_vector_field_to_pfem_entities_process.hpp"
#include "custom_processes/assign_scalar_field_to_pfem_entities_process.hpp"

//Processes

namespace Kratos
{

namespace Python
{

typedef std::vector<Flags> FlagsContainer;

void AddCustomProcessesToPython(pybind11::module &m)
{

    namespace py = pybind11;
    typedef Process ProcessBaseType;
    typedef SettleModelStructureProcess ModelStartEndMeshingProcessType;

    py::class_<RecoverVolumeLossesProcess, RecoverVolumeLossesProcess::Pointer, MesherProcess>(m, "RecoverVolumeLosses")
        .def(py::init<ModelPart &, MesherUtilities::MeshingParameters &, int>());

    py::class_<RemoveMeshNodesForFluidsProcess, RemoveMeshNodesForFluidsProcess::Pointer, MesherProcess>(m, "RemoveMeshNodesForFluids")
        .def(py::init<ModelPart &, MesherUtilities::MeshingParameters &, int>());

    py::class_<GenerateNewNodesBeforeMeshingProcess, GenerateNewNodesBeforeMeshingProcess::Pointer, MesherProcess>(m, "GenerateNewNodesBeforeMeshing")
        .def(py::init<ModelPart &, MesherUtilities::MeshingParameters &, int>());

    py::class_<SelectMeshElementsForFluidsProcess, SelectMeshElementsForFluidsProcess::Pointer, MesherProcess>(m, "SelectMeshElementsForFluids")
        .def(py::init<ModelPart &, MesherUtilities::MeshingParameters &, int>());

    py::class_<InletManagementProcess, InletManagementProcess::Pointer, MesherProcess>(m, "InletManagement")
        .def(py::init<ModelPart &, MesherUtilities::MeshingParameters &, int>());

    py::class_<SetInletProcess, SetInletProcess::Pointer, ProcessBaseType>(m, "SetInlet")
        .def(py::init<ModelPart &, int>());

    py::class_<SplitElementsProcess, SplitElementsProcess::Pointer, ProcessBaseType>(m, "SplitElementsProcess")
        .def(py::init<ModelPart &, int>());

    py::class_<SetActiveFlagProcess, SetActiveFlagProcess::Pointer, MesherProcess>(m, "SetActiveFlagProcess")
        .def(py::init<ModelPart &, bool, bool, int>());

    py::class_<SetMaterialPropertiesToFluidNodesProcess, SetMaterialPropertiesToFluidNodesProcess::Pointer, MesherProcess>(m, "SetMaterialPropertiesToFluidNodes")
        .def(py::init<ModelPart &>());

    py::class_<SetMaterialPropertiesFromFluidToRigidNodesProcess, SetMaterialPropertiesFromFluidToRigidNodesProcess::Pointer, MesherProcess>(m, "SetMaterialPropertiesFromFluidToRigidNodes")
        .def(py::init<ModelPart &, ModelPart &>());

    py::class_<SetMaterialPropertiesToSolidNodesProcess, SetMaterialPropertiesToSolidNodesProcess::Pointer, MesherProcess>(m, "SetMaterialPropertiesToSolidNodes")
        .def(py::init<ModelPart &>());

    py::class_<SetActiveFlagMesherProcess, SetActiveFlagMesherProcess::Pointer, SetActiveFlagProcess>(m, "SetActiveFlagMesherProcess")
        .def(py::init<ModelPart &, bool, bool, int>());

    py::class_<AdaptiveTimeIntervalProcess, AdaptiveTimeIntervalProcess::Pointer, ProcessBaseType>(m, "AdaptiveTimeIntervalProcess")
        .def(py::init<ModelPart &, int>());

    py::class_<ModelStartEndMeshingForFluidsProcess, ModelStartEndMeshingForFluidsProcess::Pointer, ModelStartEndMeshingProcessType>(m, "ModelMeshingForFluids")
        .def(py::init<ModelPart &, Flags, int>());

    py::class_<BuildMeshBoundaryForFluidsProcess, BuildMeshBoundaryForFluidsProcess::Pointer, MesherProcess>(m, "BuildMeshBoundaryForFluids")
        .def(py::init<ModelPart &, MesherUtilities::MeshingParameters &, int>());

    py::class_<BuildModelPartBoundaryForFluidsProcess, BuildModelPartBoundaryForFluidsProcess::Pointer, MesherProcess>(m, "BuildModelPartBoundaryForFluids")
        .def(py::init<ModelPart &, std::string, int>())
        .def("SearchConditionMasters", &BuildModelPartBoundaryForFluidsProcess::SearchConditionMasters);

    py::class_<SetDummyPropertyForRigidElementsProcess, SetDummyPropertyForRigidElementsProcess::Pointer, ProcessBaseType>(m, "SetDummyPropertyForRigidElementsProcess")
        .def(py::init<ModelPart &, unsigned int &>())
        .def("Execute", &SetDummyPropertyForRigidElementsProcess::Execute);

    //**********TRANSFER ELEMENTS TO MODEL PART*********//
    py::class_<TransferModelPartElementsProcess, TransferModelPartElementsProcess::Pointer, ProcessBaseType>(m, "TransferModelPartElementsProcess")
        .def(py::init<ModelPart &, ModelPart &>())
        .def("Execute", &TransferModelPartElementsProcess::Execute);

    py::class_<GenerateNewConditionsMesherForFluidsProcess, GenerateNewConditionsMesherForFluidsProcess::Pointer, BuildModelPartBoundaryProcess>(m, "GenerateNewConditionsForFluids")
        .def(py::init<ModelPart &, MesherUtilities::MeshingParameters &, int>());

    py::class_<LagrangianRotationProcess, LagrangianRotationProcess::Pointer, ProcessBaseType>(m, "LagrangianRotationProcess")
        .def(py::init<ModelPart &, Parameters>());

    py::class_<ComputeAveragePfemMeshParametersProcess, ComputeAveragePfemMeshParametersProcess::Pointer, MesherProcess>(m, "ComputeAveragePfemMeshParameters")
        .def(py::init<ModelPart &, MesherUtilities::MeshingParameters &, int>());

    //**********FIX AND FREE DOFS PROCESSES*********//

    py::class_<FixScalarPfemDofProcess, FixScalarPfemDofProcess::Pointer, Process>(m, "FixScalarPfemDofProcess")
        .def(py::init<ModelPart &, Parameters>())
        .def(py::init<ModelPart &, Parameters &>())
        .def(py::init<ModelPart &, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> &>())
        .def(py::init<ModelPart &, const Variable<double> &>())
        .def(py::init<ModelPart &, const Variable<int> &>())
        .def(py::init<ModelPart &, const Variable<bool> &>())
        .def("Execute", &FixScalarPfemDofProcess::Execute)

        ;

    py::class_<FreeScalarPfemDofProcess, FreeScalarPfemDofProcess::Pointer, Process>(m, "FreeScalarPfemDofProcess")
        .def(py::init<ModelPart &, Parameters>())
        .def(py::init<ModelPart &, Parameters &>())
        .def(py::init<ModelPart &, const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> &>())
        .def(py::init<ModelPart &, const Variable<double> &>())
        .def(py::init<ModelPart &, const Variable<int> &>())
        .def(py::init<ModelPart &, const Variable<bool> &>())
        .def("Execute", &FreeScalarPfemDofProcess::Execute)

        ;

    // //**********ASSIGN VALUES TO VARIABLES PROCESSES*********//

    py::class_<AssignScalarVariableToPfemEntitiesProcess, AssignScalarVariableToPfemEntitiesProcess::Pointer, Process>(m, "AssignScalarToEntitiesProcess")
        .def(py::init<ModelPart &, Parameters>())
        .def(py::init<ModelPart &, Parameters &>())
        .def("Execute", &AssignScalarVariableToPfemEntitiesProcess::Execute);

    py::class_<AssignScalarFieldToPfemEntitiesProcess, AssignScalarFieldToPfemEntitiesProcess::Pointer, AssignScalarVariableToPfemEntitiesProcess>(m, "AssignScalarFieldToEntitiesProcess")
        .def(py::init<ModelPart &, pybind11::object &, const std::string, const bool, Parameters>())
        .def(py::init<ModelPart &, pybind11::object &, const std::string, const bool, Parameters &>())
        .def("Execute", &AssignScalarFieldToPfemEntitiesProcess::Execute);

    py::class_<AssignVectorFieldToPfemEntitiesProcess, AssignVectorFieldToPfemEntitiesProcess::Pointer, AssignScalarFieldToPfemEntitiesProcess>(m, "AssignVectorFieldToEntitiesProcess")
        .def(py::init<ModelPart &, pybind11::object &, const std::string, const bool, Parameters>())
        .def(py::init<ModelPart &, pybind11::object &, const std::string, const bool, Parameters &>())
        .def("Execute", &AssignVectorFieldToPfemEntitiesProcess::Execute);

    py::class_<AssignVectorVariableToPfemConditionsProcess, AssignVectorVariableToPfemConditionsProcess::Pointer, AssignScalarVariableToPfemEntitiesProcess>(m, "AssignVectorToConditionsProcess")
        .def(py::init<ModelPart &, Parameters>())
        .def(py::init<ModelPart &, Parameters &>())
        .def(py::init<ModelPart &, const Variable<array_1d<double, 3>> &, array_1d<double, 3> &>())
        .def("Execute", &AssignVectorVariableToPfemConditionsProcess::Execute);
    ;
}

} // namespace Python.

} // Namespace Kratos
