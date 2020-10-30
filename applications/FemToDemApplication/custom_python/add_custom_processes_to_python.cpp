//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// External includes

// Project includes

#include "includes/model_part.h"
#include "fem_to_dem_application_variables.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/stress_to_nodes_process.hpp"
#include "custom_processes/compute_normalized_free_energy_on_nodes_process.h"
#include "custom_processes/damage_to_nodes_process.hpp"
#include "custom_processes/dem_after_remesh_identificator_process.hpp"
#include "custom_processes/initial_dem_skin_process.hpp"
#include "custom_processes/extend_pressure_condition_process.h"
#include "custom_processes/assign_pressure_id_process.h"
#include "custom_processes/expand_wet_nodes_process.h"
#include "custom_processes/apply_component_table_process.hpp"
#include "custom_processes/apply_double_table_process.hpp"
#include "custom_processes/generate_dem_process.h"
#include "custom_processes/update_dem_kinematics_process.h"
#include "custom_processes/transfer_nodal_forces_to_fem.h"
#include "custom_processes/compute_sand_production.h"
#include "custom_processes/regenerate_pfem_pressure_conditions_process.h"
#include "custom_processes/update_pressure_value_pfem_conditions_process.h"
#include "custom_processes/fix_free_velocity_on_nodes_process.h"
#include "custom_processes/remove_alone_DEM_elements_process.h"
#include "custom_processes/update_flag_no_remesh_femdem_boundary_process.h"
#include "custom_processes/compute_initial_volume_process.h"
#include "custom_processes/update_pressure_volume_process.h"
#include "custom_processes/generate_initial_skin_DEM_process.h"
#include "custom_processes/transfer_entities_between_model_parts_process.hpp"
#include "custom_processes/fix_scalar_dof_process.hpp"
#include "custom_processes/free_scalar_dof_process.hpp"
#include "custom_processes/assign_scalar_variable_to_entities_process.hpp"
#include "custom_processes/assign_vector_variable_to_conditions_process.hpp"
#include "custom_processes/assign_vector_field_to_entities_process.hpp"
#include "custom_processes/assign_scalar_field_to_entities_process.hpp"

namespace Kratos
{
namespace Python
{
void AddCustomProcessesToPython(pybind11::module &m)
{
    using namespace pybind11;

    typedef Process ProcessBaseType;
    typedef std::vector<Flags>  FlagsContainer;

  class_<AssignScalarVariableToEntitiesProcess, AssignScalarVariableToEntitiesProcess::Pointer, Process>(m,"AssignScalarToEntitiesProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init<ModelPart&, Parameters&>())
      .def("Execute", &AssignScalarVariableToEntitiesProcess::Execute);
      
  class_<AssignScalarFieldToEntitiesProcess, AssignScalarFieldToEntitiesProcess::Pointer, AssignScalarVariableToEntitiesProcess>(m,"AssignScalarFieldToEntitiesProcess")
      .def(init<ModelPart&, pybind11::object&, const std::string, const bool, Parameters>())
      .def(init< ModelPart&, pybind11::object&, const std::string, const bool, Parameters& >())
      .def("Execute", &AssignScalarFieldToEntitiesProcess::Execute);

  class_<AssignVectorFieldToEntitiesProcess, AssignVectorFieldToEntitiesProcess::Pointer, AssignScalarFieldToEntitiesProcess>(m,"AssignVectorFieldToEntitiesProcess")
      .def(init<ModelPart&, pybind11::object&,const std::string,const bool, Parameters>())
      .def(init< ModelPart&, pybind11::object&,const std::string,const bool, Parameters& >())
      .def("Execute", &AssignVectorFieldToEntitiesProcess::Execute);

  class_<AssignVectorVariableToConditionsProcess, AssignVectorVariableToConditionsProcess::Pointer, AssignScalarVariableToEntitiesProcess>(m,"AssignVectorToConditionsProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init< ModelPart&, Parameters& >())
      .def(init<ModelPart&, const Variable<array_1d<double,3> >&, array_1d<double,3>&>())
      .def("Execute", &AssignVectorVariableToConditionsProcess::Execute);



  class_<FixScalarDofProcess, FixScalarDofProcess::Pointer, Process>(m,"FixScalarDofProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init<ModelPart&, Parameters&>())
      .def(init<ModelPart&, const Variable<double>&>())
      .def(init<ModelPart&, const Variable<int>&>())
      .def(init<ModelPart&, const Variable<bool>&>())
      .def("Execute", &FixScalarDofProcess::Execute);


  class_<FreeScalarDofProcess, FreeScalarDofProcess::Pointer, Process>(m,"FreeScalarDofProcess")
      .def(init<ModelPart&, Parameters>())
      .def(init<ModelPart&, Parameters&>())
      .def(init<ModelPart&, const Variable<double>&>())
      .def(init<ModelPart&, const Variable<int>&>())
      .def(init<ModelPart&, const Variable<bool>&>())
      .def("Execute", &FreeScalarDofProcess::Execute);

    // Stress extrapolation to Nodes
    class_<StressToNodesProcess, StressToNodesProcess::Pointer, Process>(m, "StressToNodesProcess")
        .def(init<ModelPart &, unsigned int>())
        .def("Execute", &StressToNodesProcess::Execute);

    // Damage extrapolation to Nodes
    class_<DamageToNodesProcess, DamageToNodesProcess::Pointer, Process>(m, "DamageToNodesProcess")
        .def(init<ModelPart &, unsigned int>())
        .def("Execute", &DamageToNodesProcess::Execute);

    class_<DemAfterRemeshIdentificatorProcess, DemAfterRemeshIdentificatorProcess::Pointer, Process>(m, "DemAfterRemeshIdentificatorProcess")
        .def(init<ModelPart &, const double >())
        .def("Execute", &DemAfterRemeshIdentificatorProcess::Execute);

    class_<InitialDemSkinProcess, InitialDemSkinProcess::Pointer, Process>(m, "InitialDemSkinProcess")
        .def(init<ModelPart &>())
        .def("Execute", &InitialDemSkinProcess::Execute);

    class_<ExtendPressureConditionProcess<2>, ExtendPressureConditionProcess<2>::Pointer, Process>(m, "ExtendPressureConditionProcess2D")
        .def(init<ModelPart &>())
        .def("Execute", &ExtendPressureConditionProcess<2>::Execute);

    class_<ExtendPressureConditionProcess<3>, ExtendPressureConditionProcess<3>::Pointer, Process>(m, "ExtendPressureConditionProcess3D")
        .def(init<ModelPart&>())
        .def("Execute", &ExtendPressureConditionProcess<3>::Execute);

    class_<AssignPressureIdProcess, AssignPressureIdProcess::Pointer, Process>(m, "AssignPressureIdProcess")
        .def(init<ModelPart &>())
        .def("Execute", &AssignPressureIdProcess::Execute);

    class_<ExpandWetNodesProcess, ExpandWetNodesProcess::Pointer, Process>(m, "ExpandWetNodesProcess")
        .def(init<ModelPart &>())
        .def("Execute", &ExpandWetNodesProcess::Execute);

    // Normalized Free Energy extrapolation to Nodes
    class_<ComputeNormalizedFreeEnergyOnNodesProcess, ComputeNormalizedFreeEnergyOnNodesProcess::Pointer, Process>(m, "ComputeNormalizedFreeEnergyOnNodesProcess")
        .def(init<ModelPart &, Parameters>())
        .def("Execute", &ComputeNormalizedFreeEnergyOnNodesProcess::Execute);

    class_<ApplyComponentTableProcess, ApplyComponentTableProcess::Pointer, Process> (m, "ApplyComponentTableProcess")
        .def(init< ModelPart&, Parameters>());

    class_<ApplyDoubleTableProcess, ApplyDoubleTableProcess::Pointer, Process> (m, "ApplyDoubleTableProcess")
        .def(init< ModelPart&, Parameters>());

    class_<GenerateDemProcess, GenerateDemProcess::Pointer, Process>(m, "GenerateDemProcess")
        .def(init<ModelPart &, ModelPart &>())
        .def("Execute", &GenerateDemProcess::Execute);

    class_<UpdateDemKinematicsProcess, UpdateDemKinematicsProcess::Pointer, Process>(m, "UpdateDemKinematicsProcess")
        .def(init<ModelPart &, ModelPart &>())
        .def("Execute", &UpdateDemKinematicsProcess::Execute);

    class_<TransferNodalForcesToFem, TransferNodalForcesToFem::Pointer, Process>(m, "TransferNodalForcesToFem")
        .def(init<ModelPart &, ModelPart &>())
        .def("Execute", &TransferNodalForcesToFem::Execute);

    class_<ComputeSandProduction, ComputeSandProduction::Pointer, Process>(m, "ComputeSandProduction")
        .def(init<ModelPart &>())
        .def("Execute", &ComputeSandProduction::Execute);

    class_<RegeneratePfemPressureConditionsProcess<3>, RegeneratePfemPressureConditionsProcess<3>::Pointer, Process>(m, "RegeneratePfemPressureConditionsProcess3D")
        .def(init<ModelPart &>())
        .def("Execute", &RegeneratePfemPressureConditionsProcess<3>::Execute);

    class_<RegeneratePfemPressureConditionsProcess<2>, RegeneratePfemPressureConditionsProcess<2>::Pointer, Process>(m, "RegeneratePfemPressureConditionsProcess2D")
        .def(init<ModelPart &>())
        .def("Execute", &RegeneratePfemPressureConditionsProcess<2>::Execute);

    class_<UpdatePressureValuePfemConditionsProcess<3>, UpdatePressureValuePfemConditionsProcess<3>::Pointer, Process>(m, "UpdatePressureValuePfemConditionsProcess3D")
        .def(init<ModelPart &>())
        .def("Execute", &UpdatePressureValuePfemConditionsProcess<3>::Execute);

    class_<UpdatePressureValuePfemConditionsProcess<2>, UpdatePressureValuePfemConditionsProcess<2>::Pointer, Process>(m, "UpdatePressureValuePfemConditionsProcess2D")
        .def(init<ModelPart &>())
        .def("Execute", &UpdatePressureValuePfemConditionsProcess<2>::Execute);

    class_<FixFreeVelocityOnNodesProcess, FixFreeVelocityOnNodesProcess::Pointer, Process>(m, "FixFreeVelocityOnNodesProcess")
        .def(init<ModelPart &, const int>())
        .def("Execute", &FixFreeVelocityOnNodesProcess::Execute);

    class_<RemoveAloneDEMElementsProcess, RemoveAloneDEMElementsProcess::Pointer, Process>(m, "RemoveAloneDEMElementsProcess")
        .def(init<ModelPart &, ModelPart &>())
        .def("Execute", &RemoveAloneDEMElementsProcess::Execute);

    class_<UpdateFlagNoRemeshFemDemBoundaryProcess, UpdateFlagNoRemeshFemDemBoundaryProcess::Pointer, Process>(m, "UpdateFlagNoRemeshFemDemBoundaryProcess")
        .def(init<ModelPart &>())
        .def("Execute", &UpdateFlagNoRemeshFemDemBoundaryProcess::Execute);

    class_<ComputeInitialVolumeProcess, ComputeInitialVolumeProcess::Pointer, Process>(m, "ComputeInitialVolumeProcess")
        .def(init<ModelPart &>())
        .def("Execute", &ComputeInitialVolumeProcess::Execute);

    class_<UpdatePressureVolumeProcess, UpdatePressureVolumeProcess::Pointer, Process>(m, "UpdatePressureVolumeProcess")
        .def(init<ModelPart &>())
        .def("Execute", &UpdatePressureVolumeProcess::Execute);

    class_<GenerateInitialSkinDEMProcess, GenerateInitialSkinDEMProcess::Pointer, Process>(m, "GenerateInitialSkinDEMProcess")
        .def(init<ModelPart &, ModelPart &>())
        .def("Execute", &GenerateInitialSkinDEMProcess::Execute);

    class_<TransferEntitiesBetweenModelPartsProcess, TransferEntitiesBetweenModelPartsProcess::Pointer, Process>(m,"TransferEntitiesProcess")
        .def(init<ModelPart&, ModelPart&, const std::string>())
        .def(init<ModelPart&, ModelPart&, const std::string, const FlagsContainer&>())
        .def(init<ModelPart&, ModelPart&, const std::string, const FlagsContainer&, const FlagsContainer& >())
        .def("Execute", &TransferEntitiesBetweenModelPartsProcess::Execute);
}
} // namespace Python.
} // Namespace Kratos