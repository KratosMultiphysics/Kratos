//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
//

// Application includes
#include "custom_processes/apply_chimera_process_fractional_step.h"

namespace Kratos {

template <int TDim>
ApplyChimeraProcessFractionalStep<TDim>::ApplyChimeraProcessFractionalStep(ModelPart& rMainModelPart,
                                                                           Parameters iParameters)
    : BaseType(rMainModelPart, iParameters)
{
}

template <int TDim>
void ApplyChimeraProcessFractionalStep<TDim>::ExecuteFinalizeSolutionStep()
{
    if (BaseType::mReformulateEveryStep) {
        auto& vel_modelpart =
            BaseType::mrMainModelPart.GetSubModelPart(BaseType::mrMainModelPart.Name()+"fs_velocity_model_part");
        vel_modelpart.MasterSlaveConstraints().clear();

        auto& pre_modelpart =
            BaseType::mrMainModelPart.GetSubModelPart(BaseType::mrMainModelPart.Name()+"fs_pressure_model_part");
        pre_modelpart.MasterSlaveConstraints().clear();
    }

    BaseType::ExecuteFinalizeSolutionStep();
}

template <int TDim>
std::string ApplyChimeraProcessFractionalStep<TDim>::Info() const
{
    return "ApplyChimeraProcessFractionalStep";
}

template <int TDim>
void ApplyChimeraProcessFractionalStep<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ApplyChimeraProcessFractionalStep" << std::endl;
}

template <int TDim>
void ApplyChimeraProcessFractionalStep<TDim>::PrintData(std::ostream& rOStream) const
{
    KRATOS_INFO("ApplyChimeraProcessFractionalStep") << std::endl;
}

template <int TDim>
void ApplyChimeraProcessFractionalStep<TDim>::ApplyContinuityWithMpcs(ModelPart& rBoundaryModelPart,
                                                                      PointLocatorType& rBinLocator)
{
    // loop over nodes and find the triangle in which it falls, then do
    // interpolation
    MasterSlaveContainerVectorType velocity_ms_container_vector;
    MasterSlaveContainerVectorType pressure_ms_container_vector;
    BaseType::ReserveMemoryForConstraintContainers(
        rBoundaryModelPart, velocity_ms_container_vector);
    BaseType::ReserveMemoryForConstraintContainers(
        rBoundaryModelPart, pressure_ms_container_vector);

    const int n_boundary_nodes = static_cast<int>(rBoundaryModelPart.Nodes().size());
    for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn) {
        ModelPart::NodesContainerType::iterator i_boundary_node =
            rBoundaryModelPart.NodesBegin() + i_bn;
        Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());

        BaseType::mNodeIdToConstraintIdsMap[p_boundary_node->Id()].reserve(150);
    }

    BaseType::FormulateConstraints(rBoundaryModelPart, rBinLocator, velocity_ms_container_vector,
                                   pressure_ms_container_vector);

    BuiltinTimer mpc_add_time;
    auto& vel_modelpart =
        BaseType::mrMainModelPart.GetSubModelPart(BaseType::mrMainModelPart.Name()+"fs_velocity_model_part");
    BaseType::AddConstraintsToModelpart(vel_modelpart, velocity_ms_container_vector);
    VariableUtils().SetFlag(FS_CHIMERA_PRESSURE_CONSTRAINT, false,
                            vel_modelpart.MasterSlaveConstraints());
    VariableUtils().SetFlag(FS_CHIMERA_VELOCITY_CONSTRAINT, true,
                            vel_modelpart.MasterSlaveConstraints());
    VariableUtils().SetFlag(ACTIVE, true, vel_modelpart.MasterSlaveConstraints());

    auto& pre_modelpart =
        BaseType::mrMainModelPart.GetSubModelPart(BaseType::mrMainModelPart.Name()+"fs_pressure_model_part");
    BaseType::AddConstraintsToModelpart(pre_modelpart, pressure_ms_container_vector);
    VariableUtils().SetFlag(FS_CHIMERA_PRESSURE_CONSTRAINT, true,
                            vel_modelpart.MasterSlaveConstraints());
    VariableUtils().SetFlag(FS_CHIMERA_VELOCITY_CONSTRAINT, false,
                            vel_modelpart.MasterSlaveConstraints());
    VariableUtils().SetFlag(ACTIVE, true, pre_modelpart.MasterSlaveConstraints());

    KRATOS_INFO_IF("Adding of MPCs from containers to modelpart took         : ", BaseType::mEchoLevel > 1)
        << mpc_add_time.ElapsedSeconds() << " seconds" << std::endl;
}

// Template declarations
template class ApplyChimeraProcessFractionalStep<2>;
template class ApplyChimeraProcessFractionalStep<3>;

} // namespace Kratos.
