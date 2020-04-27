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

// Project includes
#include "custom_processes/apply_chimera_process_monolithic.h"

namespace Kratos {


template <int TDim>
    ApplyChimeraProcessMonolithic<TDim>::ApplyChimeraProcessMonolithic(ModelPart& rMainModelPart, Parameters iParameters)
        : BaseType(rMainModelPart, iParameters)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

template <int TDim>
    std::string ApplyChimeraProcessMonolithic<TDim>::Info() const
    {
        return "ApplyChimeraProcessMonolithic";
    }

template <int TDim>
    void ApplyChimeraProcessMonolithic<TDim>::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyChimeraProcessMonolithic" << std::endl;
    }

template <int TDim>
    void ApplyChimeraProcessMonolithic<TDim>::PrintData(std::ostream& rOStream) const
    {
        KRATOS_INFO("ApplyChimeraProcessMonolithic") << std::endl;
    }

template <int TDim>
    void ApplyChimeraProcessMonolithic<TDim>::ApplyContinuityWithMpcs(ModelPart& rBoundaryModelPart, PointLocatorType& rBinLocator)
    {
        // loop over nodes and find the triangle in which it falls, then do
        // interpolation
        MasterSlaveContainerVectorType master_slave_container_vector;
        BaseType::ReserveMemoryForConstraintContainers(
            rBoundaryModelPart, master_slave_container_vector);
        std::vector<int> constraints_id_vector;
        int num_constraints_required = (TDim + 1) * (rBoundaryModelPart.Nodes().size());
        BaseType::CreateConstraintIds(constraints_id_vector, num_constraints_required);

        const int n_boundary_nodes = static_cast<int>(rBoundaryModelPart.Nodes().size());
        for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn) {
            ModelPart::NodesContainerType::iterator i_boundary_node =
                rBoundaryModelPart.NodesBegin() + i_bn;
            BaseType::mNodeIdToConstraintIdsMap[i_boundary_node->Id()].reserve(150);
        }

        BaseType::FormulateConstraints(rBoundaryModelPart, rBinLocator,
                                       master_slave_container_vector,
                                       master_slave_container_vector);
        BuiltinTimer mpc_add_time;
        BaseType::AddConstraintsToModelpart(BaseType::mrMainModelPart,
                                            master_slave_container_vector);
        KRATOS_INFO_IF(
            "Adding of MPCs from containers to modelpart took         : ",
            BaseType::mEchoLevel > 1)
            << mpc_add_time.ElapsedSeconds() << " seconds" << std::endl;
    }

// Template declarations
template class ApplyChimeraProcessMonolithic<2>;
template class ApplyChimeraProcessMonolithic<3>;

} // namespace Kratos.

