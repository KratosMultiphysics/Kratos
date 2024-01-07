//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "utilities/sensitivity_utilities.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "fluid_adjoint_slip_utilities.h"

namespace Kratos
{

FluidAdjointSlipUtilities::FluidAdjointSlipUtilities(
    const IndexType Dimension,
    const IndexType BlockSize,
    const Variable<double>& rSlipDofVariable)
    : mDimension(Dimension),
      mBlockSize(BlockSize),
      mRotationTool(Dimension, mBlockSize),
      mrSlipDofVariable(rSlipDofVariable)
{
    KRATOS_TRY

    if (Dimension == 2) {
        this->mAddNodalRotationNonShapeVariableDerivativesMethod = &FluidAdjointSlipUtilities::TemplatedAddNodalRotationNonShapeVariableDerivatives<2>;
        this->mAddNodalRotationShapeVariableDerivativesMethod = &FluidAdjointSlipUtilities::TemplatedAddNodalRotationShapeVariableDerivatives<2>;
        this->mAddNodalApplySlipConditionShapeVariableDerivativesMethod = &FluidAdjointSlipUtilities::TemplatedAddNodalApplySlipConditionShapeVariableDerivatives<2>;
    } else if (Dimension == 3) {
        this->mAddNodalRotationNonShapeVariableDerivativesMethod = &FluidAdjointSlipUtilities::TemplatedAddNodalRotationNonShapeVariableDerivatives<3>;
        this->mAddNodalRotationShapeVariableDerivativesMethod = &FluidAdjointSlipUtilities::TemplatedAddNodalRotationShapeVariableDerivatives<3>;
        this->mAddNodalApplySlipConditionShapeVariableDerivativesMethod = &FluidAdjointSlipUtilities::TemplatedAddNodalApplySlipConditionShapeVariableDerivatives<3>;
    } else {
        KRATOS_ERROR << "Unsupported dimensionality requested. Only 2D and 3D "
                        "supported. [ Dimension = "
                     << Dimension << " ].\n";
    }

    KRATOS_CATCH("");
}

FluidAdjointSlipUtilities::FluidAdjointSlipUtilities(
    const IndexType Dimension,
    const IndexType BlockSize)
    : FluidAdjointSlipUtilities(Dimension, BlockSize, ADJOINT_FLUID_VECTOR_1_X)
{
}

void FluidAdjointSlipUtilities::Initialize(ModelPart& rModelPart)
{
    KRATOS_TRY

    // Find nodal neighbors for conditions
    FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> neighbor_process(
        rModelPart.GetCommunicator().GetDataCommunicator(), rModelPart,
        NEIGHBOUR_CONDITION_NODES);
    neighbor_process.Execute();
    mNodalNeighboursMap = neighbor_process.GetNeighbourIds(rModelPart.Nodes());

    // Assign condition normal derivatives to corresponding nodes
    // NORMAL_SHAPE_DERIVATIVE on conditions should already be available before
    // executing the following command. (This is done in adjoint_vmsmonolithic_solver.py)
    SensitivityUtilities::AssignEntityDerivativesToNodes<ModelPart::ConditionsContainerType>(
        rModelPart, mDimension, NORMAL_SHAPE_DERIVATIVE,
        mNodalNeighboursMap, 1.0 / mDimension, SLIP);

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_nodes = rModelPart.Nodes();
    const int number_of_nodes = r_nodes.size();

    std::vector<int> indices;
    indices.resize(number_of_nodes);
    IndexPartition<int>(number_of_nodes).for_each([&](const int Index) {
        indices[Index] = (r_nodes.begin() + Index)->Id();
    });

    const auto& r_data_communicator = r_communicator.GetDataCommunicator();
    mGlobalPointerNodalMap = GlobalPointerUtilities::RetrieveGlobalIndexedPointersMap(
        r_nodes, indices, r_data_communicator);

    KRATOS_CATCH("");
}

void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const GeometryType& rGeometry) const
{
    if (rOutput.size1() != rResidualDerivatives.size1() ||
        rOutput.size2() != rResidualDerivatives.size2()) {
        rOutput.resize(rResidualDerivatives.size1(), rResidualDerivatives.size2(), false);
    }

    rOutput.clear();

    const IndexType number_of_nodes = rGeometry.PointsNumber();

    // add residual derivative contributions
    for (IndexType a = 0; a < number_of_nodes; ++a) {
        const auto& r_node = rGeometry[a];
        const IndexType block_index = a * mBlockSize;
        if (r_node.Is(SLIP)) {
            AddNodalRotationNonShapeVariableDerivatives(rOutput, rResidualDerivatives, block_index, r_node);
            if (!r_node.IsFixed(mrSlipDofVariable)) {
                // in here we have to check whether SLIP and not INLET
                // At the inlet ADJOINT_FLUID_VECTOR_1_X is fixed.
                // There can be a situation where a node is
                // shared among both inlet and slip (corner node). In this case
                // AddNodalApplySlipConditionSlipVariableDerivatives will set all the values in the
                // column corresponding to ADJOINT_FLUID_VECTOR_1_X of the common node to zero except for the
                // nodal normal derivative contributions. If the normal is perpendicular to the X direction, then
                // this will set the diagonal in the column corresponding to this node's ADJOINT_FLUID_VECTOR_1_X dof
                // to zero. But there can be non-zero off-diagonal due to normal not being zero. Eventhough
                // ADJOINT_FLUID_VECTOR_1_X will be fixed in the block builder and solver
                // it will not recognize this column as empty, so diagonal won't be fixed by B&S to 1.0.
                // Since ADJOINT_FLUID_VECTOR_1_X of this node is fixed, B&S will set all the off-diagonals
                // to zero. This will create a zero column in the system matrix which makes it singular.
                // So this check is performed here to avoid this.
                // In here inlets are alsways assumed to be dirichlet.
                AddNodalApplySlipConditionSlipVariableDerivatives(rOutput, block_index, r_node);
            }
        } else {
            AddNodalResidualDerivatives(rOutput, rResidualDerivatives, block_index);
        }
    }
}

template<class TEntityType>
void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(
    Matrix& rOutput,
    GlobalPointersVector<NodeType>& rGPSensitivityVector,
    const Vector& rResiduals,
    const Matrix& rResidualDerivatives,
    const TEntityType& rEntity,
    const ProcessInfo& rProcessInfo) const
{
        KRATOS_TRY;

        using DofsVectorType = std::vector<Dof<double>::Pointer>;

        const IndexType residuals_size = rResiduals.size();

        KRATOS_ERROR_IF(residuals_size != rResidualDerivatives.size2())
            << "mAdjointVectors.size(): " << residuals_size
            << " incompatible with rResidualDerivatives.size1(): "
            << rResidualDerivatives.size2() << ".\n";

        // a map to store node id as key and the node order as the value
        // dofs are used in here, because if wall distance is calculated in the condition
        // then there will be derivative contributions coming to domain internal node as well
        // hence, list of dofs are used in here.
        DofsVectorType dofs;
        rEntity.GetDofList(dofs, rProcessInfo);

        // get derivative node ids
        std::vector<int> derivative_node_ids;
        for (IndexType i = 0; i < dofs.size(); ++i) {
            if (std::find(derivative_node_ids.begin(), derivative_node_ids.end(), dofs[i]->Id()) == derivative_node_ids.end()) {
                derivative_node_ids.push_back(dofs[i]->Id());
            }
        }

        std::unordered_map<IndexType, IndexType> gp_index_map;
        auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType number_of_derivative_nodes = derivative_node_ids.size();

        if (rGPSensitivityVector.size() != number_of_derivative_nodes) {
            rGPSensitivityVector.resize(number_of_derivative_nodes);
        }

        // add geometry gps
        for (IndexType i = 0; i < number_of_derivative_nodes; ++i) {
            auto& r_gp = mGlobalPointerNodalMap.find(derivative_node_ids[i])->second;
            rGPSensitivityVector(i) = r_gp;
            gp_index_map[derivative_node_ids[i]] = i;
        }

        KRATOS_DEBUG_ERROR_IF(number_of_derivative_nodes * mDimension != rResidualDerivatives.size1())
            << "Entity sensitivity matrix size mismatch. [ rResidualDerivatives.size = ( " << rResidualDerivatives.size1()
            << ", " << rResidualDerivatives.size2() << " ), required size = ( " << (number_of_derivative_nodes * mDimension)
            << ", " << residuals_size << " ) ].\n";

        // add relevant neighbour gps
        IndexType local_index = number_of_derivative_nodes;
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            if (r_node.Is(SLIP)) {
                const auto& neighbour_gps = r_node.GetValue(NEIGHBOUR_CONDITION_NODES);
                const auto& neighbour_ids = mNodalNeighboursMap.find(r_node.Id())->second;
                for (IndexType i = 0; i < neighbour_ids.size(); ++i) {
                    const auto neighbour_id = neighbour_ids[i];
                    const auto p_itr = gp_index_map.find(neighbour_id);
                    if (p_itr == gp_index_map.end()) {
                        rGPSensitivityVector.push_back(neighbour_gps(i));
                        gp_index_map[neighbour_id] = local_index++;
                    }
                }
            }
        }

        const IndexType derivatives_size = local_index * mDimension;
        if (rOutput.size1() != derivatives_size || rOutput.size2() != residuals_size) {
            rOutput.resize(derivatives_size, residuals_size, false);
        }
        rOutput.clear();

        // add residual derivative contributions
        for (IndexType a = 0; a < number_of_nodes; ++a) {
            const auto& r_node = r_geometry[a];
            const IndexType block_index = a * mBlockSize;
            if (r_node.Is(SLIP)) {
                AddNodalRotationShapeVariableDerivatives(rOutput, rResiduals, rResidualDerivatives, block_index, gp_index_map, r_node);
                AddNodalApplySlipConditionShapeVariableDerivatives(rOutput, block_index, gp_index_map, r_node);
            } else {
                AddNodalResidualDerivatives(rOutput, rResidualDerivatives, block_index);
            }
        }

        KRATOS_CATCH("");
}

template<class TEntityType>
void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(
    Matrix& rOutput,
    std::vector<IndexType>& rNodeIds,
    const Vector& rResiduals,
    const Matrix& rResidualDerivatives,
    const TEntityType& rEntity,
    const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY

    GlobalPointersVector<NodeType> r_gp_vector;
    CalculateRotatedSlipConditionAppliedShapeVariableDerivatives<TEntityType>(rOutput, r_gp_vector, rResiduals, rResidualDerivatives, rEntity, rProcessInfo);

    if (rNodeIds.size() != r_gp_vector.size()) {
        rNodeIds.resize(r_gp_vector.size());
    }

    for (IndexType i = 0; i < r_gp_vector.size(); ++i) {
        // This map needs to be changed to a map which includes nodal condition neighbours as well.
        auto it = std::find_if(mGlobalPointerNodalMap.begin(), mGlobalPointerNodalMap.end(), [&](const std::pair<int, GlobalPointer<ModelPart::NodeType>>& p) {
            return GlobalPointerComparor<ModelPart::NodeType>()(p.second, r_gp_vector(i));
        });

        if (it != mGlobalPointerNodalMap.end()) {
            rNodeIds[i] = it->first;
        } else {
            KRATOS_ERROR << "Un-identified node global pointer.\n";
        }
    }

    KRATOS_CATCH("");
}

void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedNonSlipNonShapeVariableDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const GeometryType& rGeometry) const
{
    if (rOutput.size1() != rResidualDerivatives.size1() ||
        rOutput.size2() != rResidualDerivatives.size2()) {
        rOutput.resize(rResidualDerivatives.size1(), rResidualDerivatives.size2(), false);
    }

    rOutput.clear();

    const IndexType number_of_nodes = rGeometry.PointsNumber();

    // add residual derivative contributions
    for (IndexType a = 0; a < number_of_nodes; ++a) {
        const auto& r_node = rGeometry[a];
        const IndexType block_index = a * mBlockSize;
        if (r_node.Is(SLIP)) {
            AddNodalRotationNonShapeVariableDerivatives(rOutput, rResidualDerivatives, block_index, r_node);
            // since slip condition is only based on first derivative
            // variable, make the column zero for all derivatives
            if (!r_node.IsFixed(mrSlipDofVariable)) {
                // in here we have to check whether SLIP and not INLET.
                // At the inlet ADJOINT_FLUID_VECTOR_1_X is fixed.
                // There can be a situation where a node is
                // shared among both inlet and slip (corner node). In this case
                // AddNodalApplySlipConditionSlipVariableDerivatives will set all the values in the
                // column corresponding to ADJOINT_FLUID_VECTOR_1_X of the common node to zero except for the
                // nodal normal derivative contributions. If the normal is perpendicular to the X direction, then
                // this will set the diagonal in the column corresponding to this node's ADJOINT_FLUID_VECTOR_1_X dof
                // to zero. But there can be non-zero off-diagonal due to normal not being zero. Eventhough
                // ADJOINT_FLUID_VECTOR_1_X will be fixed in the block builder and solver
                // it will not recognize this column as empty, so diagonal won't be fixed by B&S to 1.0.
                // Since ADJOINT_FLUID_VECTOR_1_X of this node is fixed, B&S will set all the off-diagonals
                // to zero. This will create a zero column in the system matrix which makes it singular.
                // So this check is performed here to avoid this.
                // In here inlets are alsways assumed to be dirichlet.
                ClearNodalResidualDerivatives(rOutput, block_index);
            }
        } else {
            AddNodalResidualDerivatives(rOutput, rResidualDerivatives, block_index);
        }
    }
}

void FluidAdjointSlipUtilities::AddNodalApplySlipConditionSlipVariableDerivatives(
    Matrix& rOutput,
    const IndexType NodeStartIndex,
    const NodeType& rNode) const
{
    KRATOS_TRY

    // Apply slip condition in primal scheme makes first momentum dof
    // fixed, making the velocity in the normal direction as rhs.

    // first clear the residual derivative
    ClearNodalResidualDerivatives(rOutput, NodeStartIndex);

    auto normal = rNode.FastGetSolutionStepValue(NORMAL);
    normal /= norm_2(normal);

    for (IndexType i = 0; i < mDimension; ++i) {
        rOutput(NodeStartIndex + i, NodeStartIndex) -= normal[i];
    }

    KRATOS_CATCH("");
}

void FluidAdjointSlipUtilities::AddNodalResidualDerivatives(
    Matrix& rOutput,
    const Matrix& rResidualDerivatives,
    const IndexType NodeStartIndex) const
{
    KRATOS_TRY

    // add non-rotated residual derivative contributions
    for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
        for (IndexType i = 0; i < mBlockSize; ++i) {
            rOutput(c, NodeStartIndex + i) += rResidualDerivatives(c, NodeStartIndex + i);
        }
    }

    KRATOS_CATCH("");
}

void FluidAdjointSlipUtilities::ClearNodalResidualDerivatives(
    Matrix& rOutput,
    const IndexType ResidualIndex) const
{
    for (IndexType c = 0; c < rOutput.size1(); ++c) {
        rOutput(c, ResidualIndex) = 0.0;
    }
}

// template instantiations
template void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(Matrix&, GlobalPointersVector<NodeType>&, const Vector&, const Matrix&, const ElementType&, const ProcessInfo&) const;
template void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(Matrix&, GlobalPointersVector<NodeType>&, const Vector&, const Matrix&, const ConditionType&, const ProcessInfo&) const;

template void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(Matrix&, std::vector<IndexType>&, const Vector&, const Matrix&, const ElementType&, const ProcessInfo&) const;
template void FluidAdjointSlipUtilities::CalculateRotatedSlipConditionAppliedShapeVariableDerivatives(Matrix&, std::vector<IndexType>&, const Vector&, const Matrix&, const ConditionType&, const ProcessInfo&) const;

} // namespace Kratos