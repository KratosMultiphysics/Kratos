//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Marc Nunez
//

// Project includes
#include "apply_far_field_process.h"
#include "utilities/openmp_utils.h"

namespace Kratos {

// Constructor for ApplyFarFieldProcess Process
ApplyFarFieldProcess::ApplyFarFieldProcess(ModelPart& rModelPart,
                                           const double ReferencePotential,
                                           const bool InitializeFlowField,
                                           const bool PerturbationField)
    : Process(),
      mrModelPart(rModelPart),
      mReferencePotential(ReferencePotential),
      mInitializeFlowField(InitializeFlowField),
      mPerturbationField(PerturbationField),
      mFreeStreamVelocity(mrModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY]) {
}

void ApplyFarFieldProcess::Execute()
{
    FindFarthestUpstreamBoundaryNode();
    AssignFarFieldBoundaryConditions();
    if (mInitializeFlowField){
        InitializeFlowField();
    }
}

void ApplyFarFieldProcess::FindFarthestUpstreamBoundaryNode()
{
    // Declaring omp variables, generating vectors of size = num_threads
    std::size_t num_threads = OpenMPUtils::GetNumThreads();
    std::vector<double> min_projections(num_threads, std::numeric_limits<double>::max());
    std::vector<std::size_t> nodes_id_list(num_threads, 0);

    // Find minimum node in the direction of the free stream velocity in openmp.
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = mrModelPart.NodesBegin() + i;
        std::size_t thread_id = OpenMPUtils::ThisThread();

        const auto& r_coordinates = it_node->Coordinates();
        const double distance_projection = inner_prod(r_coordinates, mFreeStreamVelocity);

        if (distance_projection < min_projections[thread_id]){
            min_projections[thread_id] = distance_projection;
            nodes_id_list[thread_id] = it_node->Id();
        }
    }

    // Find minimum across all threads
    std::size_t minimum_node_thread_id = 0;
    for (std::size_t i_thread = 0; i_thread<num_threads; i_thread++){
        if (min_projections[i_thread] < min_projections[minimum_node_thread_id]){
            minimum_node_thread_id = i_thread;
        }
    }

    // Storing final result in member variable
    mpReferenceNode = mrModelPart.pGetNode(nodes_id_list[minimum_node_thread_id]);
}

void ApplyFarFieldProcess::AssignFarFieldBoundaryConditions()
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Conditions().size()); i++) {
        auto it_cond = mrModelPart.ConditionsBegin() + i;
        auto& r_geometry = it_cond->GetGeometry();

        // Computing normal
        array_1d<double,3> aux_coordinates;
        r_geometry.PointLocalCoordinates(aux_coordinates, r_geometry.Center());
        const auto normal = it_cond->GetGeometry().Normal(aux_coordinates);

        const double velocity_projection = inner_prod(normal, mFreeStreamVelocity);
        if (velocity_projection < 0.0) {
            AssignDirichletFarFieldBoundaryCondition(r_geometry);
        }
        else {
            AssignNeumannFarFieldBoundaryCondition(*it_cond);
        }
    }
}

void ApplyFarFieldProcess::AssignDirichletFarFieldBoundaryCondition(Geometry<NodeType>& rGeometry)
{
    // Fixing nodes in the domain that are part of the inlet, and initializing its value
    // according to the free stream velocity.
    for (std::size_t i_node = 0; i_node < rGeometry.size(); i_node++){
        double inlet_potential = 0.0;

        if(!mPerturbationField){
            array_1d<double,3> relative_coordinates = rGeometry[i_node].Coordinates() - mpReferenceNode->Coordinates();
            inlet_potential = inner_prod(relative_coordinates, mFreeStreamVelocity);
        }

        // Checking if its the primal or the adjoint case.
        if (!rGeometry[i_node].SolutionStepsDataHas(ADJOINT_VELOCITY_POTENTIAL)) {
            rGeometry[i_node].SetLock();
            rGeometry[i_node].Fix(VELOCITY_POTENTIAL);
            rGeometry[i_node].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = inlet_potential + mReferencePotential;
            rGeometry[i_node].UnSetLock();
        }
        else {
            rGeometry[i_node].SetLock();
            rGeometry[i_node].Fix(ADJOINT_VELOCITY_POTENTIAL);
            rGeometry[i_node].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL) = 0.0;
            rGeometry[i_node].UnSetLock();
        }
    }
}

void ApplyFarFieldProcess::AssignNeumannFarFieldBoundaryCondition(Condition& rCondition)
{
    rCondition.SetValue(FREE_STREAM_VELOCITY, mFreeStreamVelocity);
}

void ApplyFarFieldProcess::InitializeFlowField()
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(root_model_part.Nodes().size()); i++) {
        auto it_node = root_model_part.NodesBegin() + i;

        const auto relative_coordinates = it_node->Coordinates() - mpReferenceNode->Coordinates();
        const double inlet_potential = inner_prod(relative_coordinates, mFreeStreamVelocity);

        it_node->FastGetSolutionStepValue(VELOCITY_POTENTIAL) = inlet_potential + mReferencePotential;
        it_node->FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = inlet_potential + mReferencePotential;
    }
}

} // namespace Kratos.
