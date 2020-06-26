//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

// System includes

// External includes

// Project includes

// Application includes
#include "surface_smoothing_process.h"

namespace Kratos
{

/* Public functions *******************************************************/
SurfaceSmoothingProcess::SurfaceSmoothingProcess(
    ModelPart& rModelPart,
    TLinearSolver::Pointer plinear_solver)
    : Process(),
      mrModelPart(rModelPart)
{
    // Generate an auxilary model part and populate it by elements of type SurfaceSmoothingElement
    CreateAuxModelPart();

    auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);

    InitializeSolutionStrategy(plinear_solver, p_builder_solver);
}

void SurfaceSmoothingProcess::CreateAuxModelPart()
{
    Model& current_model = mrModelPart.GetModel();
    if(current_model.HasModelPart( mAuxModelPartName ))
        current_model.DeleteModelPart( mAuxModelPartName );

    // Adding DISTANCE to the solution variables is not needed if it is already a solution variable of the problem
    mrModelPart.AddNodalSolutionStepVariable(DISTANCE);

    // Ensure that the nodes have distance as a DOF
    VariableUtils().AddDof<Variable<double> >(DISTANCE, mrModelPart);

    // Generate AuxModelPart
    ModelPart& r_smoothing_model_part = current_model.CreateModelPart( mAuxModelPartName );

    Element::Pointer p_smoothing_element = Kratos::make_intrusive<SurfaceSmoothingElement>();

    ConnectivityPreserveModeler modeler;
    modeler.GenerateModelPart(mrModelPart, r_smoothing_model_part, *p_smoothing_element);

    const double delta_time = mrModelPart.pGetProcessInfo()->GetValue(DELTA_TIME);
    r_smoothing_model_part.pGetProcessInfo()->SetValue(DELTA_TIME, delta_time);
}

void SurfaceSmoothingProcess::Execute()
{
    KRATOS_TRY;

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            rNode.Free(DISTANCE);
            const double distance = rNode.FastGetSolutionStepValue(DISTANCE);
            rNode.SetValue(DISTANCE, distance);
        });

    mp_solving_strategy->Solve();

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            rNode.SetValue( DISTANCE, rNode.FastGetSolutionStepValue(DISTANCE)
                - rNode.GetValue(DISTANCE) ); // Corrected distance difference
        });

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            const double x_i = rNode.X();
            const double y_i = rNode.Y();
            const double z_i = rNode.Z();

            double weight = 0.0;
            double dist_diff_avg = 0.0;

            auto& n_nodes = rNode.GetValue(NEIGHBOUR_NODES);
            for (unsigned int j = 0; j < n_nodes.size(); ++j) {
                const double dx = x_i - n_nodes[j].X();
                const double dy = y_i - n_nodes[j].Y();
                const double dz = z_i - n_nodes[j].Z();
                const double distance_ij = sqrt( dx*dx + dy*dy + dz*dz );

                if (distance_ij > 1.0e-12){
                    weight += 1.0/distance_ij;
                    dist_diff_avg += n_nodes[j].GetValue(DISTANCE)/distance_ij;
                }
            }

            if (weight > 1.0e-12)
                rNode.FastGetSolutionStepValue(DISTANCE) -= dist_diff_avg/weight;
        });

    KRATOS_CATCH("");
}

void SurfaceSmoothingProcess::ExecuteInitialize()
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void SurfaceSmoothingProcess::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void SurfaceSmoothingProcess::ExecuteInitializeSolutionStep() {
}

void SurfaceSmoothingProcess::ExecuteFinalizeSolutionStep() {
}

}; // namespace Kratos