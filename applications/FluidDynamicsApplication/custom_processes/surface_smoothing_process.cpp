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

// Project includes
#include "utilities/variable_utils.h"
#include "includes/global_pointer_variables.h"
#include "utilities/parallel_utilities.h"
#include "factories/linear_solver_factory.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
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

SurfaceSmoothingProcess::SurfaceSmoothingProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Generate an auxilary model part and populate it by elements of type SurfaceSmoothingElement
    CreateAuxModelPart();

    TLinearSolver::Pointer plinear_solver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(
        rParameters["linear_solver_settings"]);

    auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);

    InitializeSolutionStrategy(plinear_solver, p_builder_solver);
}

SurfaceSmoothingProcess::SurfaceSmoothingProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Generate an auxilary model part and populate it by elements of type SurfaceSmoothingElement
    CreateAuxModelPart();

    TLinearSolver::Pointer plinear_solver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(
        rParameters["linear_solver_settings"]);

    auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);

    InitializeSolutionStrategy(plinear_solver, p_builder_solver);
}

void SurfaceSmoothingProcess::CreateAuxModelPart()
{
    Model& current_model = mrModelPart.GetModel();
    if(current_model.HasModelPart( mAuxModelPartName ))
        current_model.DeleteModelPart( mAuxModelPartName );

    // Generate AuxModelPart
    ModelPart& r_smoothing_model_part = current_model.CreateModelPart( mAuxModelPartName );

    Element::Pointer p_smoothing_element = Kratos::make_intrusive<SurfaceSmoothingElement>();

    ConnectivityPreserveModeler modeler;
    modeler.GenerateModelPart(mrModelPart, r_smoothing_model_part, *p_smoothing_element);

    const unsigned int buffer_size = mrModelPart.GetBufferSize();
    KRATOS_ERROR_IF(buffer_size < 2) << "Buffer size should be at least 2" << std::endl;

    // Adding DISTANCE to the solution variables is not needed if it is already a solution variable of the problem
    r_smoothing_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Ensure that the nodes have distance as a DOF
    VariableUtils().AddDof<Variable<double> >(DISTANCE, r_smoothing_model_part);
}

void SurfaceSmoothingProcess::Execute()
{
    KRATOS_TRY;

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            rNode.Free(DISTANCE);
            const double distance = rNode.FastGetSolutionStepValue(DISTANCE);
            rNode.FastGetSolutionStepValue(DISTANCE, 1) = distance;
        });

    mp_solving_strategy->Solve();

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
            rNode.SetValue( DISTANCE, rNode.FastGetSolutionStepValue(DISTANCE)
                - rNode.FastGetSolutionStepValue(DISTANCE, 1) ); // Corrected distance difference
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

#ifdef KRATOS_DEBUG
                KRATOS_WARNING_IF("DistanceSmoothingProcess", distance_ij < 1.0e-12)
                    << "WARNING: Neighbouring nodes are almost coinciding" << std::endl;
#endif

                if (distance_ij > 1.0e-12){
                    weight += 1.0/distance_ij;
                    dist_diff_avg += n_nodes[j].GetValue(DISTANCE)/distance_ij;
                }
            }

#ifdef KRATOS_DEBUG
            KRATOS_WARNING_IF("DistanceSmoothingProcess", weight < 1.0e-12)
                << "WARNING: Correction is not done due to a zero weight" <<std::endl;
#endif

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

void SurfaceSmoothingProcess::Clear()
{
    Model& r_model = mrModelPart.GetModel();
    ModelPart& r_smoothing_model_part = r_model.GetModelPart( mAuxModelPartName );
    r_smoothing_model_part.Nodes().clear();
    r_smoothing_model_part.Conditions().clear();
    r_smoothing_model_part.Elements().clear();
    mp_solving_strategy->Clear();
}

void SurfaceSmoothingProcess::InitializeSolutionStrategy(
        TLinearSolver::Pointer pLinearSolver,
        BuilderSolverPointerType pBuilderAndSolver)
    {
        // Generate a linear solver strategy
        auto p_scheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();

        Model& r_model = mrModelPart.GetModel();
        ModelPart& r_smoothing_model_part = r_model.GetModelPart( mAuxModelPartName );

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        mp_solving_strategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(
            r_smoothing_model_part,
            p_scheme,
            pLinearSolver,
            pBuilderAndSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        mp_solving_strategy->Check();
    }

}; // namespace Kratos