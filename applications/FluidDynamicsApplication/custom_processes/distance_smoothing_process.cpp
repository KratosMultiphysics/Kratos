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
#include "includes/define.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "spaces/ublas_space.h"
#include "modeler/connectivity_preserve_modeler.h"
#include "utilities/variable_utils.h"
#include "includes/global_pointer_variables.h"
#include "utilities/parallel_utilities.h"
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "distance_smoothing_process.h"
#include "custom_elements/distance_smoothing_element.h"


namespace Kratos
{

/* Public functions *******************************************************/
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
DistanceSmoothingProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::DistanceSmoothingProcess(
    ModelPart& rModelPart,
    typename TLinearSolver::Pointer p_linear_solver)
    : Process(),
      mrModelPart(rModelPart),
      mrModel(rModelPart.GetModel()),
      mAuxModelPartName("smoothing_model_part")//mrModelPart.FullName()+"_Smoothing"")
{
    // Generate an auxilary model part and populate it by elements of type DistanceSmoothingElement
    CreateAuxModelPart();

    CreateSolutionStrategy(p_linear_solver);
}

template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
DistanceSmoothingProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::DistanceSmoothingProcess(
    ModelPart& rModelPart,
    Parameters Parameters)
    : DistanceSmoothingProcess(
        rModelPart,
        LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(Parameters["linear_solver_settings"])
    ){}

template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
DistanceSmoothingProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::DistanceSmoothingProcess(
    Model &rModel,
    Parameters Parameters)
    : DistanceSmoothingProcess(
        rModel.GetModelPart(Parameters["model_part_name"].GetString()),
        Parameters
    ){}

template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void DistanceSmoothingProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::CreateAuxModelPart()
{
    if(mrModel.HasModelPart( mAuxModelPartName ))
        mrModel.DeleteModelPart( mAuxModelPartName );

    // Generate AuxModelPart
    ModelPart& r_smoothing_model_part = mrModel.CreateModelPart( mAuxModelPartName );

    Element::Pointer p_smoothing_element = Kratos::make_intrusive<DistanceSmoothingElement<TDim>>();

    ConnectivityPreserveModeler modeler;
    modeler.GenerateModelPart(mrModelPart, r_smoothing_model_part, *p_smoothing_element);

    const unsigned int buffer_size = r_smoothing_model_part.GetBufferSize();
    KRATOS_ERROR_IF(buffer_size < 2) << "Buffer size should be at least 2" << std::endl;

    // Adding DISTANCE to the solution variables is not needed if it is already a solution variable of the problem
    r_smoothing_model_part.AddNodalSolutionStepVariable(DISTANCE);

    // Ensure that the nodes have distance as a DOF
    VariableUtils().AddDof<Variable<double> >(DISTANCE, r_smoothing_model_part);
}

template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void DistanceSmoothingProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::Execute()
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
                if (n_nodes[j].Is(CONTACT) == rNode.Is(CONTACT)){

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

template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void DistanceSmoothingProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::Clear()
{
    ModelPart& r_smoothing_model_part = mrModel.GetModelPart( mAuxModelPartName );
    r_smoothing_model_part.Nodes().clear();
    r_smoothing_model_part.Conditions().clear();
    r_smoothing_model_part.Elements().clear();
    mp_solving_strategy->Clear();
}

template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
void DistanceSmoothingProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::CreateSolutionStrategy(typename TLinearSolver::Pointer pLinearSolver)
{
    // Generate a linear solver strategy
    auto p_scheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();

    ModelPart& r_smoothing_model_part = mrModel.GetModelPart( mAuxModelPartName );

    const bool CalculateReactions = false;
    const bool ReformDofAtEachIteration = false;
    const bool CalculateNormDxFlag = false;

    auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(pLinearSolver);

    mp_solving_strategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(
        r_smoothing_model_part,
        p_scheme,
        p_builder_solver,
        CalculateReactions,
        ReformDofAtEachIteration,
        CalculateNormDxFlag);

    mp_solving_strategy->Initialize();
    mp_solving_strategy->SetEchoLevel(0);
    mp_solving_strategy->Check();
}

typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

template class DistanceSmoothingProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>;
template class DistanceSmoothingProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>;

}; // namespace Kratos
