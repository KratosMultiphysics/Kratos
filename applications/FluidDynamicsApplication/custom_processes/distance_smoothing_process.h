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

#ifndef KRATOS_DISTANCE_SMOOTHING_H
#define KRATOS_DISTANCE_SMOOTHING_H

// System includes
#include <string>
#include <iostream>

// Project includes
#include "processes/process.h"
#include "processes/find_elements_neighbours_process.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/find_global_nodal_elemental_neighbours_process.h"
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/global_pointer_variables.h"
#include "containers/model.h"
#include "modeler/connectivity_preserve_modeler.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/pointer_map_communicator.h"
#include "spaces/ublas_space.h"
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/distance_smoothing_element.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Utility for distance smoothing
/// Based on Tornberg, Anna-Karin, and Björn Engquist. "A finite element based level-set method
/// for multiphase flow applications." Computing and Visualization in Science 3, no. 1-2 (2000): 93-101.
/// The algorithm is improved by imposing a boundary condition and a correction step.
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) DistanceSmoothingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceSmoothingProcess
    KRATOS_CLASS_POINTER_DEFINITION(DistanceSmoothingProcess);

    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef typename SchemeType::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverPointerType;
    typedef SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DistanceSmoothingProcess(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer p_linear_solver,
        const bool AllConditionsAsBoundary = true)
        : Process(),
        mrModelPart(rModelPart),
        mrModel(rModelPart.GetModel()),
        mAuxModelPartName("smoothing_model_part"), //mrModelPart.FullName()+"_Smoothing"")
        mAllConditionsAsBoundary(AllConditionsAsBoundary),
        mAuxModelPartIsInitialized(false)
    {
        // Generate an auxilary model part and populate it by elements of type DistanceSmoothingElement
        CreateAuxModelPart();

        auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(p_linear_solver);

        CreateSolutionStrategy(p_builder_solver);
    }

    /// Constructor.
    DistanceSmoothingProcess(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer p_linear_solver,
        BuilderSolverPointerType pBuilderAndSolver,
        const bool AllConditionsAsBoundary = true)
        : Process(),
        mrModelPart(rModelPart),
        mrModel(rModelPart.GetModel()),
        mAuxModelPartName("smoothing_model_part"),
        mAllConditionsAsBoundary(AllConditionsAsBoundary),
        mAuxModelPartIsInitialized(false)
    {
        // Generate an auxilary model part and populate it by elements of type DistanceSmoothingElement
        CreateAuxModelPart();

        CreateSolutionStrategy(pBuilderAndSolver);
    }

    /// Constructor with Kratos parameters.
    DistanceSmoothingProcess(
        ModelPart& rModelPart,
        Parameters Parameters)
        : DistanceSmoothingProcess(
        rModelPart,
        LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(Parameters["linear_solver_settings"])
    ){}

    /// Constructor with Kratos model
    DistanceSmoothingProcess(
        Model& rModel,
        Parameters Parameters)
        : DistanceSmoothingProcess(
        rModel.GetModelPart(Parameters["model_part_name"].GetString()),
        Parameters
    ){}

    /// Destructor.
    ~DistanceSmoothingProcess() override
    {
        Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    void Execute() override
    {
        KRATOS_TRY;

        if(mAuxModelPartIsInitialized == false){
            CreateAuxModelPart();
        }

        auto& r_smoothing_model_part = mrModel.GetModelPart( mAuxModelPartName );

        block_for_each(r_smoothing_model_part.Nodes(), [&](Node<3>& rNode){
                rNode.Free(DISTANCE);
                const double distance = rNode.FastGetSolutionStepValue(DISTANCE);
                rNode.FastGetSolutionStepValue(DISTANCE, 1) = distance;
            });

        mp_solving_strategy->Solve();

        auto &r_communicator = r_smoothing_model_part.GetCommunicator();
        r_communicator.SynchronizeCurrentDataToMin(DISTANCE);

        block_for_each(r_smoothing_model_part.Nodes(), [&](Node<3>& rNode){
                rNode.SetValue( DISTANCE, rNode.FastGetSolutionStepValue(DISTANCE)
                    - rNode.FastGetSolutionStepValue(DISTANCE, 1) ); // Corrected distance difference
            });

        auto& r_data_comm = r_smoothing_model_part.GetCommunicator().GetDataCommunicator();
        GlobalPointersVector< Node<3 > > gp_list;

        for (int i_node = 0; i_node < static_cast<int>(r_smoothing_model_part.NumberOfNodes()); ++i_node){
            auto it_node = r_smoothing_model_part.NodesBegin() + i_node;
            GlobalPointersVector< Node<3 > >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {
                auto& global_pointer = global_pointer_list(j);
                gp_list.push_back(global_pointer);
            }
        }

        GlobalPointerCommunicator< Node<3 > > pointer_comm(r_data_comm, gp_list);

        auto coordinate_proxy = pointer_comm.Apply(
            [](GlobalPointer<Node<3> >& global_pointer) -> Point::CoordinatesArrayType
            {
                return global_pointer->Coordinates();
            }
        );

        auto distance_proxy = pointer_comm.Apply(
            [&](GlobalPointer<Node<3> >& global_pointer) -> double
            {
                return global_pointer->GetValue(DISTANCE);
            }
        );

        auto contact_proxy = pointer_comm.Apply(
            [&](GlobalPointer<Node<3> >& global_pointer) -> bool
            {
                return global_pointer->Is(CONTACT);
            }
        );

        r_communicator.GetDataCommunicator().Barrier();

        block_for_each(r_smoothing_model_part.Nodes(), [&](Node<3>& rNode){
            const auto& x_i = rNode.Coordinates();

            double weight = 0.0;
            double dist_diff_avg = 0.0;

            GlobalPointersVector< Node<3 > >& global_pointer_list = rNode.GetValue(NEIGHBOUR_NODES);

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {
                auto& global_pointer = global_pointer_list(j);
                const auto x_j = coordinate_proxy.Get(global_pointer);

                if (contact_proxy.Get(global_pointer) == rNode.Is(CONTACT)){

                    const Vector dx = x_i - x_j;

                    const double distance_ij = norm_2(dx);

#ifdef KRATOS_DEBUG
                    KRATOS_WARNING_IF("DistanceSmoothingProcess", distance_ij < 1.0e-12)
                        << "WARNING: Neighbouring nodes are almost coinciding" << std::endl;
#endif

                    if (distance_ij > 1.0e-12){
                        weight += 1.0/distance_ij;
                        dist_diff_avg += distance_proxy.Get(global_pointer)/distance_ij;
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

        r_communicator.SynchronizeCurrentDataToMin(DISTANCE);

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        if(mrModel.HasModelPart( mAuxModelPartName ))
            mrModel.DeleteModelPart( mAuxModelPartName );
        mAuxModelPartIsInitialized = false;

        mp_solving_strategy->Clear();
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Distance Smoothing Process" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "Distance Smoothing Process";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Model& mrModel;
    std::string mAuxModelPartName;
    bool mAllConditionsAsBoundary;
    bool mAuxModelPartIsInitialized;

    typename SolvingStrategyType::UniquePointer mp_solving_strategy;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Create a Solution Strategy object
     * This method creates the linear solution strategy
     * @param pBuilderAndSolver Builder and solver pointer
     */
    void CreateSolutionStrategy(BuilderSolverPointerType pBuilderAndSolver)
    {
        // Generate a linear solver strategy
        auto p_scheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();

        ModelPart& r_smoothing_model_part = mrModel.GetModelPart( mAuxModelPartName );

        const bool CalculateReactions = false;
        const bool ReformDofAtEachIteration = false;
        const bool CalculateNormDxFlag = false;

        mp_solving_strategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(
            r_smoothing_model_part,
            p_scheme,
            pBuilderAndSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        mp_solving_strategy->Initialize();
        mp_solving_strategy->SetEchoLevel(0);
        mp_solving_strategy->Check();
    }

    /**
     * @brief Initialize the process
     * Create a Model Part with the DistanceSmoothingElement
     */
    void CreateAuxModelPart()
    {
        KRATOS_TRY

        if(mrModel.HasModelPart( mAuxModelPartName ))
            mrModel.DeleteModelPart( mAuxModelPartName );

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof<Variable<double> >(DISTANCE, mrModelPart);

        // Generate AuxModelPart
        ModelPart& r_smoothing_model_part = mrModel.CreateModelPart( mAuxModelPartName );

        Element::Pointer p_smoothing_element = Kratos::make_intrusive<DistanceSmoothingElement<TDim>>();

        r_smoothing_model_part.GetNodalSolutionStepVariablesList() = mrModelPart.GetNodalSolutionStepVariablesList();

        ConnectivityPreserveModeler modeler;
        modeler.GenerateModelPart(mrModelPart, r_smoothing_model_part, *p_smoothing_element);

        const unsigned int buffer_size = r_smoothing_model_part.GetBufferSize();
        KRATOS_ERROR_IF(buffer_size < 2) << "Buffer size should be at least 2" << std::endl;

        FindGlobalNodalNeighboursProcess nodal_neighbour_process_new(r_smoothing_model_part);
        nodal_neighbour_process_new.Execute();

        const unsigned int num_dim = TDim;
        const unsigned int num_neighbouring_elements = num_dim + 1;

        FindElementalNeighboursProcess neighbour_elements_finder_new(r_smoothing_model_part, num_dim, num_neighbouring_elements);
        neighbour_elements_finder_new.Execute();

        if (mAllConditionsAsBoundary)
        {
            block_for_each(r_smoothing_model_part.Conditions(), [&](Condition& rCondition){
                rCondition.Set(CONTACT, true);
            });
        }

        mAuxModelPartIsInitialized = true;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    DistanceSmoothingProcess() = delete;

    /// Assignment operator.
    DistanceSmoothingProcess& operator=(DistanceSmoothingProcess const& rOther) = delete;

    /// Copy constructor.
    DistanceSmoothingProcess(DistanceSmoothingProcess const& rOther) = delete;

    ///@}

}; // Class DistanceSmoothingProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_DISTANCE_SMOOTHING__H
