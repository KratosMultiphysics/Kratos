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
//

#ifndef KRATOS_ED_REINITIALIZATION_PROCESS_H
#define KRATOS_ED_REINITIALIZATION_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/generic_find_elements_neighbours_process.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/find_global_nodal_elemental_neighbours_process.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"
#include "spaces/ublas_space.h"
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "modeler/connectivity_preserve_modeler.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "processes/compute_nodal_gradient_process.h"
#include "utilities/variable_utils.h" //Now necessary!
//#include "processes/compute_nodal_normal_divergence_process.h" //Not needed, already done in python
#include "utilities/element_size_calculator.h"
#include "includes/deprecated_variables.h" //For IS_STRUCTURED

// Application includes
// #include "fluid_dynamics_application_variables.h"
#include "droplet_dynamics_application_variables.h"
#include "custom_elements/elliptic_distance_reinitialization_element.h"

namespace Kratos
{
///@addtogroup DropletDynamicsApplication
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

/// Utility for level-set distance reinitialization
/// There is no need to keep distance gradient unity.Initialize
/// This process tries to keep distance zero at the interface while keeping the curvature is also penalized.
class EDReinitializationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EllipticDistanceReinitialization
    KRATOS_CLASS_POINTER_DEFINITION(EDReinitializationProcess);

    typedef UblasSpace<double, CompressedMatrix, Vector> TSparseSpace;
    typedef UblasSpace<double, Matrix, Vector> TDenseSpace;
    typedef LinearSolver<TSparseSpace, TDenseSpace > TLinearSolver;
    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef typename SchemeType::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverPointerType;
    //typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::UniquePointer BuilderSolverPointerType;

    typedef SolvingStrategy< TSparseSpace, TDenseSpace > SolvingStrategyType;

    typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable> ComputeGradientProcessType;
    typedef ComputeGradientProcessType::Pointer ComputeGradientProcessPointerType;

    // typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>    BlockBuilderAndSolverType;
    // typedef typename BlockBuilderAndSolverType::Pointer                                              BlockBuilderAndSolverPointerType; 


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    // Set AllConditionsAsBoundary = false if distance gradient should be imposed only on the conditiones priorly marked as CONTACT.
    EDReinitializationProcess(
        ModelPart& rModelPart,
        TLinearSolver::Pointer);

    /// Destructor.
    ~EDReinitializationProcess() override
    {
        Model& current_model = mrModelPart.GetModel();
        if(current_model.HasModelPart( mAuxModelPartName ))
            current_model.DeleteModelPart( mAuxModelPartName );
    }

    ///@}
    ///@name Operators
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    virtual void Clear()
    {
        Model& r_model = mrModelPart.GetModel();
        ModelPart& r_edr_model_part = r_model.GetModelPart( mAuxModelPartName );
        r_edr_model_part.Nodes().clear();
        r_edr_model_part.Conditions().clear();
        r_edr_model_part.Elements().clear();
        mp_solving_strategy->Clear();
        mpGradientCalculator->Clear();
    }

    ///@}
    ///@name Public Member Variables
    ///@{
    unsigned int Dim;// = TDim;

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
        buffer << "TESTING: EllipticDistanceReinitialization" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "EllipticDistanceReinitialization";}

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

    std::string mAuxModelPartName = "Aux_EllipticDistanceReinitialization_Model_Part";
    bool mAuxModelPartIsCreated = false;

    SolvingStrategyType::UniquePointer mp_solving_strategy;

    ComputeGradientProcessPointerType mpGradientCalculator = nullptr;
    BuilderSolverPointerType mpBlockBuilderSolver = nullptr;
    //TLinearSolver::Pointer mpLinearSolver = nullptr;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    //void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    void InitializeSolutionStrategy(
        /* TLinearSolver::Pointer pLinearSolver, */
        BuilderSolverPointerType pBuilderAndSolver)
    {
        // Generate a linear solver strategy
        auto p_scheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();

        Model& r_model = mrModelPart.GetModel();
        auto& r_edr_model_part = r_model.GetModelPart( mAuxModelPartName );

        const bool calculate_reactions = false;
        const bool reform_dof_at_each_iteration = false;
        const bool calculate_norm_dx_flag = false;

        mp_solving_strategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(
            r_edr_model_part,
            p_scheme,            
            pBuilderAndSolver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_dx_flag);

    //     mp_solving_strategy->Initialize();
    //     mp_solving_strategy->SetEchoLevel(0);
        mp_solving_strategy->Check();
    }

    /**
     * @brief Initialize the process
     * Create a Model Part with the EllipticDistanceReinitializationElement
     */
    void CreateAuxModelPart();
    // {
    //     KRATOS_TRY

    //     if(mrModel.HasModelPart( mAuxModelPartName ))
    //         mrModel.DeleteModelPart( mAuxModelPartName );

    //     // Ensure that the nodes have DISTANCE_AUX2 as a DOF
    //     VariableUtils().AddDof<Variable<double> >(DISTANCE_AUX2, mrModelPart);

    //     // Generate AuxModelPart
    //     auto& r_edr_model_part = mrModel.CreateModelPart( mAuxModelPartName );

    //     auto p_edr_element = Kratos::make_intrusive<EDReinitializationElement>();//<TDim>

    //     r_edr_model_part.GetNodalSolutionStepVariablesList() = mrModelPart.GetNodalSolutionStepVariablesList();
                
    //     ConnectivityPreserveModeler modeler;
    //     modeler.GenerateModelPart(mrModelPart, r_edr_model_part, *p_edr_element);

    //     const unsigned int buffer_size = r_edr_model_part.GetBufferSize();
    //     KRATOS_ERROR_IF(buffer_size < 2) << "Buffer size should be at least 2" << std::endl;
               
    //     mAuxModelPartIsInitialized = true;

    //     KRATOS_CATCH("")
    // }


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
    EDReinitializationProcess() = delete;

    /// Assignment operator.
    EDReinitializationProcess& operator=(EDReinitializationProcess const& rOther) = delete;

    /// Copy constructor.
    EDReinitializationProcess(EDReinitializationProcess const& rOther) = delete;

    ///@}

}; // Class EDReinitializationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_ED_REINITIALIZATION_PROCESS_H
