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

#ifndef KRATOS_CALCULATE_ELECTRIC_FORCE_PROCESS_H
#define KRATOS_CALCULATE_ELECTRIC_FORCE_PROCESS_H

// System includes
#include <string>
#include <iostream>

// Project includes
#include "processes/process.h"
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
#include "droplet_dynamics_application_variables.h"
#include "custom_elements/electrostatic_element.h"


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

template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class KRATOS_API(DROPLET_DYNAMICS_APPLICATION) CalculateElectricForceProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateElectricForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateElectricForceProcess);

    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef typename SchemeType::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverPointerType;
    typedef ImplicitSolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    // Set AllConditionsAsBoundary = false if distance gradient should be imposed only on the conditiones priorly marked as CONTACT.
    CalculateElectricForceProcess(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        const bool AllConditionsAsBoundary = true)
        : Process(),
        mrModelPart(rModelPart),
        mrModel(rModelPart.GetModel()),
        mAuxModelPartName("electric_model_part"),
        mAllConditionsAsBoundary(AllConditionsAsBoundary),
        mAuxModelPartIsInitialized(false)
    {
        // Generate an auxilary model part and populate it by elements of type DistanceSmoothingElement
        CreateAuxModelPart();

        auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(pLinearSolver);

        CreateSolutionStrategy(p_builder_solver);
    }

    /// Constructor.
    CalculateElectricForceProcess(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        BuilderSolverPointerType pBuilderAndSolver,
        const bool AllConditionsAsBoundary = true)
        : Process(),
        mrModelPart(rModelPart),
        mrModel(rModelPart.GetModel()),
        mAuxModelPartName("electric_model_part"),
        mAllConditionsAsBoundary(AllConditionsAsBoundary),
        mAuxModelPartIsInitialized(false)
    {
        // Generate an auxilary model part and populate it by elements of type DistanceSmoothingElement
        CreateAuxModelPart();

        CreateSolutionStrategy(pBuilderAndSolver);
    }

    /// Constructor with Kratos parameters.
    CalculateElectricForceProcess(
        ModelPart& rModelPart,
        Parameters Parameters)
        : CalculateElectricForceProcess(
        rModelPart,
        LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(Parameters["linear_solver_settings"]),
        Parameters["all_conditions_as_boundaries"].GetBool()
    ){}

    /// Constructor with Kratos model
    CalculateElectricForceProcess(
        Model& rModel,
        Parameters Parameters)
        : CalculateElectricForceProcess(
        rModel.GetModelPart(Parameters["model_part_name"].GetString()),
        Parameters
    ){}

    /// Destructor.
    ~CalculateElectricForceProcess() override
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

        mp_solving_strategy->Solve();
        
        auto& process_info = mrModelPart.GetProcessInfo();
        
        std::vector<array_1d<double,3>> n_Fe = { ZeroVector(3) };

        auto& r_electric_model_part = mrModel.GetModelPart( mAuxModelPartName );

        for (auto& p_e_element : r_electric_model_part.Elements()){
            p_e_element.CalculateOnIntegrationPoints(EFORCE, n_Fe , process_info);
        }

        
               
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
        buffer << "Calculate Electric Force Process" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "Calculate Electric Force Process";}

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

        auto& r_electric_model_part = mrModel.GetModelPart( mAuxModelPartName );

        const bool calculate_reactions = false;
        const bool reform_dof_at_each_iteration = false;
        const bool calculate_norm_dx_flag = false;

        mp_solving_strategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(
            r_electric_model_part,
            p_scheme,
            pBuilderAndSolver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_dx_flag);

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

        // Ensure that the nodes have EPOTENTIAL as a DOF
        VariableUtils().AddDof<Variable<double> >(EPOTENTIAL, mrModelPart);

        // Generate AuxModelPart
        auto& r_electric_model_part = mrModel.CreateModelPart( mAuxModelPartName );

        auto p_electric_element = Kratos::make_intrusive<ElectrostaticElement<TDim>>();

        r_electric_model_part.GetNodalSolutionStepVariablesList() = mrModelPart.GetNodalSolutionStepVariablesList();
                
        ConnectivityPreserveModeler modeler;
        modeler.GenerateModelPart(mrModelPart, r_electric_model_part, *p_electric_element);

        const unsigned int buffer_size = r_electric_model_part.GetBufferSize();
        KRATOS_ERROR_IF(buffer_size < 2) << "Buffer size should be at least 2" << std::endl;
               
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
    CalculateElectricForceProcess() = delete;

    /// Assignment operator.
    CalculateElectricForceProcess& operator=(CalculateElectricForceProcess const& rOther) = delete;

    /// Copy constructor.
    CalculateElectricForceProcess(CalculateElectricForceProcess const& rOther) = delete;

    ///@}

}; // Class CalculateElectricForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_CALCULATE_ELECTRIC_FORCE_PROCESS_H
