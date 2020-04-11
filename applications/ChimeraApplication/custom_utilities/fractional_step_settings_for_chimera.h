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


#ifndef KRATOS_FRACTIONAL_STEP_SETTINGS_FOR_CHIMERA_H
#define KRATOS_FRACTIONAL_STEP_SETTINGS_FOR_CHIMERA_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "processes/process.h"
#include "processes/fast_transfer_between_model_parts_process.h"
#include "containers/model.h"

// Application includes
#include "custom_utilities/solver_settings.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "custom_strategies/custom_builder_and_solver/residualbased_block_builder_and_solver_with_constraints_for_chimera.h"


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

/// Helper class to define solution strategies for FS_Strategy.
template< class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
class KRATOS_API(CHIMERA_APPLICATION) FractionalStepSettingsForChimera: public SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FractionalStepSettingsForChimera
    KRATOS_CLASS_POINTER_DEFINITION(FractionalStepSettingsForChimera);

    typedef SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
    typedef typename BaseType::StrategyType StrategyType;
    typedef typename BaseType::StrategyPointerType StrategyPointerType;
    typedef typename BaseType::ProcessPointerType ProcessPointerType;
    typedef typename BaseType::StrategyLabel StrategyLabel;
    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ResidualBasedBlockBuilderAndSolverType;
    typedef typename ResidualBasedBlockBuilderAndSolverType::Pointer ResidualBasedBlockBuilderAndSolverPointerType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FractionalStepSettingsForChimera(ModelPart& r_model_part,
                   const std::size_t ThisDomainSize,
                   const std::size_t ThisTimeOrder,
                   const bool use_slip,
                   const bool MoveMeshFlag,
                   const bool reform_dof_set):
        BaseType(r_model_part,ThisDomainSize,ThisTimeOrder,use_slip,MoveMeshFlag,reform_dof_set)
    {}

    /// Destructor.
    ~FractionalStepSettingsForChimera() = default;


    /// Default constructor.
    FractionalStepSettingsForChimera() = delete;

    /// Assignment operator.
    FractionalStepSettingsForChimera& operator=(FractionalStepSettingsForChimera const& rOther) = delete;

    /// Copy constructor.
    FractionalStepSettingsForChimera(FractionalStepSettingsForChimera const& rOther) = delete;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void SetStrategy(StrategyLabel const& rStrategyLabel,
                             typename TLinearSolver::Pointer pLinearSolver,
                             const double Tolerance,
                             const unsigned int MaxIter) override
    {
        KRATOS_TRY;

        // pointer types for solution strategy construction
        typedef typename Scheme< TSparseSpace, TDenseSpace >::Pointer SchemePointerType;
        typedef ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera<TSparseSpace, TDenseSpace, TLinearSolver > ResidualBasedBlockBuilderAndSolverWithConstraintsForChimeraType;

        // Default, fixed flags
        bool calculate_reactions = false;
        bool calculate_norm_dx_flag = true;

        ModelPart& r_model_part = BaseType::GetModelPart();

        bool use_slip = BaseType::UseSlipConditions();
        // Modification of the DofSet is managed by the fractional step strategy, not the auxiliary velocity and pressure strategies.
        bool reform_dof_set = BaseType::GetReformDofSet();
        std::size_t echo_level = BaseType::GetEchoLevel();
        std::size_t strategy_echo_level = (echo_level > 0) ? (echo_level-1) : 0;

        if ( rStrategyLabel == BaseType::Velocity )
        {
            std::string fs_vel_mp_name = r_model_part.Name()+"fs_velocity_model_part";
            ModelPart &r_fs_velocity_model_part = r_model_part.CreateSubModelPart(fs_vel_mp_name);
            FastTransferBetweenModelPartsProcess(r_fs_velocity_model_part, r_model_part).Execute();
            // Velocity Builder and Solver
            ResidualBasedBlockBuilderAndSolverPointerType p_build_and_solver =  Kratos::make_shared<ResidualBasedBlockBuilderAndSolverWithConstraintsForChimeraType>(pLinearSolver);

            SchemePointerType p_scheme;
            //initializing fractional velocity solution step
            if (use_slip)
            {
                std::size_t domain_size = BaseType::GetDomainSize();
                SchemePointerType p_temp = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticSchemeSlip< TSparseSpace, TDenseSpace>> (domain_size,domain_size);
                p_scheme.swap(p_temp);
            }
            else
            {
                SchemePointerType p_temp = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace >>();
                p_scheme.swap(p_temp);
            }

            // Strategy
            BaseType::mStrategies[rStrategyLabel] = Kratos::make_shared< ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver >>
                                                                        (r_fs_velocity_model_part, p_scheme, pLinearSolver, p_build_and_solver, calculate_reactions, reform_dof_set, calculate_norm_dx_flag);
        }
        else if ( rStrategyLabel == BaseType::Pressure )
        {
            std::string fs_pressure_mp_name = r_model_part.Name()+"fs_pressure_model_part";
            ModelPart &r_fs_pressure_model_part = r_model_part.CreateSubModelPart(fs_pressure_mp_name);
            FastTransferBetweenModelPartsProcess(r_fs_pressure_model_part, r_model_part).Execute();
            // Pressure Builder and Solver
            ResidualBasedBlockBuilderAndSolverPointerType p_build_and_solver =  Kratos::make_shared<ResidualBasedBlockBuilderAndSolverWithConstraintsForChimeraType>(pLinearSolver);
            SchemePointerType p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace >> ();

            // Strategy
            BaseType::mStrategies[rStrategyLabel] = Kratos::make_shared<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver >>
                                                                        (r_fs_pressure_model_part, p_scheme, pLinearSolver, p_build_and_solver, calculate_reactions, reform_dof_set, calculate_norm_dx_flag);
        }
        else
        {
            KRATOS_ERROR << "Error in FractionalStepSettingsForChimera: Unknown strategy label."<<std::endl;
        }

        BaseType::mTolerances[rStrategyLabel] = Tolerance;

        BaseType::mMaxIter[rStrategyLabel] = MaxIter;

        BaseType::mStrategies[rStrategyLabel]->SetEchoLevel(strategy_echo_level);

        KRATOS_CATCH("");
    }


    bool FindStrategy(StrategyLabel const& rStrategyLabel,
                              StrategyPointerType& pThisStrategy) override
    {
        pThisStrategy = BaseType::mStrategies[rStrategyLabel];
        if(pThisStrategy != nullptr)
            return true;
        return false;
    }

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
        buffer << "FractionalStepSettingsForChimera" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "FractionalStepSettingsForChimera";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{


    ///@}


}; // Class FractionalStepSettingsForChimera

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDenseSpace, class TSparseSpace, class TLinearSolver >
inline std::istream& operator >> (std::istream& rIStream,
                                  FractionalStepSettingsForChimera<TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDenseSpace, class TSparseSpace, class TLinearSolver >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FractionalStepSettingsForChimera<TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FRACTIONAL_STEP_SETTINGS_FOR_CHIMERA_H
