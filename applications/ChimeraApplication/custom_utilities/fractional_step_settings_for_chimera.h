#ifndef KRATOS_FRACTIONAL_STEP_SETTINGS_FOR_CHIMERA_H
#define KRATOS_FRACTIONAL_STEP_SETTINGS_FOR_CHIMERA_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
//#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "processes/process.h"

// Application includes
//#include "custom_processes/spalart_allmaras_turbulence_model_for_chimera.h"
#include "custom_utilities/solver_settings_for_chimera.h"
//#include "custom_strategies/custom_builder_and_solver/residualbased_block_builder_and_solver_with_mpc_chimera.h"
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
class FractionalStepSettingsForChimera: public SolverSettingsForChimera<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FractionalStepSettingsForChimera
    KRATOS_CLASS_POINTER_DEFINITION(FractionalStepSettingsForChimera);

    typedef SolverSettingsForChimera<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
    typedef typename BaseType::StrategyType StrategyType;
    typedef typename BaseType::StrategyPointerType StrategyPointerType;
    typedef typename BaseType::ProcessPointerType ProcessPointerType;

    typedef typename BaseType::StrategyLabel StrategyLabel;
    //typedef typename BaseType::TurbulenceModelLabel TurbulenceModelLabel;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FractionalStepSettingsForChimera(ModelPart& rModelPart,
                   const std::size_t ThisDomainSize,
                   const std::size_t ThisTimeOrder,
                   const bool UseSlip,
                   const bool MoveMeshFlag,
                   const bool ReformDofSet):
        BaseType(rModelPart,ThisDomainSize,ThisTimeOrder,UseSlip,MoveMeshFlag,ReformDofSet)
    {}

    /// Destructor.
    ~FractionalStepSettingsForChimera() override{}

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
                             const std::size_t MaxIter) override
    {
        KRATOS_TRY;

        // pointer types for solution strategy construcion
        typedef typename Scheme< TSparseSpace, TDenseSpace >::Pointer SchemePointerType;
//         typedef typename ConvergenceCriteria< TSparseSpace, TDenseSpace >::Pointer ConvergenceCriteriaPointerType;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;

        // Default, fixed flags
        bool CalculateReactions = false;
        bool CalculateNormDxFlag = true;

        ModelPart& rModelPart = BaseType::GetModelPart();
        bool UseSlip = BaseType::UseSlipConditions();
        // Modification of the DofSet is managed by the fractional step strategy, not the auxiliary velocity and pressure strategies.
        bool ReformDofSet = false; //BaseType::GetReformDofSet();
        std::size_t EchoLevel = BaseType::GetEchoLevel();
        std::size_t StrategyEchoLevel = (EchoLevel > 0) ? (EchoLevel-1) : 0;

        if ( rStrategyLabel == BaseType::Velocity )
        {
            // Velocity Builder and Solver
            BuilderSolverTypePointer pBuildAndSolver = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera<TSparseSpace, TDenseSpace, TLinearSolver >
                                                                                (pLinearSolver));

            SchemePointerType pScheme;
            //initializing fractional velocity solution step
            if (UseSlip)
            {
                std::size_t DomainSize = BaseType::GetDomainSize();
                SchemePointerType Temp = SchemePointerType(new ResidualBasedIncrementalUpdateStaticSchemeSlip< TSparseSpace, TDenseSpace > (DomainSize,DomainSize));
                pScheme.swap(Temp);
            }
            else
            {
                SchemePointerType Temp = SchemePointerType(new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
                pScheme.swap(Temp);
            }

            // Strategy
            BaseType::mStrategies[rStrategyLabel] = StrategyPointerType(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver >
                                                                        (rModelPart, pScheme, pLinearSolver, pBuildAndSolver, CalculateReactions, ReformDofSet, CalculateNormDxFlag));

        }
        else if ( rStrategyLabel == BaseType::Pressure )
        {
            // Pressure Builder and Solver
            BuilderSolverTypePointer pBuildAndSolver = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolverWithConstraintsForChimera<TSparseSpace, TDenseSpace, TLinearSolver > (pLinearSolver));
            SchemePointerType pScheme = SchemePointerType(new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

            // Strategy
            BaseType::mStrategies[rStrategyLabel] = StrategyPointerType(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver >
                                                                        (rModelPart, pScheme, pLinearSolver, pBuildAndSolver, CalculateReactions, ReformDofSet, CalculateNormDxFlag));
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Error in FractionalStepSettingsForChimera: Unknown strategy label.","");
        }

        BaseType::mTolerances[rStrategyLabel] = Tolerance;

        BaseType::mMaxIter[rStrategyLabel] = MaxIter;

        BaseType::mStrategies[rStrategyLabel]->SetEchoLevel(StrategyEchoLevel);

        KRATOS_CATCH("");
    }

    /* void SetTurbulenceModel(TurbulenceModelLabel const& rTurbulenceModel,
                                    typename TLinearSolver::Pointer pLinearSolver,
                                    const double Tolerance,
                                    const std::size_t MaxIter) override
    {
        KRATOS_TRY;

        BaseType::mHaveTurbulenceModel = true;

        ModelPart& rModelPart = BaseType::GetModelPart();
        std::size_t DomainSize = BaseType::GetDomainSize();
        std::size_t TimeOrder = BaseType::GetTimeOrder();
        bool ReformDofSet = BaseType::GetReformDofSet();

        if (rTurbulenceModel == BaseType::SpalartAllmaras)
        {
            BaseType::mpTurbulenceModel = ProcessPointerType( new SpalartAllmarasTurbulenceModelForChimera<TSparseSpace,TDenseSpace,TLinearSolver>
                                                              (rModelPart,pLinearSolver,DomainSize,Tolerance,MaxIter,ReformDofSet,TimeOrder));
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Error in FractionalStepSettingsForChimera: Unknown turbulence model label.","");
        }

        KRATOS_CATCH("");
    }

    void SetTurbulenceModel(ProcessPointerType pTurbulenceModel) override
    {
        BaseType::SetTurbulenceModel(pTurbulenceModel);
    } */

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

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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
    FractionalStepSettingsForChimera(){}

    /// Assignment operator.
    FractionalStepSettingsForChimera& operator=(FractionalStepSettingsForChimera const& rOther){}

    /// Copy constructor.
    FractionalStepSettingsForChimera(FractionalStepSettingsForChimera const& rOther){}


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
