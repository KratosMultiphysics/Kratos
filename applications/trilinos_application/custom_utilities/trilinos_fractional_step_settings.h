#ifndef KRATOS_TRILINOS_FRACTIONAL_STEP_SETTINGS_H
#define KRATOS_TRILINOS_FRACTIONAL_STEP_SETTINGS_H

// System includes
#include <string>
#include <iostream>


// External includes
#include "Epetra_MpiComm.h"

// Project includes
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "processes/process.h"

// Application includes
#include "custom_processes/trilinos_spalart_allmaras_turbulence_model.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"

// FluidDynamicsApplication dependences
#include "../FluidDynamicsApplication/custom_utilities/solver_settings.h"

namespace Kratos
{
///@addtogroup TrilinosApplication
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
class TrilinosFractionalStepSettings: public SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosFractionalStepSettings
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosFractionalStepSettings);

    typedef SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::StrategyType StrategyType;
    typedef typename BaseType::StrategyPointerType StrategyPointerType;
    typedef typename BaseType::ProcessPointerType ProcessPointerType;

    typedef typename BaseType::StrategyLabel StrategyLabel;
    typedef typename BaseType::TurbulenceModelLabel TurbulenceModelLabel;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    TrilinosFractionalStepSettings(Epetra_MpiComm& rComm,
                                   ModelPart& rModelPart,
                                   const unsigned int ThisDomainSize,
                                   const unsigned int ThisTimeOrder,
                                   const bool UseSlip,
                                   const bool MoveMeshFlag,
                                   const bool ReformDofSet):
        SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver>(rModelPart,ThisDomainSize,ThisTimeOrder,UseSlip,MoveMeshFlag,ReformDofSet),
        mrComm(rComm)
    {}

    /// Destructor.
    virtual ~TrilinosFractionalStepSettings(){}

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

        // pointer types for solution strategy construcion
        typedef typename Scheme< TSparseSpace, TDenseSpace >::Pointer SchemePointerType;
        //~ typedef typename ConvergenceCriteria< TSparseSpace, TDenseSpace >::Pointer ConvergenceCriteriaPointerType;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;

        // Default, fixed flags
        bool CalculateReactions = false;
        bool CalculateNormDxFlag = true;

        // Trilinos defaults
        int RowSizeGuess;
        if(this->GetDomainSize() == 2)
            RowSizeGuess = 15;
        else
            RowSizeGuess = 40;

        ModelPart& rModelPart = this->GetModelPart();
        bool ReformDofSet = this->GetReformDofSet();
        bool UseSlip = this->UseSlipConditions();
        unsigned int EchoLevel = this->GetEchoLevel();

        if ( rStrategyLabel == BaseType::Velocity )
        {
            // Velocity Builder and Solver
            BuilderSolverTypePointer pBuildAndSolver = BuilderSolverTypePointer(new TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > (mrComm,RowSizeGuess,pLinearSolver));

            SchemePointerType pScheme;
            //initializing fractional velocity solution step
            if (UseSlip)
            {
                double DomainSize = this->GetDomainSize();
                SchemePointerType Temp = SchemePointerType(new ResidualBasedIncrementalUpdateStaticSchemeSlip< TSparseSpace, TDenseSpace > (DomainSize,DomainSize));
                pScheme.swap(Temp);
            }
            else
            {
                SchemePointerType Temp = SchemePointerType(new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());
                pScheme.swap(Temp);
            }

            // Strategy
            this->mStrategies[BaseType::Velocity] = StrategyPointerType(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pLinearSolver, pBuildAndSolver, CalculateReactions, ReformDofSet, CalculateNormDxFlag));

        }
        else if ( rStrategyLabel == BaseType::Pressure )
        {
            // Pressure Builder and Solver
            BuilderSolverTypePointer pBuildAndSolver = BuilderSolverTypePointer(new TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver >(mrComm,RowSizeGuess,pLinearSolver));
            SchemePointerType pScheme = SchemePointerType(new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

            // Strategy
            this->mStrategies[BaseType::Pressure] = StrategyPointerType(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (rModelPart, pScheme, pLinearSolver, pBuildAndSolver, CalculateReactions, ReformDofSet, CalculateNormDxFlag));
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Error in TrilinosFractionalStepSettings: Unknown strategy label.","");
        }

        this->mTolerances[rStrategyLabel] = Tolerance;

        this->mMaxIter[rStrategyLabel] = MaxIter;

        this->mStrategies[rStrategyLabel]->SetEchoLevel(EchoLevel);

        KRATOS_CATCH("");
    }

    void SetTurbulenceModel(TurbulenceModelLabel const& rTurbulenceModel,
                            typename TLinearSolver::Pointer pLinearSolver,
                            const double Tolerance,
                            const unsigned int MaxIter) override
    {
        KRATOS_TRY;

        this->mHaveTurbulenceModel = true;

        ModelPart& rModelPart = this->GetModelPart();
        double DomainSize = this->GetDomainSize();
        bool ReformDofSet = this->GetReformDofSet();
        unsigned int TimeOrder = this->GetTimeOrder();

        if (rTurbulenceModel == BaseType::SpalartAllmaras)
        {
            this->mpTurbulenceModel = ProcessPointerType( new TrilinosSpalartAllmarasTurbulenceModel<TSparseSpace,TDenseSpace,TLinearSolver>(mrComm,rModelPart,pLinearSolver,DomainSize,Tolerance,MaxIter,ReformDofSet,TimeOrder));
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Error in TrilinosFractionalStepSettings: Unknown turbulence model label.","");
        }

        KRATOS_CATCH("");
    }

    void SetTurbulenceModel(ProcessPointerType pTurbulenceModel) override
    {
        BaseType::SetTurbulenceModel(pTurbulenceModel);
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
        buffer << "TrilinosFractionalStepSettings" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "TrilinosFractionalStepSettings";}

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

    Epetra_MpiComm& mrComm;

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
    TrilinosFractionalStepSettings(){}

    /// Assignment operator.
    TrilinosFractionalStepSettings& operator=(TrilinosFractionalStepSettings const& rOther){}

    /// Copy constructor.
    TrilinosFractionalStepSettings(TrilinosFractionalStepSettings const& rOther){}


    ///@}

}; // Class TrilinosFractionalStepSettings

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDenseSpace, class TSparseSpace, class TLinearSolver >
inline std::istream& operator >> (std::istream& rIStream,
                                  TrilinosFractionalStepSettings<TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDenseSpace, class TSparseSpace, class TLinearSolver >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TrilinosFractionalStepSettings<TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_FRACTIONAL_STEP_SETTINGS_H
