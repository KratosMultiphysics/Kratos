#ifndef KRATOS_SOLVER_SETTINGS_H
#define KRATOS_SOLVER_SETTINGS_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "processes/process.h"

// Application includes


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
class SolverSettings
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SolverSettings
    KRATOS_CLASS_POINTER_DEFINITION(SolverSettings);

    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> StrategyType;
    typedef typename StrategyType::Pointer StrategyPointerType;
    typedef typename Process::Pointer ProcessPointerType;
    typedef BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> TBuilderAndSolverType;

    enum StrategyLabel { Velocity, Pressure, /*EddyViscosity,*/ NumLabels };

    enum TurbulenceModelLabel { SpalartAllmaras, NumTurbModels };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SolverSettings(ModelPart& rModelPart,
                   const unsigned int ThisDomainSize,
                   const unsigned int ThisTimeOrder,
                   const bool UseSlip,
                   const bool MoveMeshFlag,
                   const bool ReformDofSet):
        mStrategies(),
        mHaveTurbulenceModel(false),
        mrModelPart(rModelPart),
        mDomainSize(ThisDomainSize),
        mTimeOrder(ThisTimeOrder),
        mEchoLevel(1),
        mUseSlip(UseSlip),
        mReformDofSet(ReformDofSet),
        mMoveMeshFlag(MoveMeshFlag)
    {}

    /// Destructor.
    virtual ~SolverSettings(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    virtual StrategyPointerType pGetStrategy(StrategyLabel const& rStrategyLabel)
    {
        return mStrategies[rStrategyLabel];
    }

    virtual void SetStrategy(StrategyLabel const& rStrategyLabel,
                             StrategyPointerType pStrategy)
    {
        mStrategies[rStrategyLabel] = pStrategy;
    }

    virtual void SetStrategy(StrategyLabel const& rStrategyLabel,
                             typename TLinearSolver::Pointer pLinearSolver,
                             const double Tolerance,
                             const unsigned int MaxIter) = 0;

    virtual void SetTurbulenceModel(TurbulenceModelLabel const& rTurbulenceModel,
                                    typename TLinearSolver::Pointer pLinearSolver,
                                    const double Tolerance,
                                    const unsigned int MaxIter) {}

    virtual void SetTurbulenceModel(ProcessPointerType pTurbulenceModel)
    {
        mpTurbulenceModel = ProcessPointerType(pTurbulenceModel);
        mHaveTurbulenceModel = true;
    }

    virtual bool GetTurbulenceModel(ProcessPointerType& pTurbulenceModel)
    {
        if( mHaveTurbulenceModel )
        {
            pTurbulenceModel.reset();
            pTurbulenceModel = ProcessPointerType(mpTurbulenceModel);
        }

        return mHaveTurbulenceModel;
    }

    virtual unsigned int GetDomainSize() const
    {
        return mDomainSize;
    }

    virtual unsigned int GetTimeOrder() const
    {
        return mTimeOrder;
    }

    virtual bool UseSlipConditions() const
    {
        return mUseSlip;
    }

    virtual bool MoveMesh() const
    {
        return mMoveMeshFlag;
    }

    virtual bool FindStrategy(StrategyLabel const& rStrategyLabel,
                              StrategyPointerType& pThisStrategy)
    {
        typename std::map<StrategyLabel,StrategyPointerType>::iterator itStrategy = mStrategies.find(rStrategyLabel);

        if ( itStrategy != mStrategies.end() )
        {
            if(itStrategy->second != nullptr)
            {
                pThisStrategy.swap(itStrategy->second);
                return true;
            } else {
                KRATOS_INFO("SolverSettingsFSStrategy")<<"Strategy for :: "<<rStrategyLabel<<" not found."<<std::endl;
                return false;
            }
        }
        else {
            KRATOS_INFO("SolverSettingsFSStrategy")<<"Strategy for :: "<<rStrategyLabel<<" not found."<<std::endl;
            return false;
        }
    }

    virtual bool FindTolerance(StrategyLabel const& rStrategyLabel,
                               double& rTolerance)
    {
        typename std::map< StrategyLabel,double >::iterator itTol = mTolerances.find(rStrategyLabel);
        if ( itTol != mTolerances.end() )
        {
            rTolerance = itTol->second;
            return true;
        }
        else
            return false;
    }

    virtual bool FindMaxIter(StrategyLabel const& rStrategyLabel,
                             unsigned int& rMaxIter)
    {
        typename std::map< StrategyLabel,unsigned int >::iterator itMaxIter = mMaxIter.find(rStrategyLabel);
        if ( itMaxIter != mMaxIter.end() )
        {
            rMaxIter = itMaxIter->second;
            return true;
        }
        else
            return false;
    }

    virtual void SetEchoLevel(unsigned int EchoLevel)
    {
        mEchoLevel = EchoLevel;
        for (typename std::map< StrategyLabel, StrategyPointerType>::iterator itStrategy = mStrategies.begin(); itStrategy != mStrategies.end(); ++itStrategy)
            (itStrategy->second)->SetEchoLevel(mEchoLevel);
    }

    virtual unsigned int GetEchoLevel()
    {
        return mEchoLevel;
    }

    bool GetReformDofSet()
    {
        return mReformDofSet;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "SolverSettings" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SolverSettings";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    ModelPart& GetModelPart()
    {
        return mrModelPart;
    }

    std::map< StrategyLabel, StrategyPointerType > mStrategies;

    std::map< StrategyLabel, double > mTolerances;

    std::map< StrategyLabel, unsigned int > mMaxIter;

    ProcessPointerType mpTurbulenceModel;

    bool mHaveTurbulenceModel;

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

    ModelPart& mrModelPart;

    unsigned int mDomainSize;

    unsigned int mTimeOrder;

    unsigned int mEchoLevel;

    bool mUseSlip;

    bool mReformDofSet;

    bool mMoveMeshFlag;


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
    SolverSettings(){}

    /// Assignment operator.
    SolverSettings& operator=(SolverSettings const& rOther){}

    /// Copy constructor.
    SolverSettings(SolverSettings const& rOther){}


    ///@}

}; // Class SolverSettings

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDenseSpace, class TSparseSpace, class TLinearSolver >
inline std::istream& operator >> (std::istream& rIStream,
                                  SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDenseSpace, class TSparseSpace, class TLinearSolver >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SolverSettings<TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SOLVER_SETTINGS_H
