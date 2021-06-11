//
//   Project Name:        KratosPFEMFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:                January 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#ifndef KRATOS_TWO_STEP_VP_SOLVER_SETTINGS_H
#define KRATOS_TWO_STEP_VP_SOLVER_SETTINGS_H

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
    ///@addtogroup PfemFluidDynamicsApplication
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

    /// Helper class to define solution strategies for TwoStepVPStrategy.
    template <class TSparseSpace,
              class TDenseSpace,
              class TLinearSolver>
    class TwoStepVPSolverSettings
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of TwoStepVPSolverSettings
        KRATOS_CLASS_POINTER_DEFINITION(TwoStepVPSolverSettings);

        typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> StrategyType;
        typedef typename StrategyType::Pointer StrategyPointerType;
        typedef typename Process::Pointer ProcessPointerType;

        enum StrategyLabel
        {
            Velocity,
            Pressure,
            /*EddyViscosity,*/ NumLabels
        };

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor.
        TwoStepVPSolverSettings(ModelPart &rModelPart,
                                const unsigned int ThisDomainSize,
                                const unsigned int ThisTimeOrder,
                                const bool ReformDofSet) : mStrategies(),
                                                           mrModelPart(rModelPart),
                                                           mDomainSize(ThisDomainSize),
                                                           mTimeOrder(ThisTimeOrder),
                                                           mEchoLevel(1),
                                                           mReformDofSet(ReformDofSet)
        {
        }

        /// Destructor.
        virtual ~TwoStepVPSolverSettings() {}

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

        ///@}
        ///@name Access
        ///@{

        virtual StrategyPointerType pGetStrategy(StrategyLabel const &rStrategyLabel)
        {
            return mStrategies[rStrategyLabel];
        }

        virtual void SetStrategy(StrategyLabel const &rStrategyLabel,
                                 StrategyPointerType pStrategy)
        {
            mStrategies[rStrategyLabel] = pStrategy;
        }

        virtual void SetStrategy(StrategyLabel const &rStrategyLabel,
                                 typename TLinearSolver::Pointer pLinearSolver,
                                 const double Tolerance,
                                 const unsigned int MaxIter) = 0;

        virtual unsigned int GetDomainSize() const
        {
            return mDomainSize;
        }

        virtual bool FindStrategy(StrategyLabel const &rStrategyLabel,
                                  StrategyPointerType &pThisStrategy)
        {
            typename std::map<StrategyLabel, StrategyPointerType>::iterator itStrategy = mStrategies.find(rStrategyLabel);

            if (itStrategy != mStrategies.end())
            {
                pThisStrategy.swap(itStrategy->second);
                return true;
            }
            else
                return false;
        }

        virtual bool FindTolerance(StrategyLabel const &rStrategyLabel,
                                   double &rTolerance)
        {
            typename std::map<StrategyLabel, double>::iterator itTol = mTolerances.find(rStrategyLabel);
            if (itTol != mTolerances.end())
            {
                rTolerance = itTol->second;
                return true;
            }
            else
                return false;
        }

        virtual bool FindMaxIter(StrategyLabel const &rStrategyLabel,
                                 unsigned int &rMaxIter)
        {
            typename std::map<StrategyLabel, unsigned int>::iterator itMaxIter = mMaxIter.find(rStrategyLabel);
            if (itMaxIter != mMaxIter.end())
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
            for (typename std::map<StrategyLabel, StrategyPointerType>::iterator itStrategy = mStrategies.begin(); itStrategy != mStrategies.end(); ++itStrategy)
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
            buffer << "TwoStepVPSolverSettings";
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "TwoStepVPSolverSettings"; }

        /// Print object's data.
        virtual void PrintData(std::ostream &rOStream) const {}

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

        ModelPart &GetModelPart()
        {
            return mrModelPart;
        }

        std::map<StrategyLabel, StrategyPointerType> mStrategies;

        std::map<StrategyLabel, double> mTolerances;

        std::map<StrategyLabel, unsigned int> mMaxIter;

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

        ModelPart &mrModelPart;

        unsigned int mDomainSize;

        unsigned int mEchoLevel;

        unsigned int mTimeOrder;

        bool mReformDofSet;

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
        TwoStepVPSolverSettings() {}

        /// Assignment operator.
        TwoStepVPSolverSettings &operator=(TwoStepVPSolverSettings const &rOther) {}

        /// Copy constructor.
        TwoStepVPSolverSettings(TwoStepVPSolverSettings const &rOther) {}

        ///@}

    }; // Class TwoStepVPSolverSettings

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    template <class TDenseSpace, class TSparseSpace, class TLinearSolver>
    inline std::istream &operator>>(std::istream &rIStream,
                                    TwoStepVPSolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> &rThis)
    {
        return rIStream;
    }

    /// output stream function
    template <class TDenseSpace, class TSparseSpace, class TLinearSolver>
    inline std::ostream &operator<<(std::ostream &rOStream,
                                    const TwoStepVPSolverSettings<TSparseSpace, TDenseSpace, TLinearSolver> &rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_VP_SOLVER_SETTINGS_H
