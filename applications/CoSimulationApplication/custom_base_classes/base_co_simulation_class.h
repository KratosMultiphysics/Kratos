//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

#if !defined(KRATOS_CO_SIMULATION_BASE_CLASS_H_INCLUDED)
#define KRATOS_CO_SIMULATION_BASE_CLASS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class CoSimulationBaseClass : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    KRATOS_CLASS_POINTER_DEFINITION(CoSimulationBaseClass);
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationBaseClass(ModelPart &modelPart) : BaseType(modelPart)
    {
    }

    ~CoSimulationBaseClass()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    //// Overriding the methods of solving strategy
    virtual void Predict() override
    {
    }

    virtual void Initialize() override
    {
    }

    virtual double Solve() override
    {
        return 0.0;
    }

    virtual void Clear() override
    {
    }

    virtual bool IsConverged() override
    {
        return mIsConverged;
    }

    virtual void CalculateOutputData() override
    {
    }

    virtual void InitializeSolutionStep() override
    {
    }

    virtual void FinalizeSolutionStep() override
    {
    }

    virtual bool SolveSolutionStep() override
    {
        return true;
    }

    virtual void SetEchoLevel(int Level) override
    {
        mEchoLevel = Level;
    }

    virtual int GetEchoLevel() override
    {
        return mEchoLevel;
    }

    virtual void SetRebuildLevel(int Level) override
    {
    }

    virtual int GetRebuildLevel() override
    {
        return 0;
    }

    virtual double GetResidualNorm() override
    {
        return 0.0;
    }

    virtual int Check() override
    {
        std::cout<<"CoSimulationBaseClass :: Check "<<std::endl;
        return 0;
    }

    /// Methods specific for Co-Simulation

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

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
    int mEchoLevel;
    bool mIsConverged;

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

    ///@}
};
}

#endif