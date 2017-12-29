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

#if !defined(KRATOS_CO_SIMULATION_BASE_COUPLING_STRATEGY_H_INCLUDED)
#define KRATOS_CO_SIMULATION_BASE_COUPLING_STRATEGY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "solving_strategies/strategies/solving_strategy.h"
#include "includes/model_part.h"

// Application includes
#include "custom_base_classes/base_co_simulation_application.h"
#include "custom_base_classes/base_co_simulation_convergence_acceleration_scheme.h"

namespace Kratos
{

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class CoSimulationBaseCouplingStrategy : public CoSimulationBaseClass<TSparseSpace, TDenseSpace, TLinearSolver>
{

  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(CoSimulationBaseCouplingStrategy);
    typedef CoSimulationBaseClass<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::Pointer CoSimBaseClassPointerType;
    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationBaseCouplingStrategy(CoSimBaseClassPointerType app1, CoSimBaseClassPointerType app2, CoSimulationBaseConvergenceAccelerationScheme& iConvAccelerator) : BaseType(*(ModelPart::Pointer(new ModelPart("cosim")))), 
                                                                                                                                                                        mpApplicationOne(app1), 
                                                                                                                                                                        mpApplicationTwo(app2),
                                                                                                                                                                        mrConvAccelerator(iConvAccelerator)
    {
    }

    virtual ~CoSimulationBaseCouplingStrategy()
    {
    }
    ///@}
    ///@name Operators
    ///@{

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
        return true;
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
        return 0;
    }

    /// Methods specific for Co-Simulation
    
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

    CoSimBaseClassPointerType mpApplicationOne;
    CoSimBaseClassPointerType mpApplicationTwo;
    CoSimulationBaseConvergenceAccelerationScheme& mrConvAccelerator;

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

    ///@}

}; // End class

} // End Kratos namespace

#endif // KRATOS_CO_SIMULATION_BASE_COUPLING_STRATEGY_H_INCLUDED  defined
