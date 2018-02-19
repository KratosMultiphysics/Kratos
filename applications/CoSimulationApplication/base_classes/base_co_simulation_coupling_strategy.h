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
#include "base_classes/base_co_simulation_solver.h"
#include "base_classes/base_co_simulation_convergence_acceleration_scheme.h"

namespace Kratos
{

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class CoSimulationBaseCouplingStrategy : public CoSimulationBaseSolver
{

  public:
    ///@name Type Definitions
    ///@{
    typedef CoSimulationBaseSolver BaseType;
    typedef typename CoSimulationBaseSolver::Pointer CoSimBaseClassPointerType;
    typedef typename CoSimulationBaseConvergenceAccelerationScheme::Pointer CoSimulationBaseConvergenceAccelerationSchemePointerType;
    KRATOS_CLASS_POINTER_DEFINITION(CoSimulationBaseCouplingStrategy);
    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationBaseCouplingStrategy(std::string iName, CoSimBaseClassPointerType app1, 
                                     CoSimBaseClassPointerType app2, 
                                     CoSimulationBaseConvergenceAccelerationSchemePointerType iConvAccelerator=nullptr) 
                                     : BaseType(iName), mpApplicationOne(app1), 
                                       mpApplicationTwo(app2),
                                       mpConvAccelerator(iConvAccelerator)
    {
    }

    virtual ~CoSimulationBaseCouplingStrategy()
    {
    }
    ///@}
    ///@name Operators
    ///@{


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
    CoSimulationBaseConvergenceAccelerationSchemePointerType mpConvAccelerator;

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
