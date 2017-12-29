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

#if !defined(KRATOS_COSIMULATION_BASE_APPLICATION_H_INCLUDED)
#define KRATOS_COSIMULATION_BASE_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "custom_base_classes/base_co_simulation_application_io.h"

namespace Kratos
{

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class CoSimulationBaseApplication : public CoSimulationBaseClass<TSparseSpace, TDenseSpace, TLinearSolver>
{

  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(CoSimulationBaseApplication);
    typedef CoSimulationBaseClass<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationBaseApplication(Parameters iParameters) : BaseType(*(ModelPart::Pointer(new ModelPart("cosim")))),
                                                                                                         mParameters(iParameters)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    virtual ModelPart &GetModelPart()
    {
        return mrModelPart;
    }

    /////////////////////////////////////////////////
    /// Methods specific for Co-Simulation
    /////////////////////////////////////////////////

    /// Data synchronization methods
    virtual void SynchronizeInputData()
    {
    }

    virtual void SynchronizeOutputData()
    {
    }

    /// model part synchronization methods
    virtual void ImportModelPart()
    {
    }    

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

    ///@}
    ///@name Member Variables
    ///@{
    ModelPart mrModelPart;
    Parameters mParameters;
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

#endif // KRATOS_CO_SIMULATION_APPLICATION_H_INCLUDED  defined
