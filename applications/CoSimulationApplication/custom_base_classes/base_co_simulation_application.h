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
    CoSimulationBaseApplication(CoSimulationBaseIo &iIo, Parameters iParameters) : BaseType(*(ModelPart::Pointer(new ModelPart("cosim")))),
                                                                                                         mParameters(iParameters),
                                                                                                         mrIo(iIo)
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

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    // This will update the data field on the mrModelPart :: Done via IO
    virtual void SynchronizeInputData() // TODO: Check how to give scalar and vectorial data field names as arguments
    {
    }

    virtual void SynchronizeOutputData() // TODO: This will take the name of the data field and the other application from where the data is to be obtained.
    {
    }

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
    CoSimulationBaseIo &mrIo;
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
