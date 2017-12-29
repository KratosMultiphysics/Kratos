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

#if !defined(KRATOS_CO_SIMULATION_BASE_APPLICATION_IO_H_INCLUDED)
#define KRATOS_CO_SIMULATION_BASE_APPLICATION_IO_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/kratos_parameters.h"

// Application includes

namespace Kratos
{

class CoSimulationBaseIo
{

  public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationBaseIo(Parameters iParameters) : mParameters(iParameters)
    {
        mpModelPart = nullptr;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
    ModelPart::Pointer mpModelPart;
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

#endif // KRATOS_CO_SIMULATION_BASE_APPLICATION_IO_H_INCLUDED  defined
