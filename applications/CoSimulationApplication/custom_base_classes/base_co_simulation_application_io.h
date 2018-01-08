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
#include <stdio.h>

// Project includes
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_utilities/co_simulation_io_utilities.h"

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
    /* 
     * This function should read/take the specified output datafields of a
     * solver available. Once it reads/takes the datafield it SHOULD call the function
     * MakeDataFieldNotAvailable to let the other solvers that the solver took its input
     */
    virtual void SynchronizeInputData()
    {
    }

    /* 
     * This function should write/put out the specified output datafields of a 
     * solver available. Once this function write/put out the datafield it SHOULD call
     * function MakeDataFieldAvailable so as to notifiy other solvers the field is available
     */
    virtual void SynchronizeOutputData()
    {
    }

    /// model part synchronization methods
    virtual void ImportModelPart()
    {
    }

    /// Check if an input datafield is available.
    /// This is done by checking if a file with the name of the data field exists or not.
    /// This is also for synchronizing between different solvers.
    virtual bool IsDataFieldAvailable(std::string iDataFieldName)
    {
        return CoSimulation_IsDataFieldAvailable(iDataFieldName.c_str());
    }

    bool MakeDataFieldAvailable(std::string iFileName)
    {
        std::ofstream outputFile(iFileName.c_str());
        return outputFile.is_open();
    }

    bool MakeDataFieldNotAvailable(std::string iFileName)
    {
        return (remove(iFileName.c_str()) != 0);
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
    bool FileExists(std::string iFileName)
    {
        return (access(iFileName.c_str(), F_OK) != -1);
    }

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
