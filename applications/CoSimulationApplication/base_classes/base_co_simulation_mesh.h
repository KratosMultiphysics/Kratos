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

#if !defined(KRATOS_CO_SIMULATION_BASE_MESH_H_INCLUDED)
#define KRATOS_CO_SIMULATION_BASE_MESH_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <stdio.h>

// Application includes
#include "custom_utilities/co_simulation_include.h"

class CoSimulationMesh
{

  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(CoSimulationMesh);        

    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationMesh(std::string iName, unsigned int iSize) : mName(iName), mSize(iSize)
    {
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
    virtual std::string Name(){ return mName; }

    /// Data synchronization methods
    /* 
     * This function should read/take the specified output datafields of a
     * solver available. Once it reads/takes the datafield it SHOULD call the function
     * MakeDataFieldNotAvailable to let the other solvers that the solver took its input
     */
    virtual void SetSize(const unsigned int iSize)
    {
    }

    virtual void SetMeshData(std::string iMeshName, DataLocationOnMesh iLocationOnMesh)
    {
        mMeshName = iMeshName;
    }

    virtual void ExportDetails(std::ofstream& iFile)
    {
    }

    virtual size_t GetDataFieldSize(){
        return 0;
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
    const std::string mName;
    unsigned int mSize;
    DataLocationOnMesh mLocation;
    std::string mMeshName;
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


#endif // KRATOS_CO_SIMULATION_BASE_DATAFIELD_H_INCLUDED  defined
