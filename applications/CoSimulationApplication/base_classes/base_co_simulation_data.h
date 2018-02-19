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

#if !defined(KRATOS_CO_SIMULATION_BASE_DATAFIELD_H_INCLUDED)
#define KRATOS_CO_SIMULATION_BASE_DATAFIELD_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>

// Application includes
#include "custom_utilities/co_simulation_include.h"

template <typename Type>
class CoSimulationData
{

  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(CoSimulationData);

    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationData(std::string iName, unsigned int iSize) : mData(nullptr), mName(iName), mSize(iSize)
    {
        mLocationOnMesh = ON_NONE;
        mMeshName = "NONE";
    }

    virtual ~CoSimulationData()
    {
        if(mData != nullptr)
        {
            delete [] mData;
        }
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
    virtual std::string Name() { return mName; }

    /* 
     * This function should read/take the specified output datafields of a
     * solver available. Once it reads/takes the datafield it SHOULD call the function
     * MakeDataFieldNotAvailable to let the other solvers that the solver took its input
     */
    virtual void SetSize(const unsigned int iSize)
    {
    }

    virtual void Clear()
    {
        if(mData != nullptr)
        {
            delete [] mData;
        }
    }

    virtual void SetMeshData(std::string iMeshName, DataLocationOnMesh iLocationOnMesh)
    {
        mMeshName = iMeshName;
        mLocationOnMesh = iLocationOnMesh;
    }

    virtual void ExportDetails(std::ofstream &iFile)
    {
        iFile << "Name : " << mName << std::endl;
        iFile << "Size : " << mSize << std::endl;
        if (mLocationOnMesh != ON_NONE)
            iFile << "LocationOnMesh : " << mLocationOnMesh << std::endl;
        if (mMeshName != "NONE")
            iFile << "Mesh : " << mMeshName << std::endl;
    }

    virtual void ImportDetails(std::ifstream &iFile)
    {
        if (iFile)
        {
            std::string line;
            while (getline(iFile, line))
            {
                std::istringstream is_line(line);
                std::string key;
                if (std::getline(is_line, key, ':'))
                {
                    std::string value;
                    if (std::getline(is_line, value))
                    {
                        if (key == "Name ")
                        {
                            mName = value;
                        }
                        if (key == "Size ")
                        {
                          is_line >> mSize;
                        }
                        if (key == "LocationOnMesh ")
                        {
                            is_line >> value;
                        }
                        if (key == "Mesh ")
                        {
                            is_line >> value;
                        }
                    }
                }
            }
        }
        else
        {
            assert(false);
        }
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
    Type *mData;
    std::string mName;
    unsigned int mSize;
    DataLocationOnMesh mLocationOnMesh;
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
