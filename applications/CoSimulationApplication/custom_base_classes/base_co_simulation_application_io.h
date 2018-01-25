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

#if !defined(KRATOS_CO_SIMULATION_BASE_IO_H_INCLUDED)
#define KRATOS_CO_SIMULATION_BASE_IO_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <stdio.h>

// Application includes
#include "custom_utilities/co_simulation_include.h"
#include "custom_base_classes/base_co_simulation_data.h"
#include "custom_base_classes/base_co_simulation_mesh.h"

class CoSimulationBaseIo
{

  public:
    ///@name Type Definitions
    ///@{
    typedef std::shared_ptr<CoSimulationBaseIo> Pointer;
    typedef CoSimulationData<double>::Pointer DataPointerType;
    typedef CoSimulationMesh::Pointer MeshPointerType;    

    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationBaseIo()
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

    /// Data synchronization methods
    /* 
     * This function should read/take the specified output datafields of a
     * solver available. Once it reads/takes the datafield it SHOULD call the function
     * MakeDataFieldNotAvailable to let the other solvers that the solver took its input
     */
    virtual void ExportData(DataPointerType iData, std::string iFrom, std::string iTo)
    {
    }

    /* 
     * This function should write/put out the specified output datafields of a 
     * solver available. Once this function write/put out the datafield it SHOULD call
     * function MakeDataFieldAvailable so as to notifiy other solvers the field is available
     */
    virtual void ImportData(DataPointerType iData, std::string iFrom, std::string iTo)
    {
    }

    virtual void ExportMesh(MeshPointerType iData, std::string iFrom, std::string iTo)
    {
    }

    virtual void ImportMesh(MeshPointerType iData, std::string iFrom, std::string iTo)
    {
    }

    virtual void MakeDataAvailable(DataPointerType iData, std::string iFrom, std::string iTo)
    {
        std::string AvailFileName = (dot + slash + dot + iTo + slash + "DATA" + dot + iData->Name() + dot + availExtension);
        std::ofstream outputFile(AvailFileName.c_str());
        if(outputFile.is_open())
        {
            iData -> PrintDetails(outputFile);
        }
    }

    virtual void MakeMeshAvailable(MeshPointerType iMesh, std::string iFrom, std::string iTo)
    {
        std::string AvailFileName = (dot + slash + dot + iTo + slash + "MESH" + dot + iMesh->Name() + dot + availExtension);
        std::ofstream outputFile(AvailFileName.c_str());
        if(outputFile.is_open())
        {
            iMesh -> PrintDetails(outputFile);
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
    bool MakeDataFieldNotAvailable(std::string iFileName)
    {
        return (remove(iFileName.c_str()) != 0);
    }

    virtual void ReadDataFieldDetails(DataPointerType iData, std::string iAvailFileName)
    {
    }

    /// Check if an input datafield is available.
    /// This is done by checking if a file with the name of the data field exists or not.
    /// This is also for synchronizing between different solvers.
    virtual bool IsDataAvailable(DataPointerType iData, std::string iFor)
    {
        std::string AvailFileName = (dot + slash + dot + iFor + slash + "DATA" + dot + iData->Name() + dot + availExtension);
        while (CoSimulation_FileExists(AvailFileName.c_str()))
        { // file not available
            std::cout << "Data :: " << iData->Name() << " Not available .. waiting ... " << std::endl;
            CoSimulation_Wait(1);
        }
        CoSimulation_Wait(2);

        // Once we know the data is available we read the details of the data
        ReadDataFieldDetails(iData, AvailFileName);

        return true;
    }

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

#endif // KRATOS_CO_SIMULATION_BASE_IO_H_INCLUDED  defined
