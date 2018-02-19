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
#include <map>

// Project includes

// Application includes
#include "base_classes/base_co_simulation_application_io.h"
#include "base_classes/base_co_simulation_data.h"
#include "base_classes/base_co_simulation_mesh.h"

class CoSimulationBaseSolver
{

  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(CoSimulationBaseSolver);
    typedef typename CoSimulationData<double>::Pointer DataPointerType;
    typedef typename CoSimulationBaseIo::Pointer BaseIoPointerType;
    typedef typename CoSimulationMesh::Pointer MeshPointerType;
    ///@}
    ///@name Life Cycle
    ///@{
    CoSimulationBaseSolver(std::string iName) : mName(iName)
    {
    }

    virtual ~CoSimulationBaseSolver()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    std::string Name() { return this->mName; }

    virtual void Initialize(std::string iMode)
    {
        mpIo->Initialize(this->mName, iMode);
    }

    /////////////////////////////////////////////////
    /// Methods specific for Co-Simulation
    /////////////////////////////////////////////////

    virtual DataPointerType AddCoSimulationData(std::string iName, size_t iSize)
    {
        DataPointerType data = DataPointerType(new CoSimulationData<double>(iName, iSize));
        mCoSimulationDataMap[iName] = data;
        return data;
    }

    virtual DataPointerType AddCoSimulationData(std::string iName, std::string iMeshName)
    {
        size_t size = mCoSimulationMeshMap[iMeshName]->GetDataFieldSize();
        DataPointerType data = DataPointerType(new CoSimulationData<double>(iName, size));
        mCoSimulationDataMap[iName] = (data);
        return data;
    }

    virtual void SetIo(BaseIoPointerType iIo)
    {
        mpIo = iIo;
    }

    /// Data synchronization methods
    virtual void ExportData(DataPointerType iData, Pointer iToApp)
    {
        mpIo->ExportData(iData, this->Name(), iToApp->Name());
    }

    virtual void ImportData(DataPointerType iData, Pointer iFromApp)
    {

        if (mCoSimulationDataMap.count(iData->Name())) // If its available overrite it
        {
        }
        else // If not create it first then fill it.
        {
        }

        mpIo->ImportData(iData, iFromApp->Name(), this->Name());
    }

    virtual void MakeDataAvailable(DataPointerType iData, Pointer iToApp)
    {
        mpIo->MakeDataAvailable(iData, this->Name(), iToApp->Name());
    }

    virtual void ExportMesh(MeshPointerType iMesh, Pointer iToApp)
    {
        mpIo->ExportMesh(iMesh, this->Name(), iToApp->Name());
    }

    virtual void ImportMesh(MeshPointerType iMesh, Pointer iFromApp)
    {
        mpIo->ImportMesh(iMesh, iFromApp->Name(), this->Name());
    }

    virtual void MakeMeshAvailable(MeshPointerType iMesh, Pointer iToApp)
    {
        mpIo->MakeMeshAvailable(iMesh, this->Name(), iToApp->Name());
    }

    // Control functions
    virtual void SolveSolutionStep()
    {
    }

    virtual void ExecuteBeforeSolve()
    {
    }

    virtual void ExecuteAfterSolve()
    {
    }

    virtual void Finalize()
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
    std::string mName;
    BaseIoPointerType mpIo;
    std::map<std::string, DataPointerType> mCoSimulationDataMap;
    std::map<std::string, MeshPointerType> mCoSimulationMeshMap;
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

#endif // KRATOS_CO_SIMULATION_APPLICATION_H_INCLUDED  defined
