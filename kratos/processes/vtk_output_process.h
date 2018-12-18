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
//

#if !defined(VTK_OUTPUT_PROCESS_H)
#define VTK_OUTPUT_PROCESS_H
// System includes
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip> // for std::setprecision
#include <map>
// For MPI-parallel output
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif
#include "includes/kratos_parameters.h"
#include "includes/io.h"
#include "process.h"
// project includes

namespace Kratos
{
/** \brief VtkOutputProcess
* A simple class that has functionality to write vtk output
*/
class KRATOS_API(KRATOS_CORE) VtkOutputProcess : public Process
{

    typedef ProcessInfo ProcessInfoType;

  public:
    /// Pointer definition of VtkOutputProcess
    KRATOS_CLASS_POINTER_DEFINITION(VtkOutputProcess);

    ///@name Life Cycle
    ///@{

    /**
		Creates a VtkOutputProcess data object
		*/
    VtkOutputProcess(ModelPart &rModelPart, Parameters rParameters);
    /// Destructor.
    virtual ~VtkOutputProcess();

    ///@}
    ///@name Operators
    ///@{

    /// Execute method is used to execute the Process algorithms.
    void Execute() override {}

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override;

    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override;

    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override;

    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override;

    /// this function is designed for being called after ExecuteInitialize ONCE to
    /// verify that the input is correct.
    int Check() override;

    ///@}

    ///@name Access
    ///@{

    ///@

    ///@name Static Operations
    ///@{
    void GetInfo() const;

    ///@}
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << " VtkOutputProcess object " << std::endl;
    }

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    ModelPart &mrModelPart;
    std::string mCaseName;
    std::string mOutputFilename;

    Parameters mrOutputSettings;
    unsigned int mDefaultPrecision;
    std::map<int, int> mKratosIdToVtkId;
    unsigned int mVtkCellListSize;
    unsigned int mStep;
    bool mDoneTest;
    bool mShouldSwap;
    bool mOutputSubModelParts;

    ///@}
    ///@name Operators
    ///@{
    void CreateMapFromKratosIdToVTKId(ModelPart &rModelPart);

    unsigned int DetermineVtkCellListSize(ModelPart &rModelPart);

    void Initialize(ModelPart &rModelPart);

    void WriteHeader(ModelPart &rModelPart);

    void WriteMesh(ModelPart &rModelPart);

    void WriteNodes(ModelPart &rModelPart);

    void WriteConditionsAndElements(ModelPart &rModelPart);

    void WriteConditionAndElementTypes(ModelPart &rModelPart);

    void WriteNodalResultsAsPointData(ModelPart &rModelPart);

    void WriteElementData(ModelPart &rModelPart);

    void WriteHeaderBinary(ModelPart &rModelPart);

    void WriteMeshBinary(ModelPart &rModelPart);

    void WriteNodesBinary(ModelPart &rModelPart);

    void WriteConditionsAndElementsBinary(ModelPart &rModelPart);

    void WriteConditionAndElementTypesBinary(ModelPart &rModelPart);

    void WriteNodalResultsAsPointDataBinary(ModelPart &rModelPart);

    void WriteElementDataBinary(ModelPart &rModelPart);

    void PrintOutputModelPart(ModelPart &modelPart);

    void PrintOutput();

    std::string GetOutputFileName(ModelPart &rModelPart);

    void ForceBigEndian(unsigned char *bytes);
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VtkOutputProcess);
    }

    virtual void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VtkOutputProcess);
    }

    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif // VTK_OUTPUT_PROCESS_H
