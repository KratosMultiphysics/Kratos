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
/** \brief VtkOutput
* A simple class that has functionality to write vtk output
*/
class KRATOS_API(KRATOS_CORE) VtkOutput : public IO
{

    typedef ProcessInfo ProcessInfoType;

  public:
    /// Pointer definition of VtkOutput
    KRATOS_CLASS_POINTER_DEFINITION(VtkOutput);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor by passing a vector of Master and slave dofs and corresponding Matrix and constant vector
     * @param rModelPart The modelpart which is used for output
     * @param rParameters Parameters including settings for the output
     */
    VtkOutput(ModelPart &rModelPart, Parameters rParameters);
    /// Destructor.
    virtual ~VtkOutput();

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Prints mrModelPart in VTK format together with the results
     */
    void PrintOutput();

    ///@}

    ///@name Access
    ///@{

    ///@

    ///@name Static Operations
    ///@{
    void GetInfo() const;

    ///@}
    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " VtkOutput object " << std::endl;
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

    /**
     * @brief Print the given rModelPart as VTK file together with the requested results
     * @param rModelPart modelpart which is beging output
     */
    void WriteModelPart(const ModelPart &rModelPart);

    /**
     * @brief Create a map with kratos nodeId as key and VTK nodeId as value. This require for VTK that the node numbers are in sequence.
     * @param rModelPart modelpart which is beging output
     */
    void CreateMapFromKratosIdToVTKId(const ModelPart &rModelPart);

    /**
     * @brief Calculate the total number of cells which are in the provided rModelPart. = num_elements + num_conditions
     *          It is necessary to be known prior to output
     * @param rModelPart modelpart which is beging output
     */
    unsigned int DetermineVtkCellListSize(const ModelPart &rModelPart);

    /**
     * @brief Initialize function for the class
     * @param rModelPart modelpart which is beging output
     */
    void Initialize(const ModelPart &rModelPart);

    /**
     * @brief Write the VTK header for the output of given rModelPart. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteHeader(const ModelPart &rModelPart);

    /**
     * @brief Write the mesh from rModelPart. Nodes, Elements or/and Conditions. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteMesh(const ModelPart &rModelPart);

    /**
     * @brief Write the nodes in the rModelPart. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteNodes(const ModelPart &rModelPart);

    /**
     * @brief Write the elements and conditions in rModelPart. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteConditionsAndElements(const ModelPart &rModelPart);

    /**
     * @brief Write the types for elements and conditions in rModelPart. This is specific for VTK format. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteConditionAndElementTypes(const ModelPart &rModelPart);

    /**
     * @brief Write the results on the nodes. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteNodalResultsAsPointData(const ModelPart &rModelPart);

    /**
     * @brief Write the results/flags on the elements of rModelPart for example : ACTIVE. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteElementData(const ModelPart &rModelPart);

    /**
     * @brief Write the VTK header for the output of given rModelPart. In Binary format
     * @param rModelPart modelpart which is beging output
     */
    void WriteHeaderBinary(const ModelPart &rModelPart);

    /**
     * @brief Write the mesh from rModelPart. Nodes, Elements or/and Conditions. In Binary format
     * @param rModelPart modelpart which is beging output
     */
    void WriteMeshBinary(const ModelPart &rModelPart);

    /**
     * @brief Write the nodes in the rModelPart. In Binary format
     * @param rModelPart modelpart which is beging output
     */
    void WriteNodesBinary(const ModelPart &rModelPart);

    /**
     * @brief Write the elements and conditions in rModelPart. In Binary format
     * @param rModelPart modelpart which is beging output
     */
    void WriteConditionsAndElementsBinary(const ModelPart &rModelPart);

    /**
     * @brief Write the types for elements and conditions in rModelPart. This is specific for VTK format. In Binary format
     * @param rModelPart modelpart which is beging output
     */
    void WriteConditionAndElementTypesBinary(const ModelPart &rModelPart);

    /**
     * @brief Write the results on the nodes. In Binary format
     * @param rModelPart modelpart which is beging output
     */
    void WriteNodalResultsAsPointDataBinary(const ModelPart &rModelPart);

    /**
     * @brief Write the results/flags on the elements of rModelPart for example : ACTIVE. In Binary format
     * @param rModelPart modelpart which is beging output
     */
    void WriteElementDataBinary(const ModelPart &rModelPart);

    /**
     * @brief Get the output file name based on the provided settings and the MPI rank
     * @param rModelPart modelpart which is beging output
     */
    std::string GetOutputFileName(const ModelPart &rModelPart);

    /**
     * @brief Only used in the Binary format output. This function forces the big endian format for the input binary stream
     * @param rModelPart modelpart which is beging output
     */
    void ForceBigEndian(unsigned char *bytes);
    ///@}
    ///@name Serialization
    ///@{

    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif // VTK_OUTPUT_PROCESS_H
