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

    enum FileFormat{
        VTK_ASCII,
        VTK_BINARY
    };

    enum WriteDataType{
        VTK_SCALAR,
        VTK_VECTOR
    };

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    ModelPart &mrModelPart;
    std::string mCaseName;
    std::string mOutputFilename;
    VtkOutput::FileFormat mFileFormat;

    Parameters mrOutputSettings;
    unsigned int mDefaultPrecision;
    std::map<int, int> mKratosIdToVtkId;
    unsigned int mVtkCellListSize;
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
     * @brief Initialize function for the class
     * @param rModelPart modelpart which is beging output
     */
    void Initialize(const ModelPart &rModelPart);

    /**
     * @brief Write the VTK header for the output of given rModelPart. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteHeader(const ModelPart &rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the mesh from rModelPart. Nodes, Elements or/and Conditions. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteMesh(const ModelPart &rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the nodes in the rModelPart. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteNodes(const ModelPart &rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the elements and conditions in rModelPart. In ASCII format
     *        IMPORTANT : Need to write them together because of the CELLS block in VTK format
     * @param rModelPart modelpart which is beging output
     */
    void WriteConditionsAndElements(const ModelPart &rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the types for elements and conditions in rModelPart. This is specific for VTK format. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteConditionAndElementTypes(const ModelPart &rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the results on the nodes. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteNodalResultsAsPointData(const ModelPart &rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the results/flags on the elements of rModelPart for example : ACTIVE. In ASCII format
     * @param rModelPart modelpart which is beging output
     */
    void WriteElementData(const ModelPart &rModelPart, std::ofstream& rFileStream);

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
     * @brief Write the scalar value to the file provided, takes care of binary and ascii formats
     * @template TData The type of data to be written to the file stream rFileStream
     * @param rData data to be written
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TData>
    void WriteScalarDataToFile(const TData& rData, std::ofstream& rFileStream);


    /**
     * @brief Write the vector values to the file provided, takes care of binary and ascii formats
     * @template TData The type of data to be written to the file stream rFileStream
     * @param rData data to be written
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TData>
    void WriteVectorDataToFile(const TData& rData, std::ofstream& rFileStream);

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
