//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala, Philipp Bucher
//
//

#if !defined( KRATOS_VTK_OUTPUT_H_INCLUDED )
#define KRATOS_VTK_OUTPUT_H_INCLUDED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/io.h"


namespace Kratos
{
/** \brief VtkOutput
* A simple class that has functionality to write vtk output
* @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*/
class KRATOS_API(KRATOS_CORE) VtkOutput : public IO
{
public:
    /// Pointer definition of VtkOutput
    KRATOS_CLASS_POINTER_DEFINITION(VtkOutput);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor by passing a ModelPart and Kratos-Parameters
     * @param rModelPart The modelpart which is used for output
     * @param Parameters Parameters including settings for the output
     */
    VtkOutput(ModelPart& rModelPart, Parameters Parameters);

    /// Destructor.
    virtual ~VtkOutput() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Prints mrModelPart in VTK format together with the results
     */
    void PrintOutput();

    ///@}
    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << " VtkOutput object " << std::endl;
    }

    enum class FileFormat {
        VTK_ASCII,
        VTK_BINARY
    };

protected:
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;
    VtkOutput::FileFormat mFileFormat;

    Parameters mOutputSettings;
    unsigned int mDefaultPrecision;
    std::unordered_map<int, int> mKratosIdToVtkId;
    bool mShouldSwap = false;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Print the given rModelPart as VTK file together with the requested results
     * @param rModelPart modelpart which is beging output
     * @param IsSubModelPart whether the modelpart is to be treated as a submodelpart
     * this is only relevant for the file-name
     */
    void WriteModelPartToFile(const ModelPart& rModelPart, const bool IsSubModelPart);

    /**
     * @brief Get the output file name based on the provided settings and the MPI rank
     * @param rModelPart modelpart which is beging output
     */
    std::string GetOutputFileName(const ModelPart& rModelPart, const bool IsSubModelPart);

    /**
     * @brief Initialize function for the class
     * @param rModelPart modelpart which is beging output
     */
    void Initialize(const ModelPart& rModelPart);

    /**
     * @brief Create a map with kratos nodeId as key and VTK nodeId as value. This require for VTK that the node numbers are in sequence.
     * @param rModelPart modelpart which is beging output
     */
    void CreateMapFromKratosIdToVTKId(const ModelPart& rModelPart);

    /**
     * @brief Write the VTK header for the output of given rModelPart.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteHeaderToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const;

    /**
     * @brief Write the mesh from rModelPart. Nodes, Elements or/and Conditions.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteMeshToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const;

    /**
     * @brief Write the nodes in the rModelPart.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteNodesToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const;

    /**
     * @brief Write the elements and conditions in rModelPart.
     *        IMPORTANT : Need to write them together because of the CELLS block in VTK format
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteConditionsAndElementsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const;

    /**
     * @brief Calculate the total number of cells which are in the provided rModelPart. = num_elements + num_conditions
     *          It is necessary to be known prior to output
     * @template TContainerType type of container.
     * @param rContainer the container which is beging output
     */
    template<typename TContainerType>
    unsigned int DetermineVtkCellListSize(const TContainerType& rContainer) const;

    /**
     * @brief Write the element/condition WriteConnectivity provided the container they are in
     * @template TEntity Element/Condition
     * @param rContainer The container containing elements/conditions
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TContainerType>
    void WriteConnectivity(const TContainerType& rContainer, std::ofstream& rFileStream) const;

    /**
     * @brief Write the element/condition cell types provided the container they are in
     * @template TEntity Element/Condition
     * @param rContainer The container containing elements/conditions
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TContainerType>
    void WriteCellType(const TContainerType& rContainer, std::ofstream& rFileStream) const;

    /**
     * @brief Write the results on the nodes.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteNodalResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the results/flags on the elements of rModelPart.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteElementResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the results/flags on the conditions of rModelPart.
     * @param rModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteConditionResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream);

    /**
     * @brief Write the results of rNodes. Synchronization is necessary because both local
     * and ghost-node-values are printed in MPI and can overlap!
     * @param rVariableName name of the result to be written.
     * @param rNodes the nodes which is beging output
     * @param IsHistoricalValue whether the values are historical or not
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteNodalContainerResults(
        const std::string& rVariableName,
        const ModelPart::NodesContainerType& rNodes,
        const bool IsHistoricalValue,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the variable results of rContainer (Elements or Conditions).
     * @template TContainerType The type of container of the entity on which the results are to be written
     * @param rVariableName name of the result to be written.
     * @param rContainer the container which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WriteGeometricalContainerResults(const std::string& rVariableName,
                                          const TContainerType& rContainer,
                                          std::ofstream& rFileStream) const;

    /**
     * @brief Writes scalar results of rNodes. Wraps the necessary synchronization-calls
     * @param rNodes the nodes which is beging output
     * @param rVariable Variable of the result to be written.
     * @param IsHistoricalValue whether the values are historical or not
     * @param rFileStream the file stream to which data is to be written.
     */
    template<class TVarType>
    void WriteNodalScalarValues(
        const ModelPart::NodesContainerType& rNodes,
        const TVarType& rVariable,
        const bool IsHistoricalValue,
        std::ofstream& rFileStream) const;

    /**
     * @brief Writes vector results of rNodes. Wraps the necessary synchronization-calls
     * @param rNodes the nodes which is beging output
     * @param rVariable Variable of the result to be written.
     * @param IsHistoricalValue whether the values are historical or not
     * @param rFileStream the file stream to which data is to be written.
     */
    template<class TVarType>
    void WriteNodalVectorValues(
        const ModelPart::NodesContainerType& rNodes,
        const TVarType& rVariable,
        const bool IsHistoricalValue,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the scalar-historical variable results of rContainer.
     * @template TContainerType The type of container of the entity on which the results are to be written
     * @template TVarType The type of Variable of the entity on which the results are to be written
     * @param rContainer the container which is beging output
     * @param rVariable Variable of the result to be written.
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType, class TVarType>
    void WriteScalarSolutionStepVariable(
        const TContainerType& rContainer,
        const TVarType& rVariable,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the vector-historical variable results of rContainer.
     * @template TContainerType The type of container of the entity on which the results are to be written
     * @template TVarType The type of Variable of the entity on which the results are to be written
     * @param rContainer the container which is beging output
     * @param rVariable Variable of the result to be written.
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType, class TVarType>
    void WriteVectorSolutionStepVariable(
        const TContainerType& rContainer,
        const TVarType& rVariable,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the scalar-nonhistorical variable results of rContainer.
     * @template TContainerType The type of container of the entity on which the results are to be written
     * @template TVarType The type of Variable of the entity on which the results are to be written
     * @param rContainer the container which is beging output
     * @param rVariable Variable of the result to be written.
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType, class TVarType>
    void WriteScalarContainerVariable(
        const TContainerType& rContainer,
        const TVarType& rVariable,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the vector-nonhistorical variable results of rContainer.
     * @template TContainerType The type of container of the entity on which the results are to be written
     * @template TVarType The type of Variable of the entity on which the results are to be written
     * @param rContainer the container which is beging output
     * @param rVariable Variable of the result to be written.
     * @param VtkDataType type of vtk data
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType, class TVarType>
    void WriteVectorContainerVariable(
        const TContainerType& rContainer,
        const TVarType& rVariable,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the scalar value to the file provided, takes care of binary and ascii formats
     * @template TData The type of data to be written to the file stream rFileStream
     * @param rData data to be written
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TData>
    void WriteScalarDataToFile(const TData& rData, std::ofstream& rFileStream) const;

    /**
     * @brief Write the vector values to the file provided, takes care of binary and ascii formats
     * @template TData The type of data to be written to the file stream rFileStream
     * @param rData data to be written
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TData>
    void WriteVectorDataToFile(const TData& rData, std::ofstream& rFileStream) const;

    /**
     * @brief Only used in the binary format output.
     * This function forces the big endian format for the input binary stream
     * @param pBytes bytes on which the big endian format is to be applied
     */
    void ForceBigEndian(unsigned char* pBytes) const;

    ///@}
};

} // namespace Kratos

#endif // KRATOS_VTK_OUTPUT_H_INCLUDED
