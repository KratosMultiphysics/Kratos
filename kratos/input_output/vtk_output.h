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
//  Collaborator:    Vicente Mataix Ferrandiz
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
#include "processes/integration_values_extrapolation_to_nodes_process.h"

namespace Kratos
{
/** \brief VtkOutput
* A simple class that has functionality to write vtk output
* @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*/
class KRATOS_API(KRATOS_CORE) VtkOutput : public IO
{
public:

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Pointer definition of VtkOutput
    KRATOS_CLASS_POINTER_DEFINITION(VtkOutput);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor by passing a ModelPart and Kratos-Parameters
     * @param rModelPart The modelpart which is used for output
     * @param Parameters Parameters including settings for the output
     */
    explicit VtkOutput(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})" )
        );

    /// Destructor.
    virtual ~VtkOutput() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    static Parameters GetDefaultParameters();

    /**
     * @brief Prints mrModelPart in VTK format together with the results
     */
    void PrintOutput();

    ///@}

    /// Turn back information as a string.
    std::string Info() const override
    {
        return " VtkOutput object ";
    }

    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " VtkOutput object " << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    enum class FileFormat {
        VTK_ASCII,
        VTK_BINARY
    };

protected:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;                        /// The main model part to post process
    VtkOutput::FileFormat mFileFormat;             /// The file format considered

    Parameters mOutputSettings;                    /// The configuration parameters
    unsigned int mDefaultPrecision;                /// The default precision
    std::unordered_map<int, int> mKratosIdToVtkId; /// The map storing the relationship between the Kratos ID and VTK ID
    bool mShouldSwap = false;                      /// If it should swap

    // pointer to object of the extrapolation from gauss point to nodes process
    IntegrationValuesExtrapolationToNodesProcess::UniquePointer mpGaussToNodesProcess;


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Interpolates the gauss point results on to the node using IntegrationValuesExtrapolationToNodesProcess
     */
    void PrepareGaussPointResults();



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
     * @tparam TContainerType type of container.
     * @param rContainer the container which is beging output
     */
    template<typename TContainerType>
    unsigned int DetermineVtkCellListSize(const TContainerType& rContainer) const;

    /**
     * @brief Write the element/condition WriteConnectivity provided the container they are in
     * @tparam TEntity Element/Condition
     * @param rContainer The container containing elements/conditions
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TContainerType>
    void WriteConnectivity(const TContainerType& rContainer, std::ofstream& rFileStream) const;

    /**
     * @brief Write the element/condition cell types provided the container they are in
     * @tparam TEntity Element/Condition
     * @param rContainer The container containing elements/conditions
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TContainerType>
    void WriteCellType(const TContainerType& rContainer, std::ofstream& rFileStream) const;

    /**
     * @brief It checks if the variable is compatible with the VTK format
     * @param rVariableName name of the result to be written.
     */
    bool IsCompatibleVariable(const std::string& rVariableName) const;

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
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @param rVariableName name of the result to be written.
     * @param rContainer the container which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WriteGeometricalContainerResults(const std::string& rVariableName,
                                          const TContainerType& rContainer,
                                          std::ofstream& rFileStream) const;

    /**
     * @brief Write the variable GP results of rContainer (Elements or Conditions).
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @param rVariableName name of the result to be written.
     * @param rContainer the container which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WriteGeometricalContainerIntegrationResults(const std::string& rVariableName,
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
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @tparam TVarType The type of Variable of the entity on which the results are to be written
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
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @tparam TVarType The type of Variable of the entity on which the results are to be written
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
     * @brief Write the flag results of rContainer.
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @param rContainer the container which is beging output
     * @param Flag The flag to be considered to be written
     * @param rFlagName The name of the flag that will appear on the post file
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WriteFlagContainerVariable(
        const TContainerType& rContainer,
        const Flags Flag,
        const std::string& rFlagName,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the scalar-nonhistorical variable results of rContainer.
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @tparam TVarType The type of Variable of the entity on which the results are to be written
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
     * @brief Write the scalar GP variable results of rContainer.
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @tparam TVarType The type of Variable of the entity on which the results are to be written
     * @param rContainer the container which is beging output
     * @param rVariable Variable of the result to be written.
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType, class TVarType>
    void WriteIntegrationScalarContainerVariable(
        const TContainerType& rContainer,
        const Variable<TVarType>& rVariable,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the vector-nonhistorical variable results of rContainer.
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @tparam TVarType The type of Variable of the entity on which the results are to be written
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
     * @brief Write the vector-GP variable results of rContainer.
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @tparam TVarType The type of Variable of the entity on which the results are to be written
     * @param rContainer the container which is beging output
     * @param rVariable Variable of the result to be written.
     * @param VtkDataType type of vtk data
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType, class TVarType>
    void WriteIntegrationVectorContainerVariable(
        const TContainerType& rContainer,
        const Variable<TVarType>& rVariable,
        std::ofstream& rFileStream) const;

    /**
     * @brief Write the scalar value to the file provided, takes care of binary and ascii formats
     * @tparam TData The type of data to be written to the file stream rFileStream
     * @param rData data to be written
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TData>
    void WriteScalarDataToFile(const TData& rData, std::ofstream& rFileStream) const
    {
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) {
            rFileStream << rData;
        } else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY) {
            TData data = rData;
            ForceBigEndian(reinterpret_cast<unsigned char *>(&data));
            rFileStream.write(reinterpret_cast<char *>(&data), sizeof(TData));
        }
    }

    /**
     * @brief Write the vector values to the file provided, takes care of binary and ascii formats
     * @tparam TData The type of data to be written to the file stream rFileStream
     * @param rData data to be written
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TData>
    void WriteVectorDataToFile(const TData& rData, std::ofstream& rFileStream) const
    {
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) {
            for (const auto& r_data_comp : rData) {
                rFileStream << r_data_comp << " ";
            }
        } else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY) {
            for (const auto& r_data_comp : rData ) {
                float data_comp_local = (float)r_data_comp; // should not be const or a reference for enforcing big endian
                ForceBigEndian(reinterpret_cast<unsigned char *>(&data_comp_local));
                rFileStream.write(reinterpret_cast<char *>(&data_comp_local), sizeof(float));
            }
        }
    }

    /**
     * @brief Only used in the binary format output.
     * This function forces the big endian format for the input binary stream
     * @param pBytes bytes on which the big endian format is to be applied
     */
    void ForceBigEndian(unsigned char* pBytes) const;

    ///@}

private:
    ///@name Operations
    ///@{

    /**
     * @brief Prints the Properties Id as an integer variable in each element/condition
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @param rContainer the container which is being output
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WritePropertiesIdsToFile(
        const TContainerType& rContainer,
        std::ofstream& rFileStream) const;

    /**
     * @brief Prints the Ids of the container entities as an integer variable in entity (e.g. node, element, condition)
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @param rContainer the container which is being output
     * @param DataName name of the data in the vtk file
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WriteIdsToFile(
        const TContainerType& rContainer,
        const std::string DataName,
        std::ofstream& rFileStream) const;


    /**
     * @brief Print the given rModelPart as VTK file together with the requested results (Only for model parts without nodes)
     * @param rModelPart modelpart which is beging output
     */
    void WriteModelPartWithoutNodesToFile(ModelPart& rModelPart);

    ///@}
};

} // namespace Kratos

#endif // KRATOS_VTK_OUTPUT_H_INCLUDED
