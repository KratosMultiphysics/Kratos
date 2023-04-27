//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/io.h"

namespace Kratos
{
/** \brief ParticleVtkOutput
* A simple class that has functionality to write vtk output
* @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*/
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) ParticleVtkOutput : public IO
{
public:

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Pointer definition of ParticleVtkOutput
    KRATOS_CLASS_POINTER_DEFINITION(ParticleVtkOutput);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor by passing a ModelPart and Kratos-Parameters
     * @param rMPMModelPart The modelpart which is used for output
     * @param Parameters Parameters including settings for the output
     */
    explicit ParticleVtkOutput(
        ModelPart& rMPMModelPart,
        Parameters ThisParameters = Parameters(R"({})" )
        );

    /// Destructor.
    virtual ~ParticleVtkOutput() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters
     */
    static Parameters GetDefaultParameters();

    /**
     * @brief Prints mrMPMModelPart in VTK format together with the results
     */
    void PrintOutput(const std::string& rOutputFilename = "");

    ///@}

    /// Turn back information as a string.
    std::string Info() const override
    {
        return " ParticleVtkOutput object ";
    }

    /**
     * @brief Prints information about the class
     * @param rOStream ostream object where output is printed
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " ParticleVtkOutput object " << std::endl;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    enum class FileFormat {
        VTK_ASCII,
        VTK_BINARY
    };

    enum class OutputEntities {
        ELEMENTS,
        CONDITIONS
    };

protected:
    ///@name Member Variables
    ///@{

    ModelPart& mrMPMModelPart;                          /// The main model part to post process
    ParticleVtkOutput::FileFormat mFileFormat;          /// The file format (ascii or binary) considered
    ParticleVtkOutput::OutputEntities mOutputEntities;  /// The entity (element or condition) to print

    Parameters mOutputSettings;                         /// The configuration parameters
    unsigned int mDefaultPrecision;                     /// The default precision
    bool mShouldSwap = false;                           /// If it should swap

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Print the given rMPMModelPart as VTK file together with the requested results
     * @param rMPMModelPart modelpart which is beging output
     * @param IsSubModelPart whether the modelpart is to be treated as a submodelpart
     * this is only relevant for the file-name
     */
    void WriteMPMModelPartToFile(
        const ModelPart& rMPMModelPart,
        const bool IsSubModelPart,
        const std::string& rOutputFilename
        );

    /**
     * @brief Get the output file name based on the provided settings and the MPI rank
     * @param rMPMModelPart modelpart which is beging output
     */
    std::string GetOutputFileName(
        const ModelPart& rMPMModelPart,
        const bool IsSubModelPart,
        const std::string& rOutputFilename
        );

    /**
     * @brief Write the VTK header for the output of given rMPMModelPart.
     * @param rMPMModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteHeaderToFile(
        const ModelPart& rMPMModelPart,
        std::ofstream& rFileStream
        ) const;

    /**
     * @brief Write the mesh from rMPMModelPart: material point Elements or Conditions
     * @param rMPMModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteMPMMeshToFile(
        const ModelPart& rMPMModelPart,
        std::ofstream& rFileStream
        ) const;

    /**
     * @brief Write the material point elements from rMPMModelPart.
     * @param rMPMModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteMPMElementsToFile(
        const ModelPart& rMPMModelPart,
        std::ofstream& rFileStream
        ) const;

    /**
     * @brief Write the material point conditions from rMPMModelPart.
     * @param rMPMModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteMPMConditionsToFile(
        const ModelPart& rMPMModelPart,
        std::ofstream& rFileStream
        ) const;

    /**
     * @brief Write the results/flags on the elements of rMPMModelPart.
     * @param rMPMModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteMPMResultsToFile(
        const ModelPart& rMPMModelPart,
        std::ofstream& rFileStream
        );

    /**
     * @brief Write the results/flags on the elements of rMPMModelPart.
     * @param rMPMModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteMPMElementResultsToFile(
        const ModelPart& rMPMModelPart,
        std::ofstream& rFileStream
        );

    /**
     * @brief Write the results/flags on the conditions of rMPMModelPart.
     * @param rMPMModelPart modelpart which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    void WriteMPMConditionResultsToFile(
        const ModelPart& rMPMModelPart,
        std::ofstream& rFileStream
        );

    /**
     * @brief It checks if the variable is compatible with the VTK format
     * @param rVariableName name of the result to be written.
     */
    bool IsCompatibleVariable(const std::string& rVariableName) const;

    /**
     * @brief Write the variable results of rContainer (Elements or Conditions).
     * @tparam TContainerType The type of container of the entity on which the results are to be written
     * @param rVariableName name of the result to be written.
     * @param rContainer the container which is beging output
     * @param rFileStream the file stream to which data is to be written.
     */
    template<typename TContainerType>
    void WriteGeometricalContainerResults(
        const std::string& rVariableName,
        const TContainerType& rContainer,
        std::ofstream& rFileStream
        ) const;

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
        std::ofstream& rFileStream
        ) const;

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
        const Variable<TVarType>& rVariable,
        std::ofstream& rFileStream
        ) const;

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
        const Variable<TVarType>& rVariable,
        std::ofstream& rFileStream
        ) const;

    /**
     * @brief Write the vector values to the file provided, takes care of binary and ascii formats
     * @tparam TData The type of data to be written to the file stream rFileStream
     * @param rData data to be written
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TData>
    void WriteVectorDataToFile(
        const TData& rData,
        std::ofstream& rFileStream
        ) const
    {
        if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
            for (const auto& r_data_comp : rData) {
                rFileStream << float(r_data_comp) << " ";
            }
        } else if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_BINARY) {
            for (const auto& r_data_comp : rData ) {
                float data_comp_local = (float)r_data_comp; // should not be const or a reference for enforcing big endian
                ForceBigEndian(reinterpret_cast<unsigned char *>(&data_comp_local));
                rFileStream.write(reinterpret_cast<char *>(&data_comp_local), sizeof(float));
            }
        }
    }

    /**
     * @brief Write the scalar value to the file provided, takes care of binary and ascii formats
     * @tparam TData The type of data to be written to the file stream rFileStream
     * @param rData data to be written
     * @param rFileStream the file stream to which data is to be written.
     */
    template <typename TData>
    void WriteScalarDataToFile(const TData& rData, std::ofstream& rFileStream) const
    {
        if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
            rFileStream << rData;
        } else if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_BINARY) {
            TData data = rData;
            ForceBigEndian(reinterpret_cast<unsigned char *>(&data));
            rFileStream.write(reinterpret_cast<char *>(&data), sizeof(TData));
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
        std::ofstream& rFileStream
        ) const;

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
        const std::string& DataName,
        std::ofstream& rFileStream
        ) const;

    ///@}
};

} // namespace Kratos
