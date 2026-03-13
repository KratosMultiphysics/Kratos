//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <map>
#include <string>
#include <variant>

// External includes

// Project includes
#include "containers/flags.h"
#include "containers/nd_data.h"
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/global_variables.h"
#include "includes/io.h"
#include "includes/model_part.h"
#include "tensor_adaptors/tensor_adaptor.h"
#include "utilities/xml_utilities/xml_data_element_wrapper.h"

namespace Kratos {

/**
 * @class VtuOutput
 * @brief Handles VTU (Visualization Toolkit Unstructured grid) output for Kratos ModelParts.
 * @author Suneth Warnakulasuriya
 * @ingroup KratosCore
 *
 * This class provides functionality to export simulation data from a Kratos ModelPart to VTU files,
 * supporting both ASCII, BINARY, RAW and COMPRESSED_RAW with zlib compression. It allows users to register @ref Variable , @ref Flags ,
 * and @ref TensorAdaptor for output, and supports writing data on nodes, elements, conditions, and integration points.
 * The output can be configured to include submodel parts and supports parallel execution (MPI).
 *
 * @section vtu_output_usage Usage
 * This output put process can be used at any point in the simulation (static or dynamic). But, upon construction of an object of
 * @ref VtuOutput, all the variables, flags, integration point variables and tensor adaptors should be added
 * using @ref AddFlag, @ref AddVariable, @ref AddIntegrationPointVariable or @ref AddTensorAdaptor . Once the
 * @ref PrintOutput is called, no more data fields can be added because Vtk library does not support varying fields.
 *
 * But, already added @ref TensorAdaptor objects may be replaced with new @ref TensorAdaptor objects after calling @ref PrintOutput
 * by using @ref ReplaceTensorAdaptor method.
 *
 * A convenience method called @ref EmplaceTensorAdaptor is there to add the @ref TensorAdaptor if it is not existing and if the call is
 * before the first call of @ref PrintOutput. Otherwise, it will try replacing an existing @ref TensorAdaptor.
 *
 */
class KRATOS_API(KRATOS_CORE) VtuOutput : public IO
{
public:
    ///@name Type definitions
    ///@{

    // Index type definition
    using IndexType = std::size_t;

    // Variant defining supported container (i.e. @ref PointerVectorSet ) type pointers.
    using SupportedContainerPointerType = std::variant<
                                            ModelPart::NodesContainerType::Pointer,
                                            ModelPart::ConditionsContainerType::Pointer,
                                            ModelPart::ElementsContainerType::Pointer
                                        >;

    // Variant defining supported cell type @ref PointerVectorSte which can be visualized with cell data with vtk library.
    using CellContainerPointerType = std::variant<
                                        ModelPart::ConditionsContainerType::Pointer,
                                        ModelPart::ElementsContainerType::Pointer
                                    >;

    // Variant defining all the @ref Variable types which are supported to be visualized at point or cell data with vtk library.
    using SupportedVariablePointerType = std::variant<
                                                Variable<int> const *,
                                                Variable<double> const *,
                                                Variable<array_1d<double, 3>> const *,
                                                Variable<array_1d<double, 4>> const *,
                                                Variable<array_1d<double, 6>> const *,
                                                Variable<array_1d<double, 9>> const *,
                                                Variable<Vector> const *,
                                                Variable<Matrix> const *
                                            >;

    // Variant definining all the supported @ref TensorAdaptor which can be used to pass data fields to point or cell data using the vtk library.
    using SupportedTensorAdaptorPointerType = std::variant<
                                                    TensorAdaptor<bool>::Pointer,
                                                    TensorAdaptor<int>::Pointer,
                                                    TensorAdaptor<double>::Pointer
                                                >;


    // A list of maps of name and a @p T for different locations ( @ref Globals::DataLocation ) which are used for output at later time.
    template<class T>
    using DataList = std::array<std::map<std::string, T>, static_cast<std::size_t>(Kratos::Globals::DataLocation::NumberOfDataLocations)>;

    KRATOS_CLASS_POINTER_DEFINITION(VtuOutput);

    ///@}
    ///@name Public enums
    ///@{

    /// Enumerations for the output writer format.
    enum WriterFormat
    {
        ASCII,                          ///< ASCII format.
        BINARY,                         ///< Binary format.
        RAW,                            ///< Raw format. All data is appended to one stream.
        COMPRESSED_RAW                  ///< Data is first compressed with zlib and the appended to one stream.
    };

    struct UnstructuredGridData
    {
        bool                                    UsePointsForDataFieldOutput;                // If true, then the mpPoints represents a container in the model part, otherwise, mpPoints refers to a container which is on the fly generated.
        ModelPart *                             mpModelPart;                                // Model part associated with the unstructured grid data.
        ModelPart::NodesContainerType::Pointer  mpPoints;                                   // Points to be used in the unstructured grid.
        std::optional<CellContainerPointerType> mpCells;                                    // Cells to be used in the unstructured grid.
        std::map<std::string, SupportedTensorAdaptorPointerType> mMapOfPointTensorAdaptors; // Point data tensor adaptors.
        std::map<std::string, SupportedTensorAdaptorPointerType> mMapOfCellTensorAdaptors;  // Cell data tensor adaptors.
    };

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Vtu Output
     *
     * @param rModelPart                Model part to be used in the vtu output.
     * @param IsInitialConfiguration    Which nodal positions to be written out.
     * @param OutputFormat              The format of the output.
     * @param Precision                 Precision of the double values to be used when writing the doubles as ASCII.
     * @param OutputSubModelParts       To consider all the submodel parts recursively for output.
     * @param WriteIds                  To write ids under the name "KRATOS_ID" to be visualized.
     * @param EchoLevel                 Echo level for to print information.
     */
    VtuOutput(
        ModelPart& rModelPart,
        const bool IsInitialConfiguration = true,
        const WriterFormat OutputFormat = WriterFormat::COMPRESSED_RAW,
        const IndexType Precision = 9,
        const bool OutputSubModelParts = false,
        const bool WriteIds = false,
        const IndexType EchoLevel = 0);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Adds a flag to the VTU output.
     *
     * This method registers a flag @p rFlagVariable to be included in the VTU output file.
     * The flag will be associated with the specified data location (e.g., node, element, condition, refer @ref Globals::DataLocation )
     * given by @p DataLocation . This allows the Flag's values to be written to the output file
     * when @ref PrintOutput is called with the name @p rFlagName.
     *
     * @throws std::runtime_error if @ref AddFlag is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p rFlagName more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p rFlagName in one of the compatible data locations with @p DataLocation .
     * @throws std::runtime_error if @p DataLocation is not NodeNonHistorical, Condition, or Element.
     *
     * @param rFlagName The name of the flag to be added.
     * @param rFlagVariable The flag variable to be added.
     * @param DataLocation The location where the flag data is stored (e.g., NodeNonHistorical, Condition, Element).
     *
     * @see @ref FlagsTensorAdaptor
     */
    void AddFlag(
        const std::string& rFlagName,
        const Flags& rFlagVariable,
        const Globals::DataLocation& DataLocation);

    /**
     * @brief Adds a variable to the output for a specified data location.
     *
     * Registers the given @p rVariable to be included in the VTU output at the specified
     * @p DataLocation (e.g., NodeHistorical, NodeNonHistorical, Condition, or Element, refer @ref Globals::DataLocation).
     * This allows the variable's values to be written to the output file when @ref PrintOutput is called.
     *
     * @throws std::runtime_error if @ref AddVariable is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p rVariable more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p rVariable in one of the compatible data locations with @p DataLocation .
     * @throws std::runtime_error if @p DataLocation is not NodeHistorical, NodeNonHistorical, Condition, or Element.
     *
     * @tparam TDataType        Data type of the variable.
     * @param rVariable         Variable to be added.
     * @param DataLocation      The location in the model where the variable is defined (e.g., NodeHistorical, NodeNonHistorical, Condition, or Element).
     *
     * @see @ref VariableTensorAdaptor
     * @see @ref HistoricalVariableTensorAdaptor
     */
    template<class TDataType>
    void AddVariable(
        const Variable<TDataType>& rVariable,
        const Globals::DataLocation& DataLocation);

    /**
     * @brief Adds a variable to be output at integration points.
     *
     * This method registers a @p rVariable for output at the integration points of the mesh elements or conditions
     * specified by @p DataLocation  (refer @ref Globals::DataLocation ). This allows the entities variable's integration point values to
     * be written to the output file when @ref PrintOutput is called.
     *
     * @throws std::runtime_error if @ref AddIntegrationPointVariable is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p rVariable more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p rVariable in one of the compatible data locations with @p DataLocation .
     * @throws std::runtime_error if @p DataLocation is not Condition or Element.
     *
     * @tparam TDataType        Data type of the variable.
     * @param rVariable         Variable to be added.
     * @param DataLocation      Specifies the location (e.g., Condition or Element) where the variable is used to calculate integration point values.
     */
    template<class TDataType>
    void AddIntegrationPointVariable(
        const Variable<TDataType>& rVariable,
        const Globals::DataLocation& DataLocation);

    /**
     * @brief Adds a tensor adaptor to the output.
     *
     * Registers a @p pTensorAdaptor with the specified @p TensorAdaptorName for the the vtu output. this does not
     * copy the tensor adaptor. When @ref PrintOutput is called, the values of the @p pTensorAdaptor at this point will
     * be written to vtu.
     *
     * @throws std::runtime_error if @ref AddTensorAdaptor is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p rTensorAdaptorName more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p rTensorAdaptorName in one of the compatible data locations represented by the container of @p pTensorAdaptor .
     * @throws std::runtime_error if @p pTensorAdaptor does not corresponds to any of the containers which is written by this @ref VtuOutput .
     *
     * @tparam TTensorAdaptorPointerType    Pointer type of the tensor adaptor.
     * @param rTensorAdaptorName            The name to associate with the tensor adaptor.
     * @param pTensorAdaptor                Pointer to the tensor adaptor to be added.
     */
    template<class TTensorAdaptorPointerType>
    void AddTensorAdaptor(
        const std::string& rTensorAdaptorName,
        TTensorAdaptorPointerType pTensorAdaptor);

    /**
     * @brief Updates the tensor adaptor associated with the given name.
     *
     * This method assigns a new tensor adaptor to the specified name, allowing
     * for dynamic changes to the tensor adaptor used in the @ref VtuOutput process.
     * This also copies the internal data of @p pTensorAdaptor to vtu output's internal data.
     *
     * @throws std::runtime_error if there is no field name same as @p rTensorAdaptorName in one of the compatible data locations represented by the container of @p pTensorAdaptor .
     * @throws std::runtime_error if @p pTensorAdaptor does not corresponds to any of the containers which is written by this @ref VtuOutput .
     *
     * @tparam TTensorAdaptorPointerType    Pointer type of the tensor adaptor.
     * @param rTensorAdaptorName            The name identifying the tensor adaptor to update.
     * @param pTensorAdaptor                Pointer to the new tensor adaptor to be associated.
     */
    template<class TTensorAdaptorPointerType>
    void ReplaceTensorAdaptor(
        const std::string& rTensorAdaptorName,
        TTensorAdaptorPointerType pTensorAdaptor);

    /**
     * @brief Inserts a @ref TensorAdaptor into the internal storage with the specified name.
     *
     * This method:
     *  - If PrintOutput is not called at all, then calls @ref AddTensorAdaptor
     *  - If PrintOutput is at least once called, then calls @ref ReplaceTensorAdaptor
     *
     * @tparam TTensorAdaptorPointerType    Pointer type of the tensor adaptor.
     * @param rTensorAdaptorName            The name to associate with the @ref TensorAdaptor.
     * @param pTensorAdaptor                Pointer to the @ref TensorAdaptor to be stored.
     */
    template<class TTensorAdaptorPointerType>
    void EmplaceTensorAdaptor(
        const std::string& rTensorAdaptorName,
        TTensorAdaptorPointerType pTensorAdaptor);

    /**
    * @brief Returns the model part.
    * @return The constant reference to the model part.
    */
    const ModelPart& GetModelPart() const;

    /**
     * @brief Get the list of output containers
     * @details This method returns containers of the model part which are used when writing fields.
     * @return std::vector<SupportedContainerPointerType>   List of container which are used in when writing fields.
     */
    std::vector<SupportedContainerPointerType> GetOutputContainerList() const;

    /**
     * @brief Prints the vtu output data to a file with the specified filename prefix.
     *
     * Will create the folder with the provided @p rOutputFilenamePrefix . In this folder, following files will be created.
     * - One vtu file per timestep, per model part's container (Including sub model parts if @p OutputSubModelParts is true) (in MPI, per rank as well).
     * - If it find any gauss point fields, then there will be a vtu file per time step, per model part's container (in MPI, per rank as well).
     * - If this is called in MPI, then this will create one pvtu file per timestep, per model part's container (will be having information of all the corresponding rank vtu files)
     * - If this is called in MPI, then this will create one pvtu file per timestep, per model part's gauss point container (will be having information of all the corresponding rank vtu files)
     *
     * Finally it will create one pvd ( @p rOutputFilenamePrefix .pvd ) file for all the model parts, all the timesteps all the gauss point container linking
     * - In MPI the pvtu files created.
     * - In shared memory parallelism the vtu files created.
     *
     * @param rOutputFilenamePrefix The prefix to use for the output filename.
     */
    void PrintOutput(const std::string& rOutputFilenamePrefix);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    bool mIsPVDFileHeaderWritten;

    ModelPart& mrModelPart;

    const Globals::Configuration mConfiguration;

    const IndexType mEchoLevel;

    const WriterFormat mOutputFormat;

    const IndexType mPrecision;

    DataList<Flags const *> mFlags;

    DataList<SupportedVariablePointerType> mVariables;

    DataList<SupportedVariablePointerType> mIntegrationPointVariables;

    std::vector<UnstructuredGridData> mUnstructuredGridDataList;

    ///@}
    ///@name Private operations
    ///@{

    template<class TXmlElementDataWrapperCreateFunctor, class TXmlElementDataWrapperAppendFunctor>
    void WriteData(
        std::vector<std::pair<std::string, std::string>>& rPVDFileNameInfo,
        TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor,
        TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor,
        UnstructuredGridData& rUnstructuredGridData,
        const std::string& rOutputPrefix,
        const IndexType Step) const;

    template<class TXmlElementDataWrapperCreateFunctor, class TXmlElementDataWrapperAppendFunctor>
    std::pair<std::string, std::string> WriteUnstructuredGridData(
        TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor,
        TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor,
        UnstructuredGridData& rUnstructuredGridData,
        const std::string& rOutputPrefix,
        const IndexType Step) const;

    template<class TXmlElementDataWrapperCreateFunctor, class TXmlElementDataWrapperAppendFunctor>
    std::pair<std::string, std::string> WriteIntegrationPointData(
        TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor,
        TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor,
        UnstructuredGridData& rUnstructuredGridData,
        const std::string& rOutputPrefix,
        const IndexType Step) const;

    ///@}
};
} // namespace Kratos
