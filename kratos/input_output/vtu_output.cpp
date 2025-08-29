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

// System includes
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <filesystem>

// External includes

// Project includes
#include "includes/data_communicator.h"
#include "input_output/base_64_encoded_output.h"
#include "input_output/vtk_definitions.h"
#include "tensor_adaptors/flags_tensor_adaptor.h"
#include "tensor_adaptors/historical_variable_tensor_adaptor.h"
#include "tensor_adaptors/node_position_tensor_adaptor.h"
#include "tensor_adaptors/variable_tensor_adaptor.h"
#include "utilities/container_io_utils.h"
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/string_utilities.h"
#include "utilities/xml_utilities/xml_ascii_nd_data_element.h"
#include "utilities/xml_utilities/xml_base64_binary_nd_data_element.h"
#include "utilities/xml_utilities/xml_elements_array.h"
#include "utilities/variable_utils.h"
#include "includes/variables.h"

// Include base h
#include "vtu_output.h"
namespace Kratos
{

namespace
{

DenseVector<unsigned int> GetSynchronizedShape(
    const DenseVector<unsigned int>& rShape,
    const DataCommunicator& rDataCommunicator)
{
    std::vector<unsigned int> local_shape(rShape.begin(), rShape.end());
    const auto& max_shape = rDataCommunicator.MaxAll(local_shape);
    DenseVector<unsigned int> result(rShape.size());
    std::copy(max_shape.begin() + 1, max_shape.end(), result.begin() + 1);
    result[0] = rShape[0];
    return result;
}

struct XmlAsciiNDDataElementWrapper
{
    IndexType mPrecision;

    DataCommunicator const * mpDataCommunicator;

    template<class NDDataTypePointer>
    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataTypePointer pNDData)
    {
        using data_type = typename BareType<decltype(*pNDData)>::DataType;
        return Kratos::make_shared<XmlAsciiNDDataElement<data_type>>(
            rDataArrayName, pNDData,
            GetSynchronizedShape(pNDData->Shape(), *mpDataCommunicator), mPrecision);
    }
};

struct XmlBase64BinaryNDDataElementWrapper
{
    DataCommunicator const * mpDataCommunicator;

    template<class NDDataTypePointer>
    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataTypePointer pNDData)
    {
        using data_type = typename BareType<decltype(*pNDData)>::DataType;
        return Kratos::make_shared<XmlBase64BinaryNDDataElement<data_type>>(
            rDataArrayName, pNDData,
            GetSynchronizedShape(pNDData->Shape(), *mpDataCommunicator));
    }
};

std::string GetEndianness()
{
    int i = 0x0001;

    if (*reinterpret_cast<char*>(&i) != 0) {
        return "LittleEndian";
    } else {
        return "BigEndian";
    }
}

template<class... TLists>
void CheckDataArrayName(
    const std::string& rName,
    TLists&... rLists)
{
    const bool found_existing_name = !(... && (rLists.find(rName) == rLists.end()));
    KRATOS_ERROR_IF(found_existing_name)
        << "Found an existing data array with the same name = \"" << rName << "\".\n";
}

template <class TContainerType>
NDData<int>::Pointer GetOffsets(const TContainerType& rContainer)
{
    auto p_offsets = Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, rContainer.size()));
    auto span = p_offsets->ViewData();

    IndexType total_offset = 0;
    auto data_itr = span.begin();
    for (const auto& r_entity : rContainer) {
        total_offset += r_entity.GetGeometry().size();
        *(data_itr++) = total_offset;
    }

    return p_offsets;
}

template<class TContainerType>
NDData<unsigned char>::Pointer GetGeometryTypes(const TContainerType& rContainer)
{
    auto p_geometry_types = Kratos::make_shared<NDData<unsigned char>>(DenseVector<unsigned int>(1, rContainer.size()));
    auto span = p_geometry_types->ViewData();

    IndexPartition<IndexType>(rContainer.size()).for_each([&span, &rContainer](const IndexType Index) {
        const auto p_itr = VtkDefinitions::KratosVtkGeometryTypes.find((rContainer.begin() + Index)->GetGeometry().GetGeometryType());
        if (p_itr != VtkDefinitions::KratosVtkGeometryTypes.end()) {
            *(span.begin() + Index) = static_cast<unsigned char>(p_itr->second);
        } else {
            KRATOS_ERROR << "Element with id " << (rContainer.begin() + Index)->Id() << " has unsupported geometry.";
        }
    });
    return p_geometry_types;
}

template<class TContainerType>
NDData<int>::Pointer GetConnectivities(
    const NDData<int>& rOffsets,
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    const auto offsets_span = rOffsets.ViewData();

    auto p_connectivities = Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, offsets_span.back()));
    auto connectivities_span = p_connectivities->ViewData();

    IndexPartition<IndexType>(rContainer.size()).for_each([&connectivities_span, &offsets_span, &rContainer, &rKratosVtuIndicesMap](const IndexType Index) {
        const auto& r_geometry = (rContainer.begin() + Index)->GetGeometry();
        auto entity_data_begin_itr = connectivities_span.begin() + offsets_span[Index] - r_geometry.size();

        for (const auto& r_node : r_geometry) {
            const auto p_itr = rKratosVtuIndicesMap.find(r_node.Id());
            if (p_itr != rKratosVtuIndicesMap.end()) {
                *(entity_data_begin_itr++) = p_itr->second;
            } else {
                KRATOS_ERROR << "Node with id " << r_node.Id() << " not found in nodes list.";
            }
        }
    });
    return p_connectivities;
}

void CopyAttributes(
    const XmlElement& rSource,
    XmlElement& rDestination)
{
    for (const auto& [attribute, value] : rSource.GetAttributes()) {
        rDestination.AddAttribute(attribute, value);
    }
}

template<class TXmlDataElementWrapper>
XmlElement::Pointer CreateCellsXmlElement(
    std::optional<VtuOutput::CellContainerPointerType> pCells,
    const std::unordered_map<IndexType, IndexType>& rIndicesMap,
    TXmlDataElementWrapper& rXmlDataElementWrapper)
{
    auto p_cells_xml_element = Kratos::make_shared<XmlElementsArray>("Cells");

    if (pCells.has_value()) {
        std::visit([&rIndicesMap, &p_cells_xml_element, &rXmlDataElementWrapper](auto p_container) {
            auto p_offsets = GetOffsets(*p_container);
            p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("connectivity", GetConnectivities(*p_offsets, *p_container, rIndicesMap)));
            p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("offsets", p_offsets));
            p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("types", GetGeometryTypes(*p_container)));
        }, pCells.value());
    } else {
        // this means non of the ranks has any cells to output.
        // so we output blank datasets
    }

    return p_cells_xml_element;
}

template<class TTensorAdaptorType, class TContainerPointerType, class TMapType, class TXmlDataElementWrapper, class... TArgs>
void AddFieldsFromTensorAdaptor(
    XmlElementsArray& rXmlElement,
    TContainerPointerType pContainer,
    const TMapType& rMap,
    TXmlDataElementWrapper& rXmlDataElementWrapper,
    TArgs&&... rArgs)
{
    KRATOS_TRY

    for (const auto& r_pair : rMap) {
        using data_type = BareType<decltype(r_pair.second)>;
        if constexpr(std::is_same_v<data_type, Flags>) {
            // the map is of type flags
            TTensorAdaptorType tensor_adaptor(pContainer, *r_pair.second, rArgs...);
            tensor_adaptor.CollectData();
            rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, tensor_adaptor.pGetStorage()));
        } else {
            std::visit([&](const auto p_variable) {
                using primitive_data_type = typename DataTypeTraits<typename BareType<decltype(*p_variable)>::Type>::PrimitiveType;
                if constexpr(std::is_same_v<primitive_data_type, int>) {
                    // we only support int variable, and there are no TensorAdaptors to
                    // collect data from integer variables, we do it manually here.
                    auto p_nd_data = Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, pContainer->size()));
                    const auto& shape = p_nd_data->Shape();
                    if constexpr (std::is_same_v<TTensorAdaptorType, HistoricalVariableTensorAdaptor>) {
                        ContainerIOUtils::CopyToContiguousArray<int>(
                            *pContainer, p_nd_data->ViewData(),
                            shape.data().begin(), shape.data().begin() + 1,
                            [p_variable](int& rValue, const Node& rNode) {
                                rValue = rNode.FastGetSolutionStepValue(*p_variable);
                            });
                    } else {
                        ContainerIOUtils::CopyToContiguousArray<int>(
                            *pContainer, p_nd_data->ViewData(),
                            shape.data().begin(), shape.data().begin() + 1,
                            [p_variable](int& rValue, const auto& rNode) {
                                rValue = rNode.GetValue(*p_variable);
                            });
                    }
                    rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, p_nd_data));
                } else {
                    TTensorAdaptorType tensor_adaptor(pContainer, p_variable, rArgs...);
                    tensor_adaptor.CollectData();
                    rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, tensor_adaptor.pGetStorage()));
                }
            }, r_pair.second);
        }
    }

    KRATOS_CATCH("");
}

template<class TXmlDataElementWrapper>
void AddFields(
    XmlElementsArray& rXmlElement,
    const std::map<std::string, VtuOutput::FieldPointerType>& rMap,
    TXmlDataElementWrapper& rXmlDataElementWrapper)
{
    for (const auto& r_pair : rMap) {
        std::visit([&rXmlElement, &r_pair, &rXmlDataElementWrapper](auto pNDData) {
            rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, pNDData));
        }, r_pair.second);
    }
}

void FillIndicesMap(
    VtuOutput::IndicesMap& rIndicesMap,
    const ModelPart::NodesContainerType& rNodes)
{
    rIndicesMap.clear();
    IndexType vtu_index = 0;
    for (const auto& r_node : rNodes) {
        rIndicesMap[r_node.Id()] = vtu_index++;
    }
}

template<class TEntityContainerType>
ModelPart::NodesContainerType::Pointer GetNodesContainer(TEntityContainerType& rContainer)
{
    std::vector<ModelPart::NodeType::Pointer> temp_nodes;
    temp_nodes.reserve(rContainer.size() * 20);
    for (auto& r_entity : rContainer) {
        auto& r_geometry = r_entity.GetGeometry();
        for (auto p_itr = r_geometry.ptr_begin(); p_itr != r_geometry.ptr_end(); ++p_itr) {
            temp_nodes.push_back(*p_itr);
        }
    }

    auto p_nodes_container = Kratos::make_shared<ModelPart::NodesContainerType>();
    p_nodes_container->insert(temp_nodes.begin(), temp_nodes.end());
    return p_nodes_container;
}

void AddModelPartData(
    std::vector<VtuOutput::ModelPartData>& rOutput,
    ModelPart& rModelPart)
{
    const bool has_elements   = rModelPart.GetCommunicator().GlobalNumberOfElements() > 0;
    const bool has_conditions = rModelPart.GetCommunicator().GlobalNumberOfConditions() > 0;

    if (has_elements) {
        // Model part has elements. Hence add a separate output
        // for elements.

        // indices map to be filled
        auto p_indices_map = Kratos::make_shared<VtuOutput::IndicesMap>();

        // now check if it has proper nodes
        if (rModelPart.NumberOfNodes() > 0) {
            FillIndicesMap(*p_indices_map, rModelPart.Nodes());
            VtuOutput::ModelPartData model_part_data{true, &rModelPart, rModelPart.pNodes(), rModelPart.pElements(), p_indices_map};
            rOutput.push_back(model_part_data);
        } else {
            // create the nodes container.
            auto p_nodes = GetNodesContainer(rModelPart.Elements());
            FillIndicesMap(*p_indices_map, *p_nodes);
            VtuOutput::ModelPartData model_part_data{false, &rModelPart, p_nodes, rModelPart.pElements(), p_indices_map};
            rOutput.push_back(model_part_data);
        }
    }

    if (has_conditions) {
        // Model part has conditions. Hence add a separate output
        // for conditions.

        // indices map to be filled
        auto p_indices_map = Kratos::make_shared<VtuOutput::IndicesMap>();

        if (!has_elements && rModelPart.NumberOfNodes() > 0) {
            FillIndicesMap(*p_indices_map, rModelPart.Nodes());
            VtuOutput::ModelPartData model_part_data{true, &rModelPart, rModelPart.pNodes(), rModelPart.pConditions(), p_indices_map};
            rOutput.push_back(model_part_data);
        } else {
            // either this model part also contains elements, or it does not
            // contain nodes. In either case, the nodes list given by the rModelPart
            // does not reflect the actual nodes used by the conditions.
            // In order to avoid writing nodes, which are not used by the conditions,
            // this will use a new nodes list.

            auto p_nodes = GetNodesContainer(rModelPart.Conditions());
            FillIndicesMap(*p_indices_map, *p_nodes);
            VtuOutput::ModelPartData model_part_data{false, &rModelPart, p_nodes, rModelPart.pConditions(), p_indices_map};
            rOutput.push_back(model_part_data);
        }
    }

    if (!has_elements && !has_conditions) {
        // Model part does not have either conditions or elements.
        // Hence, only adding the nodes.
        VtuOutput::ModelPartData model_part_data{true, &rModelPart, rModelPart.pNodes(), std::nullopt, Kratos::make_shared<VtuOutput::IndicesMap>()};
        rOutput.push_back(model_part_data);
    }
}

void AddModelPartDataRecursively(
    std::vector<VtuOutput::ModelPartData>& rOutput,
    ModelPart& rModelPart)
{
    // add the current model part data to the output.
    AddModelPartData(rOutput, rModelPart);

    // now recursively add all the sub model part data.
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        AddModelPartDataRecursively(rOutput, r_sub_model_part);
    }
}

template<class... TArgs>
std::string GetName(
    const std::variant<TArgs...>& rArg)
{
    return std::visit([](const auto& rArg) { return rArg->Name(); }, rArg);
}

template<class TKeyType, class T>
const T& GetUnorderedMapValue(
    const TKeyType& rKey,
    const std::unordered_map<TKeyType, T>& rMap)
{
    auto itr = rMap.find(rKey);
    if (itr != rMap.end()) {
        return itr->second;
    } else {
        static const T default_value{};
        return default_value;
    }
}

template<class TMapType, class TCellPointerType>
const TMapType& GetContainerMap(
    const std::unordered_map<Globals::DataLocation, TMapType>& rMap,
    TCellPointerType pCellPointer)
{
    return std::visit([&rMap](auto pContainer) -> const TMapType& {
        using container_type = BareType<decltype(*pContainer)>;
        if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType>) {
            return GetUnorderedMapValue(Globals::DataLocation::Condition, rMap);
        } else if constexpr(std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
            return GetUnorderedMapValue(Globals::DataLocation::Element, rMap);
        } else {
            KRATOS_ERROR << "Unsupported container type.";
            return GetUnorderedMapValue(Globals::DataLocation::Element, rMap);
        }
    }, pCellPointer);
}

std::string WritePartitionedUnstructuredGridData(
    XmlElementsArray& rPointDataElement,
    XmlElementsArray& rCellDataElement,
    const std::string& rOutputVtuFileName,
    const DataCommunicator& rDataCommunicator)
{
    // get list of file names
    std::stringstream list_of_file_names;

    // TODO: May be we want to check a rank which has some entities (not empty ranks)
    //       Then write on that rank.
    const int writing_rank = 0;

    list_of_file_names << rOutputVtuFileName << "\n";
    if (rDataCommunicator.Rank() == writing_rank) {
        for (int rank = 0; rank < rDataCommunicator.Size(); ++rank) {
            if (rank != writing_rank) {
                std::string msg;
                rDataCommunicator.Recv(msg, rank);
                list_of_file_names << msg;
            }
        }
    } else {
        rDataCommunicator.Send(list_of_file_names.str(), writing_rank);
    }

    rDataCommunicator.Barrier();

    // remove the rank from the rOutputVtuFileName.
    const auto& p_vtu_file_name = rOutputVtuFileName.substr(0, rOutputVtuFileName.rfind("_"))  + ".pvtu";

    if (rDataCommunicator.Rank() == writing_rank) {
        // create the pvtu file
        XmlElementsArray p_vtu_file_element("VTKFile");
        p_vtu_file_element.AddAttribute("type", "PUnstructuredGrid");
        p_vtu_file_element.AddAttribute("version", "0.1");
        p_vtu_file_element.AddAttribute("byte_order", GetEndianness());

        // create the unstructured grid
        auto p_unstructured_grid_element = Kratos::make_shared<XmlElementsArray>("PUnstructuredGrid");
        p_unstructured_grid_element->AddAttribute("GhostLevel", "0");
        p_vtu_file_element.AddElement(p_unstructured_grid_element);

        // ppoints_element
        auto p_points_element = Kratos::make_shared<XmlElementsArray>("PPoints");
        p_unstructured_grid_element->AddElement(p_points_element);

        // position element
        auto p_position_element = Kratos::make_shared<XmlElementsArray>("PDataArray");
        p_position_element->AddAttribute("type", "Float64");
        p_position_element->AddAttribute("Name", "Position");
        p_position_element->AddAttribute("NumberOfComponents", "3");
        p_points_element->AddElement(p_position_element);

        // pcells element
        auto p_cells_element = Kratos::make_shared<XmlElementsArray>("PCells");
        p_unstructured_grid_element->AddElement(p_cells_element);

        // connectivity element
        auto p_connectivity_element = Kratos::make_shared<XmlElementsArray>("PDataArray");
        p_connectivity_element->AddAttribute("type", "Int32");
        p_connectivity_element->AddAttribute("Name", "connectivity");
        p_connectivity_element->AddAttribute("NumberOfComponents", "1");
        p_cells_element->AddElement(p_connectivity_element);

        // offsets element
        auto p_offsets_element = Kratos::make_shared<XmlElementsArray>("PDataArray");
        p_offsets_element->AddAttribute("type", "Int32");
        p_offsets_element->AddAttribute("Name", "offsets");
        p_offsets_element->AddAttribute("NumberOfComponents", "1");
        p_cells_element->AddElement(p_offsets_element);

        // types element
        auto p_types_element = Kratos::make_shared<XmlElementsArray>("PDataArray");
        p_types_element->AddAttribute("type", "Int32");
        p_types_element->AddAttribute("Name", "types");
        p_types_element->AddAttribute("NumberOfComponents", "1");
        p_cells_element->AddElement(p_types_element);

        // ppoint data element
        auto p_point_data_element = Kratos::make_shared<XmlElementsArray>("PPointData");
        p_unstructured_grid_element->AddElement(p_point_data_element);

        // now add the point data fields
        for (const auto& p_element : rPointDataElement.GetElements()) {
            auto p_current_element = Kratos::make_shared<XmlElementsArray>("PDataArray");
            CopyAttributes(*p_element, *p_current_element);
            p_point_data_element->AddElement(p_current_element);
        }

        // pcell data element
        auto p_cell_data_element = Kratos::make_shared<XmlElementsArray>("PCellData");
        p_unstructured_grid_element->AddElement(p_cell_data_element);

        // now add the cell data fields
        for (const auto& p_element : rCellDataElement.GetElements()) {
            auto p_current_element = Kratos::make_shared<XmlElementsArray>("PDataArray");
            CopyAttributes(*p_element, *p_current_element);
            p_cell_data_element->AddElement(p_current_element);
        }

        // now add the piece elements
        const auto& r_file_names = StringUtilities::SplitStringByDelimiter(list_of_file_names.str(), '\n');
        for (const auto& r_file_name : r_file_names) {
            auto piece = Kratos::make_shared<XmlElementsArray>("Piece");
            // since we are writing to the same folder the pvtu files
            piece->AddAttribute(
                "Source", std::filesystem::relative(
                              std::filesystem::absolute(r_file_name),
                              std::filesystem::absolute(p_vtu_file_name).parent_path())
                              .string());
            p_unstructured_grid_element->AddElement(piece);
        }

        // writing to file
        std::ofstream output_file;
        output_file.open(p_vtu_file_name, std::ios::out | std::ios::trunc);
        p_vtu_file_element.Write(output_file);
    }

    return p_vtu_file_name;
}

template<class TShapeType>
NDData<double>::Pointer GetNodalNDData(
    ModelPart::NodesContainerType& rNodes,
    const TShapeType& rDataShape)
{
    // construct the nd_data shape
    DenseVector<unsigned int> nd_data_shape(rDataShape.size() + 1);
    // vtu output always write data for local and ghost nodes.
    nd_data_shape[0] = rNodes.size();
    std::copy(rDataShape.begin(), rDataShape.end(), nd_data_shape.begin() + 1);

    auto p_nd_data = Kratos::make_shared<NDData<double>>(nd_data_shape);
    const auto nd_data_span = p_nd_data->ViewData();

    // now read the variable data back to the nd_data
    IndexPartition<IndexType>(nd_data_shape[0]).for_each([&rNodes, &nd_data_span](const auto Index) {
        const auto& r_values = (rNodes.begin() + Index)->GetValue(TENSOR_ADAPTOR_SYNC);
        const auto data_begin_index = Index * r_values.size();
        for (IndexType i = 0; i < r_values.size(); ++i) {
            nd_data_span[data_begin_index + i] = r_values[i];
        }
    });

    return p_nd_data;
}

template<class TMapType>
void PrintDataLocationData(
    std::ostream& rOStream,
    const std::string& rMapType,
    const std::unordered_map<Globals::DataLocation, TMapType>& rMap)
{
    rOStream << "List of " << rMapType << "s:";
    for (const auto& [data_location, map] : rMap) {
        switch (data_location) {
            case Globals::DataLocation::NodeHistorical:
                rOStream << "\n\tNode historical " << rMapType << "s:";
                break;
            case Globals::DataLocation::NodeNonHistorical:
                rOStream << "\n\tNode " << rMapType << "s:";
                break;
            case Globals::DataLocation::Condition:
                rOStream << "\n\tCondition " << rMapType << "s:";
                break;
            case Globals::DataLocation::Element:
                rOStream << "\n\tElement " << rMapType << "s:";
                break;
            default:
                rOStream << "\n\tUnsupported " << rMapType << "s:";
                break;
        }

        for ([[maybe_unused]]const auto& [name, variable] : map) {
            rOStream << "\n\t\t" << name;
        }
    }
}

} // namespace

VtuOutput::VtuOutput(
    ModelPart& rModelPart,
    const bool IsInitialConfiguration,
    const WriterFormat OutputFormat,
    const IndexType Precision,
    const bool OutputSubModelParts)
    : mrModelPart(rModelPart),
      mIsInitialConfiguration(IsInitialConfiguration),
      mOutputFormat(OutputFormat),
      mPrecision(Precision)
{
    if (OutputSubModelParts) {
        // Collect all model part data recursively.
        AddModelPartDataRecursively(mListOfModelPartData, rModelPart);
    } else {
        // Only collect the data from the Passed model part
        AddModelPartData(mListOfModelPartData, rModelPart);
    }
}

const ModelPart& VtuOutput::GetModelPart() const
{
    return mrModelPart;
}

void VtuOutput::AddFlag(
    const std::string& rFlagName,
    const Flags& rFlagVariable,
    Globals::DataLocation DataLocation)
{
    switch (DataLocation) {
        case Globals::DataLocation::NodeNonHistorical:
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            mFlags[DataLocation][rFlagName] = &rFlagVariable;
            break;
        default:
            KRATOS_ERROR << "Flags can be only added to NodeNonHistorical, Condition, and Element data locations.";
            break;
    }
}

void VtuOutput::AddVariable(
    SupportedVariablePointerType pVariable,
    Globals::DataLocation DataLocation)
{
    switch (DataLocation) {
        case Globals::DataLocation::NodeHistorical:
        case Globals::DataLocation::NodeNonHistorical:
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            mVariables[DataLocation][GetName(pVariable)] = pVariable;
            break;
        default:
            KRATOS_ERROR << "Variables can be only added to NodeHistorical, NodeNonHistorical, Condition, and Element data locations.";
            break;
    }
}

void VtuOutput::AddIntegrationPointVariable(
    SupportedVariablePointerType pVariable,
    Globals::DataLocation DataLocation)
{
    switch (DataLocation) {
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            mIntegrationPointVariables[DataLocation][GetName(pVariable)] = pVariable;
            break;
        default:
            KRATOS_ERROR << "Integration point variables can be only added to Condition, and Element data locations.";
            break;
    }
}

void VtuOutput::ClearFlags(Globals::DataLocation DataLocation)
{
    switch (DataLocation) {
        case Globals::DataLocation::NodeNonHistorical:
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            mFlags[DataLocation].clear();
            break;
        default:
            KRATOS_ERROR << "Flags can be only cleared on NodeNonHistorical, Condition, and Element data locations.";
            break;
    }
}

void VtuOutput::ClearVariables(Globals::DataLocation DataLocation)
{
    switch (DataLocation) {
        case Globals::DataLocation::NodeHistorical:
        case Globals::DataLocation::NodeNonHistorical:
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            mVariables[DataLocation].clear();
            break;
        default:
            KRATOS_ERROR << "Variables can be only cleared on NodeHistorical, NodeNonHistorical, Condition, and Element data locations.";
            break;
    }
}

void VtuOutput::ClearIntegrationPointVariables(Globals::DataLocation DataLocation)
{
    switch (DataLocation) {
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            mIntegrationPointVariables[DataLocation].clear();
            break;
        default:
            KRATOS_ERROR << "Integration point variables can be only cleared on Condition, and Element data locations.";
            break;
    }
}

void VtuOutput::ClearPointFields()
{
    for (auto& r_model_part_data : mListOfModelPartData) {
        r_model_part_data.mPointFields.clear();
    }
}

void VtuOutput::ClearCellFields()
{
    for (auto& r_model_part_data : mListOfModelPartData) {
        r_model_part_data.mCellFields.clear();
    }
}

void VtuOutput::AddContainerExpression(
    const std::string& rExpressionName,
    SupportedContainerExpressionPointerType pContainerExpression)
{
    std::visit([this, &rExpressionName](auto p_container_expression) {
        const auto& r_expression = p_container_expression->GetExpression();

        const auto number_of_data_components = r_expression.GetItemComponentCount();
        const auto& data_shape = r_expression.GetItemShape();

        auto& r_expression_container = p_container_expression->GetContainer();

        using expression_container_type = BareType<decltype(r_expression_container)>;

        for (IndexType i_model_part_data = 0; i_model_part_data < this->mListOfModelPartData.size(); ++i_model_part_data) {
            auto& r_model_part_data = this->mListOfModelPartData[i_model_part_data];

            // get the model part
            auto p_model_part = r_model_part_data.mpModelPart;

            // expressions are created from local mesh entities.

            if constexpr(std::is_same_v<expression_container_type, ModelPart::NodesContainerType>) {
                if (&r_expression_container == &p_model_part->GetCommunicator().LocalMesh().Nodes()) {
                    // since expressions are created from local mesh nodes
                    // we need to synchronize the values so that ghost mesh nodes
                    // will be correctly filled.

                    // first clear the TENSOR_ADAPTOR_SYNC.
                    VariableUtils().SetNonHistoricalVariableToZero(TENSOR_ADAPTOR_SYNC, p_model_part->Nodes());

                    // secondly fill in the local nodal values to the temporary variable TENSOR_ADAPTOR_SYNC
                    IndexPartition<IndexType>(r_expression_container.size()).for_each(Vector{number_of_data_components}, [&r_expression_container, &r_expression, number_of_data_components](const auto Index, auto& rTLS) {
                        const auto data_begin_index = Index * number_of_data_components;
                        auto& r_node = *(r_expression_container.begin() + Index);
                        for (IndexType i = 0; i < number_of_data_components; ++i) {
                            rTLS[i] = r_expression.Evaluate(Index, data_begin_index, i);
                        }
                        r_node.SetValue(TENSOR_ADAPTOR_SYNC, rTLS);
                    });

                    // the p_nd_data is now correctly filled. Add it to the
                    // map and then exit the for loop since, the given container expression
                    // is already found.
                    p_model_part->GetCommunicator().SynchronizeVariable(TENSOR_ADAPTOR_SYNC);
                    r_model_part_data.mPointFields[rExpressionName] = GetNodalNDData(*r_model_part_data.mpPoints, data_shape);
                    return true;
                }
            } else {
                if (r_model_part_data.mpCells.has_value()) {
                    return std::visit([this, &r_model_part_data, &rExpressionName, &r_expression, &r_expression_container, &data_shape, number_of_data_components](auto p_model_part_container) {
                        using model_part_container = BareType<decltype(*p_model_part_container)>;

                        if constexpr(std::is_same_v<expression_container_type, model_part_container>) {
                            if (&r_expression_container == &*p_model_part_container) {
                                // found a correct container. Just copy the data

                                // construct the nd_data shape
                                DenseVector<unsigned int> nd_data_shape(data_shape.size() + 1);
                                // vtu output always write data for local and ghost nodes.
                                nd_data_shape[0] = r_expression_container.size();
                                std::copy(data_shape.begin(), data_shape.end(), nd_data_shape.begin() + 1);

                                auto p_nd_data = Kratos::make_shared<NDData<double>>(nd_data_shape);
                                const auto nd_data_span = p_nd_data->ViewData();

                                IndexPartition<IndexType>(r_expression_container.size()).for_each([&nd_data_span, &r_expression, number_of_data_components](auto Index){
                                    const auto data_begin_index = Index * number_of_data_components;
                                    for (IndexType i = 0; i < number_of_data_components; ++i) {
                                        nd_data_span[data_begin_index + i] = r_expression.Evaluate(Index, data_begin_index, i);
                                    }
                                });

                                r_model_part_data.mCellFields[rExpressionName] = p_nd_data;
                                return true;
                            }
                        }
                        return false;
                    }, r_model_part_data.mpCells.value());
                }
            }
        }

        KRATOS_ERROR
            << "The container in the ContainerExpression is not referring to any of the containers "
            << "written by this VTU output [ container expression name = " << rExpressionName
            << ", tensor_adaptor = " << *p_container_expression << " ]\n"
            << *this;

    }, pContainerExpression);
}

void VtuOutput::AddTensorAdaptor(
    const std::string& rTensorAdaptorName,
    SupportedTensorAdaptorPointerType pTensorAdaptor,
    const bool Copy)
{
    std::visit([this, &rTensorAdaptorName, Copy](auto p_tensor_adaptor) {
        using tensor_adaptor_type = BareType<decltype(*p_tensor_adaptor)>;
        const auto shape = p_tensor_adaptor->Shape();
        const auto number_of_data_components = std::accumulate(shape.begin() + 1, shape.end(), 1u, std::multiplies<unsigned int>{});

        const bool ta_added = std::visit([this, &rTensorAdaptorName, &p_tensor_adaptor, &shape, number_of_data_components, Copy](auto p_ta_container){
            using ta_container_type = BareType<decltype(*p_ta_container)>;

            // first check to which model part output this tenor adaptor belongs to
            for (IndexType i = 0; i < this->mListOfModelPartData.size(); ++i) {
                auto& r_model_part_data = this->mListOfModelPartData[i];
                auto p_model_part = r_model_part_data.mpModelPart;

                if constexpr(std::is_same_v<ta_container_type, ModelPart::NodesContainerType>) {
                    // tensor adaptor is having a nodal container. Now we check if
                    // it is of normal model part nodes, local nodes, ghost nodes or interface nodes

                    auto& r_communicator = p_model_part->GetCommunicator();

                    KRATOS_ERROR_IF(&*(p_ta_container) == &r_communicator.GhostMesh().Nodes())
                        << "TensorAdaptors having containers referring to ghost nodes container is not allowed [ ghost nodes model part = "
                        << p_model_part->FullName() << ", tensor_adaptor = " << *p_tensor_adaptor << " ].\n";

                    KRATOS_ERROR_IF(&*(p_ta_container) == &r_communicator.InterfaceMesh().Nodes())
                        << "TensorAdaptors having containers referring to interface nodes container is not allowed [ interface nodes model part = "
                        << p_model_part->FullName() << ", tensor_adaptor = " << *p_tensor_adaptor << " ].\n";

                    if (&*(p_ta_container) == &p_model_part->Nodes()) {
                        // here we don't have to do anything. add the tensor adaptor's
                        // NDData
                        r_model_part_data.mPointFields[rTensorAdaptorName] = tensor_adaptor_type(*p_tensor_adaptor, Copy).pGetStorage();
                        return true;
                    } else if (&*(p_ta_container) == &p_model_part->GetCommunicator().LocalMesh().Nodes()) {
                        // the tensor adaptor is having a container referring to the local mesh nodes.
                        // now we have to construct the nodal values for all the nodes including ghost mesh
                        // nodes.

                        // first clear the TENSOR_ADAPTOR_SYNC.
                        VariableUtils().SetNonHistoricalVariableToZero(TENSOR_ADAPTOR_SYNC, p_model_part->Nodes());

                        auto span = p_tensor_adaptor->ViewData();

                        // secondly fill in the local nodal values to the temporary variable TENSOR_ADAPTOR_SYNC
                        IndexPartition<IndexType>(p_ta_container->size()).for_each(Vector{number_of_data_components}, [&p_ta_container, &span, number_of_data_components](const auto Index, auto& rTLS) {
                            const auto data_begin_index = Index * number_of_data_components;
                            std::copy(span.begin() + data_begin_index, span.begin() + data_begin_index + number_of_data_components, rTLS.begin());
                            (p_ta_container->begin() + Index)->SetValue(TENSOR_ADAPTOR_SYNC, rTLS);
                        });

                        // the p_nd_data is now correctly filled. Add it to the
                        // map and then exit the for loop since, the given container expression
                        // is already found.
                        std::vector<unsigned int> data_shape(shape.begin() + 1, shape.end());
                        p_model_part->GetCommunicator().SynchronizeVariable(TENSOR_ADAPTOR_SYNC);
                        r_model_part_data.mPointFields[rTensorAdaptorName] = GetNodalNDData(*r_model_part_data.mpPoints, data_shape);
                        return true;
                    }
                } else if constexpr(std::is_same_v<ta_container_type, ModelPart::ConditionsContainerType>) {
                    // ghost, local and model part containers are the same. so we check only the model part container
                    if (r_model_part_data.mpCells.has_value() &&
                        std::holds_alternative<ModelPart::ConditionsContainerType::Pointer>(r_model_part_data.mpCells.value()) &&
                        &*p_ta_container == &p_model_part->Conditions()) {
                            r_model_part_data.mCellFields[rTensorAdaptorName] = tensor_adaptor_type(*p_tensor_adaptor, Copy).pGetStorage();
                            return true;
                    }
                } else if constexpr(std::is_same_v<ta_container_type, ModelPart::ElementsContainerType>) {
                    // ghost, local and model part containers are the same. so we check only the model part container
                    if (r_model_part_data.mpCells.has_value() &&
                        std::holds_alternative<ModelPart::ElementsContainerType::Pointer>(r_model_part_data.mpCells.value()) &&
                        &*p_ta_container == &p_model_part->Elements()) {
                            r_model_part_data.mCellFields[rTensorAdaptorName] = tensor_adaptor_type(*p_tensor_adaptor, Copy).pGetStorage();
                            return true;
                    }
                }
            }

            // hasn't found a valid container.
            return false;
        }, p_tensor_adaptor->GetContainer());

        KRATOS_ERROR_IF_NOT(ta_added)
            << "The container in the TensorAdaptor is not referring to any of the containers "
            << "written by this VTU output [ tensor adaptor name = " << rTensorAdaptorName
            << ", tensor_adaptor = " << *p_tensor_adaptor << " ]\n"
            << *this;

    }, pTensorAdaptor);
}

template<class TXmlDataElementWrapper>
std::pair<std::string, std::string> VtuOutput::WriteUnstructuredGridData(
    const std::string& rOutputFileNamePrefix,
    ModelPartData& rModelPartData,
    TXmlDataElementWrapper& rXmlDataElementWrapper) const
{
    auto p_nodes = rModelPartData.mpPoints;

    // create the vtk file
    XmlElementsArray vtk_file_element("VTKFile");
    vtk_file_element.AddAttribute("type", "UnstructuredGrid");
    vtk_file_element.AddAttribute("version", "0.1");
    vtk_file_element.AddAttribute("byte_order", GetEndianness());

    // create the unstructured grid
    auto unstructured_grid_element = Kratos::make_shared<XmlElementsArray>("UnstructuredGrid");
    vtk_file_element.AddElement(unstructured_grid_element);

    // create the piece element
    auto piece_element = Kratos::make_shared<XmlElementsArray>("Piece");
    // adding number of points
    piece_element->AddAttribute("NumberOfPoints", std::to_string(p_nodes->size()));
    // adding number of cells
    piece_element->AddAttribute(
        "NumberOfCells",
        std::to_string(rModelPartData.mpCells.has_value()
                           ? std::visit([](auto v) { return v->size(); },
                                        rModelPartData.mpCells.value())
                           : 0u));
    unstructured_grid_element->AddElement(piece_element);

    // create the position element
    NodePositionTensorAdaptor node_position_tensor_adaptor(
        p_nodes, mIsInitialConfiguration ? Globals::Configuration::Initial
                                         : Globals::Configuration::Current);
    node_position_tensor_adaptor.CollectData();

    // create the points element
    auto points_element = Kratos::make_shared<XmlElementsArray>("Points");
    points_element->AddElement(rXmlDataElementWrapper.Get("Position", node_position_tensor_adaptor.pGetStorage()));
    piece_element->AddElement(points_element);

    // create the cells element
    auto cells_element = CreateCellsXmlElement(rModelPartData.mpCells, *rModelPartData.mpIndicesMap, rXmlDataElementWrapper);
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
    piece_element->AddElement(point_data_element);

    if (rModelPartData.UsePointsForDataFieldOutput) {
        // generate and add point field data
        AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeNonHistorical, mFlags), rXmlDataElementWrapper);
        AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeNonHistorical, mVariables), rXmlDataElementWrapper);
        AddFieldsFromTensorAdaptor<HistoricalVariableTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeHistorical, mVariables), rXmlDataElementWrapper);
        AddFields(*point_data_element, rModelPartData.mPointFields, rXmlDataElementWrapper);
    }

    // create cell data
    auto cell_data_element = Kratos::make_shared<XmlElementsArray>("CellData");
    piece_element->AddElement(cell_data_element);

    // generate and add cell field data
    if (rModelPartData.mpCells.has_value()) {
        std::visit([this, &rModelPartData, &cell_data_element, &rXmlDataElementWrapper](auto p_container){
            AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*cell_data_element, p_container, GetContainerMap(this->mFlags, rModelPartData.mpCells.value()), rXmlDataElementWrapper);
            AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*cell_data_element, p_container, GetContainerMap(this->mVariables, rModelPartData.mpCells.value()), rXmlDataElementWrapper);
        }, rModelPartData.mpCells.value());
    }
    AddFields(*cell_data_element, rModelPartData.mCellFields, rXmlDataElementWrapper);

    std::stringstream output_vtu_file_name;
    output_vtu_file_name << rOutputFileNamePrefix << "/" << rModelPartData.mpModelPart->FullName();

    // identify suffix with the entity type.
    std::string suffix;
    if (rModelPartData.mpCells.has_value()) {
        suffix = std::visit([](auto p_container) {
            using container_type = BareType<decltype(*p_container)>;
            return "_" + ModelPart::Container<container_type>::GetEntityName() + "s";
        }, rModelPartData.mpCells.value());
    } else {
        suffix = "_nodes";
    }

    const std::string pvd_data_set_name = rModelPartData.mpModelPart->FullName() + suffix;
    const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();

    // append with the step value and rank and extension
    output_vtu_file_name << suffix << "_" << mrModelPart.GetProcessInfo()[STEP]
                         << (r_data_communicator.IsDistributed()
                                 ? "_" + std::to_string(r_data_communicator.Rank())
                                 : "")
                         << ".vtu";

    // write the vtu file.
    std::ofstream output_file;
    output_file.open(output_vtu_file_name.str(), std::ios::out | std::ios::trunc);
    vtk_file_element.Write(output_file);
    output_file.close();

    // if it is run on a distributed system, create the pvtu file.
    if (r_data_communicator.IsDistributed()) {
        return std::make_pair(pvd_data_set_name,
                              WritePartitionedUnstructuredGridData(
                                  *point_data_element, *cell_data_element,
                                  output_vtu_file_name.str(), r_data_communicator));
    }

    // return the final file name for
    return std::make_pair(pvd_data_set_name, output_vtu_file_name.str());
}

template<class TXmlDataElementWrapper>
std::pair<std::string, std::string> VtuOutput::WriteIntegrationPointData(
    const std::string& rOutputFileNamePrefix,
    ModelPartData& rModelPartData,
    TXmlDataElementWrapper& rXmlDataElementWrapper) const
{
    if (!rModelPartData.mpCells.has_value()) {
        // nothing to do here.
        return std::make_pair("", "");
    }

    const auto& integration_point_vars = GetContainerMap(
        mIntegrationPointVariables, rModelPartData.mpCells.value());

    if (integration_point_vars.empty()) {
        // nothing to do here.
        return std::make_pair("", "");
    }

    // create the vtk file
    XmlElementsArray vtk_file_element("VTKFile");
    vtk_file_element.AddAttribute("type", "UnstructuredGrid");
    vtk_file_element.AddAttribute("version", "0.1");
    vtk_file_element.AddAttribute("byte_order", GetEndianness());

    // create the unstructured grid
    auto unstructured_grid_element = Kratos::make_shared<XmlElementsArray>("UnstructuredGrid");
    vtk_file_element.AddElement(unstructured_grid_element);

    DenseVector<std::size_t> offsets;

    const auto total_gauss_points = std::visit([&offsets](auto p_container) {
        // resize the offsets
        offsets.resize(p_container->size(), false);

        // first entity gp offset is zero.
        offsets[0] = 0;

        // now compute the offsets for each entity. This allows
        // having different number of gps in different entities.
        // which is the case if we have a model part with mixed type
        // of elements.
        for (IndexType i = 1; i < p_container->size(); ++i) {
            const auto& r_entity = *(p_container->begin() + i);
            const auto number_of_gps = r_entity.GetGeometry().IntegrationPointsNumber(r_entity.GetIntegrationMethod());
            offsets[i] = offsets[i - 1] + number_of_gps;
        }

        // now add the number of gps of the first entity.
        const auto& r_entity = *(p_container->begin());
        return offsets[offsets.size() - 1] + r_entity.GetGeometry().IntegrationPointsNumber(r_entity.GetIntegrationMethod());
    }, rModelPartData.mpCells.value());

    // create the piece element
    auto piece_element = Kratos::make_shared<XmlElementsArray>("Piece");
    piece_element->AddAttribute("NumberOfPoints", std::to_string(total_gauss_points));
    piece_element->AddAttribute("NumberOfCells", "0");
    unstructured_grid_element->AddElement(piece_element);

    // construct the gauss point position data
    DenseVector<unsigned int> gauss_point_nd_data_shape(2);
    gauss_point_nd_data_shape[0] = total_gauss_points;
    gauss_point_nd_data_shape[1] = 3;
    auto gauss_point_positions = Kratos::make_shared<NDData<double>>(gauss_point_nd_data_shape);
    const auto span = gauss_point_positions->ViewData();

    std::visit([&span, &offsets](auto p_container){
        IndexPartition<IndexType>(p_container->size()).for_each(array_1d<double, 3>{}, [&span, &p_container, &offsets](const auto Index, auto& rTLS) {
            const auto& r_entity = *(p_container->begin() + Index);
            const auto number_of_gauss_points = r_entity.GetGeometry().IntegrationPointsNumber(r_entity.GetIntegrationMethod());
            for (IndexType i = 0; i < number_of_gauss_points; ++i) {
                r_entity.GetGeometry().GlobalCoordinates(rTLS, i);
                std::copy(rTLS.begin(), rTLS.end(), span.begin() + offsets[Index] * 3 + i * 3);
            }

        });
    }, rModelPartData.mpCells.value());

    // create the gauss points element
    auto points_element = Kratos::make_shared<XmlElementsArray>("Points");
    points_element->AddElement(rXmlDataElementWrapper.Get("Position", gauss_point_positions));
    piece_element->AddElement(points_element);

    auto cells_element = Kratos::make_shared<XmlElementsArray>("Cells");
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
    piece_element->AddElement(point_data_element);

    // add the gauss point data
    bool is_gauss_point_data_available = false;
    for (const auto& r_pair : integration_point_vars) {
        std::visit([&offsets, &point_data_element, &rXmlDataElementWrapper, &rModelPartData, &is_gauss_point_data_available, total_gauss_points](auto pVariable, auto pContainer) {
            // type information of the variable
            using data_type = typename BareType<decltype(*pVariable)>::Type;
            using data_type_traits = DataTypeTraits<data_type>;
            using primitive_data_type = typename data_type_traits::PrimitiveType;

            if constexpr(std::is_same_v<primitive_data_type, double>) {
                // here we cannot use GaussPointVariableTensorAdaptor because
                // it does not allow recording gauss point information if different
                // entities have different number of gauss points.
                // hence we do the computation here manually.

                std::vector<data_type> output;

                // first we need to find out the shape of the gauss point data
                std::vector<unsigned int> local_shape(data_type_traits::Dimension, 0u);
                if (!pContainer->empty()) {
                    pContainer->front().CalculateOnIntegrationPoints(
                        *pVariable, output, rModelPartData.mpModelPart->GetProcessInfo());

                    if (!output.empty()) {
                        // if there are available gauss point information
                        data_type_traits::Shape(output.front(), local_shape.data(),
                                                local_shape.data() + local_shape.size());
                    }
                }

                // now do the communication between ranks to get the correct size
                const auto& max_local_shape = rModelPartData.mpModelPart->GetCommunicator().GetDataCommunicator().MaxAll(local_shape);

                // now we construct the nd_data_shape
                DenseVector<unsigned int> nd_data_shape(local_shape.size() + 1);
                std::copy(max_local_shape.begin(), max_local_shape.end(), nd_data_shape.begin() + 1);
                nd_data_shape[0] = total_gauss_points;

                const auto total_number_of_components = std::accumulate(max_local_shape.begin(), max_local_shape.end(), 1u, std::multiplies<unsigned int>{});

                auto p_gauss_data = Kratos::make_shared<NDData<double>>(nd_data_shape);
                auto span = p_gauss_data->ViewData();

                if (span.size() > 0) {
                    is_gauss_point_data_available = true;
                    IndexPartition<IndexType>(pContainer->size()).for_each(output, [&span, &pContainer, &pVariable, &rModelPartData, &offsets, total_number_of_components](const auto Index, auto& rTLS) {
                        auto& r_entity = *(pContainer->begin() + Index);
                        r_entity.CalculateOnIntegrationPoints(*pVariable, rTLS, rModelPartData.mpModelPart->GetProcessInfo());
                        DataTypeTraits<std::vector<data_type>>::CopyToContiguousData(span.begin() + offsets[Index] * total_number_of_components, rTLS);
                    });
                    point_data_element->AddElement(rXmlDataElementWrapper.Get(pVariable->Name(), p_gauss_data));
                }
            }
        }, r_pair.second, rModelPartData.mpCells.value());
    }

    const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();

    if (r_data_communicator.OrReduceAll(is_gauss_point_data_available)) {
        std::stringstream output_vtu_file_name;
        output_vtu_file_name
            << rOutputFileNamePrefix << "/" << rModelPartData.mpModelPart->FullName()
            << "_gauss_" << mrModelPart.GetProcessInfo()[STEP]
            << (r_data_communicator.IsDistributed()
                    ? "_" + std::to_string(r_data_communicator.Rank())
                    : "")
            << ".vtu";

        std::ofstream output_file;
        output_file.open(output_vtu_file_name.str(), std::ios::out | std::ios::trunc);
        vtk_file_element.Write(output_file);
        output_file.close();

        if (r_data_communicator.IsDistributed()) {
            return std::make_pair(
                rModelPartData.mpModelPart->FullName() + "_gauss",
                WritePartitionedUnstructuredGridData(
                    *point_data_element, *Kratos::make_shared<XmlElementsArray>(""),
                    output_vtu_file_name.str(), r_data_communicator));
        }

        return std::make_pair(rModelPartData.mpModelPart->FullName() + "_gauss",
                              output_vtu_file_name.str());
    }
    else {
        return std::make_pair("", "");
    }
}

void VtuOutput::PrintOutput(const std::string& rOutputFileNamePrefix)
{
    KRATOS_TRY

    const auto& r_process_info = mrModelPart.GetProcessInfo();

    // check if the step and time variables are there. These
    // are must for the *.pvd file and step file writing.
    KRATOS_ERROR_IF_NOT(r_process_info.Has(STEP))
        << "STEP variable is not found in the process info of " << mrModelPart.FullName() << ".\n";

    KRATOS_ERROR_IF_NOT(r_process_info.Has(TIME))
        << "TIME variable is not found in the process info of " << mrModelPart.FullName() << ".\n";

    // Add the step info.
    mStepInfo.push_back(std::make_pair(r_process_info[STEP], r_process_info[TIME]));

    std::filesystem::create_directories(rOutputFileNamePrefix);

    std::vector<std::pair<std::string, std::string>> pvd_file_name_info;

    for (auto& r_model_part_data : mListOfModelPartData) {
        switch (mOutputFormat) {
            case ASCII:
            {
                XmlAsciiNDDataElementWrapper data_element_wrapper{mPrecision, &mrModelPart.GetCommunicator().GetDataCommunicator()};
                // first write the unstructured grid data
                pvd_file_name_info.push_back(WriteUnstructuredGridData(
                    rOutputFileNamePrefix, r_model_part_data, data_element_wrapper));

                // now write the gauss point info
                pvd_file_name_info.push_back(WriteIntegrationPointData(
                    rOutputFileNamePrefix, r_model_part_data, data_element_wrapper));
                break;
            }
            case BINARY:
            {
                XmlBase64BinaryNDDataElementWrapper data_element_wrapper{&mrModelPart.GetCommunicator().GetDataCommunicator()};
                // first write the unstructured grid data
                pvd_file_name_info.push_back(WriteUnstructuredGridData(
                    rOutputFileNamePrefix, r_model_part_data, data_element_wrapper));

                // now write the gauss point info
                pvd_file_name_info.push_back(WriteIntegrationPointData(
                    rOutputFileNamePrefix, r_model_part_data, data_element_wrapper));
                break;
            }
        }
    }

    // now generate the *.pvd file
    if (mrModelPart.GetCommunicator().MyPID() == 0) {
        // Single pvd file links all the vtu files from sum-model parts
        // partitioned model_parts and time step vtu files together.

        // creates the pvd file element
        XmlElementsArray pvd_file_element("VTKFile");
        pvd_file_element.AddAttribute("type", "Collection");
        pvd_file_element.AddAttribute("version", "1.0");
        pvd_file_element.AddAttribute("byte_order", GetEndianness());

        // creates the collection element
        auto collection_element = Kratos::make_shared<XmlElementsArray>("Collection");
        pvd_file_element.AddElement(collection_element);

        // now iterate through all the time steps and correctly write
        // the file names for each time step.
        for (auto& [current_step, current_time] : mStepInfo) {
            IndexType local_index = 0;
            for (IndexType i = 0; i < pvd_file_name_info.size(); ++i) {
                if (pvd_file_name_info[i].second != "") {
                    auto current_element = Kratos::make_shared<XmlElementsArray>("DataSet");

                    // write the time with the specified precision.
                    std::stringstream str_time;
                    str_time << std::scientific << std::setprecision(mPrecision) << current_time;

                    // linking file name. This file name is the vtu file name of
                    // the current step. We need to adjust the filename accordingly
                    // to get the file name corresponding to current_step. When we
                    // were writing the file names, the step was appended. So here we
                    // remove the appended step.
                    const auto vtu_file_name = pvd_file_name_info[i].second;
                    auto current_vtu_file_name = vtu_file_name.substr(0, vtu_file_name.rfind("_"));
                    current_vtu_file_name += "_" + std::to_string(current_step);
                    // now we need to add back the correct extension
                    current_vtu_file_name += vtu_file_name.substr(vtu_file_name.rfind("."));

                    current_element->AddAttribute("timestep", str_time.str());
                    current_element->AddAttribute("name", pvd_file_name_info[i].first);
                    current_element->AddAttribute("part", std::to_string(local_index++));
                    current_element->AddAttribute(
                        "file", std::filesystem::relative(
                                    std::filesystem::absolute(current_vtu_file_name),
                                    std::filesystem::absolute(rOutputFileNamePrefix).parent_path())
                                    .string());
                    collection_element->AddElement(current_element);
                }
            }
        }

        std::ofstream output_file;
        output_file.open(rOutputFileNamePrefix + ".pvd", std::ios::out | std::ios::trunc);
        pvd_file_element.Write(output_file);
        output_file.close();
    }

    KRATOS_CATCH("");
}

std::string VtuOutput::Info() const
{
    std::stringstream info;
    info << "VtuOutput: " << mrModelPart.FullName() << " [ writer = ";
    switch (mOutputFormat) {
        case ASCII:
            info << "ASCII";
            break;
        case BINARY:
            info << "BINARY";
            break;
    }
    info << ", precision = " << mPrecision << ", is initial configuration = "
         << (mIsInitialConfiguration ? "yes" : "no") << " ]";
    return info.str();
}

void VtuOutput::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void VtuOutput::PrintData(std::ostream& rOStream) const
{
    PrintDataLocationData(rOStream, "flag", mFlags);
    rOStream << "\n";
    PrintDataLocationData(rOStream, "variable", mVariables);
    rOStream << "\n";
    PrintDataLocationData(rOStream, "integration variable", mIntegrationPointVariables);
    rOStream << "\n";

    rOStream << "List of model part info:";
    for (const auto& r_model_part_data : mListOfModelPartData) {
        rOStream << "\n\tModel part: \"" << r_model_part_data.mpModelPart->FullName() << "\"";
        if (r_model_part_data.mpCells.has_value()) {
            std::visit([&rOStream](auto pContainer){
                using container_type = BareType<decltype(*pContainer)>;
                rOStream << " with " << pContainer->size() << " " << ModelPart::Container<container_type>::GetEntityName() << "s";
            }, r_model_part_data.mpCells.value());
        } else {
            rOStream << " with " << r_model_part_data.mpModelPart->NumberOfNodes() << " nodes";
        }
        rOStream << ", used for point fields = " << (r_model_part_data.UsePointsForDataFieldOutput ? "yes" : "no");

        if (r_model_part_data.UsePointsForDataFieldOutput) {
            rOStream << "\n\t\t" << "Point fields:";
            for (const auto& [name, field] : r_model_part_data.mPointFields) {
                std::visit([&rOStream, &name](auto pNDData){
                    rOStream << "\n\t\t\t" << name << ": " << *pNDData;
                }, field);
            }
        }

        rOStream << "\n\t\t" << "Cell fields:";
        for (const auto& [name, field] : r_model_part_data.mCellFields) {
            std::visit([&rOStream, &name](auto pNDData){
                rOStream << "\n\t\t\t" << name << ": " << *pNDData;
            }, field);
        }
    }
}

} // namespace Kratos
