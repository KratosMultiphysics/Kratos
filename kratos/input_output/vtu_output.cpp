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

/// ----------------------------------------- NDDataElementWrappers ----------------------------------------- ///
struct XmlAsciiNDDataElementWrapper
{
    IndexType mPrecision;

    template<class NDDataTypePointer>
    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataTypePointer pNDData)
    {
        using data_type = typename BareType<decltype(*pNDData)>::DataType;
        return Kratos::make_shared<XmlAsciiNDDataElement<data_type>>(
            rDataArrayName, pNDData, pNDData->Shape(), mPrecision);
    }
};

struct XmlBase64BinaryNDDataElementWrapper
{
    template<class NDDataTypePointer>
    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataTypePointer pNDData)
    {
        using data_type = typename BareType<decltype(*pNDData)>::DataType;
        return Kratos::make_shared<XmlBase64BinaryNDDataElement<data_type>>(
            rDataArrayName, pNDData, pNDData->Shape());
    }
};

/// ----------------------------------------- NDDataElementWrappers ----------------------------------------- ///

std::string GetEndianness()
{
    int i = 0x0001;

    if (*reinterpret_cast<char*>(&i) != 0) {
        return "LittleEndian";
    } else {
        return "BigEndian";
    }
}

template<class... T>
void CheckDataArrayName(
    const std::string& rName,
    const std::vector<Globals::DataLocation>& rLocations,
    const VtuOutput::DataMap<T>&... rMaps)
{
    for (const auto& r_location : rLocations) {
        const bool found_existing_name =
            !(... && (rMaps.find(r_location) == rMaps.end() ||
                      rMaps.find(r_location)->second.find(rName) ==
                          rMaps.find(r_location)->second.end()));
        KRATOS_ERROR_IF(found_existing_name)
                << "Found an existing data array with the same name = \"" << rName << "\".\n";
    }
}

void CheckDataArrayName(
    const std::string& rName,
    const Globals::DataLocation& rLocation,
    const VtuOutput::UnstructuredGridData& rUnstructuredGridData)
{
    bool found_existing_name = false;
    switch (rLocation) {
        case Globals::DataLocation::NodeHistorical:
        case Globals::DataLocation::NodeNonHistorical:
            found_existing_name = rUnstructuredGridData.mPointFields.find(rName) != rUnstructuredGridData.mPointFields.end();
            break;
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            found_existing_name = rUnstructuredGridData.mCellFields.find(rName) != rUnstructuredGridData.mCellFields.end();
            break;
        default:
            KRATOS_ERROR << "Unsupported data location type.";
    }
    KRATOS_ERROR_IF(found_existing_name)
            << "Found an existing data array with the same name = \"" << rName << "\".\n";
}

void CheckDataArrayName(
    const std::string& rName,
    const Globals::DataLocation& rLocation,
    const std::vector<VtuOutput::UnstructuredGridData>& rListOfUnstructuredGridData)
{
    for (const auto& r_model_part_data : rListOfUnstructuredGridData) {
        switch (rLocation) {
            case Globals::DataLocation::NodeHistorical:
            case Globals::DataLocation::NodeNonHistorical:
                if (r_model_part_data.UsePointsForDataFieldOutput) {
                    CheckDataArrayName(rName, rLocation, r_model_part_data);
                }
                break;
            default:
                CheckDataArrayName(rName, rLocation, r_model_part_data);
        }
    }
}

std::string GetEntityName(const std::optional<VtuOutput::CellContainerPointerType>& pCellContainer)
{
    if (pCellContainer.has_value()) {
        return std::visit([](auto p_cell_container) {
            using container_type = BareType<decltype(*p_cell_container)>;
            return ModelPart::Container<container_type>::GetEntityName();
        }, pCellContainer.value());
    } else {
        return "node";
    }
}

template<class NDDataPointerType>
NDDataPointerType GetNonIgnoredNDData(
    const std::set<IndexType>& rIgnoredIndices,
    NDDataPointerType pNDData)
{
    if (rIgnoredIndices.empty()) {
        return pNDData;
    } else {
        const auto& shape = pNDData->Shape();
        const auto number_of_components = std::accumulate(shape.begin() + 1, shape.end(), 1u, std::multiplies<unsigned int>{});

        using data_type = typename BareType<decltype(*pNDData)>::DataType;

        DenseVector<unsigned int> nd_shape(pNDData->Shape());
        nd_shape[0] = shape[0] - rIgnoredIndices.size();
        auto p_destination_nd_data = Kratos::make_shared<NDData<data_type>>(nd_shape);

        auto origin_span = pNDData->ViewData();
        auto destination_span = p_destination_nd_data->ViewData();

        IndexType origin_offset{}, destination_offset{};
        for (auto itr = rIgnoredIndices.begin(); itr != rIgnoredIndices.end(); ++itr) {
            std::copy(origin_span.begin() + origin_offset, origin_span.begin() + (*itr) * number_of_components, destination_span.begin() + destination_offset);
            destination_offset += (*itr) * number_of_components - origin_offset;
            origin_offset = (*itr + 1) * number_of_components;
        }
        std::copy(origin_span.begin() + origin_offset, origin_span.end(), destination_span.begin() + destination_offset);

        return p_destination_nd_data;
    }
}

template<class TContainerType>
NDData<unsigned char>::Pointer GetGeometryTypes(
    std::set<IndexType>& rIgnoredIndices,
    const TContainerType& rContainer,
    const IndexType EchoLevel)
{
    auto p_geometry_types = Kratos::make_shared<NDData<unsigned char>>(DenseVector<unsigned int>(1, rContainer.size()));
    auto span = p_geometry_types->ViewData();

    DenseVector<IndexType> ignored_indices(rContainer.size(), 0);

    IndexPartition<IndexType>(rContainer.size()).for_each([&span, &rContainer, &ignored_indices, EchoLevel](const IndexType Index) {

        const auto p_itr = VtkDefinitions::KratosVtkGeometryTypes.find((rContainer.begin() + Index)->GetGeometry().GetGeometryType());
        if (p_itr != VtkDefinitions::KratosVtkGeometryTypes.end()) {
            *(span.begin() + Index) = static_cast<unsigned char>(p_itr->second);
        } else {
            ignored_indices[Index] = 1;
            KRATOS_WARNING_IF("VtuOutput", EchoLevel > 1)
                << "Skipping unsupported geometry type in "
                << ModelPart::Container<TContainerType>::GetEntityName()
                << " with id " << (rContainer.begin() + Index)->Id() << ".\n";
        }
    });

    for (IndexType i = 0; i < rContainer.size(); ++i) {
        if (ignored_indices[i] == 1) {
            rIgnoredIndices.insert(i);
        }
    }
    return p_geometry_types;
}

template <class TContainerType>
NDData<int>::Pointer GetOffsets(
    const std::set<IndexType>& rIgnoredIndices,
    const TContainerType& rContainer)
{
    auto p_offsets = Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, rContainer.size() - rIgnoredIndices.size()));
    auto span = p_offsets->ViewData();

    int total_offset = 0;
    auto data_itr = span.begin();
    IndexType ignored_indices_offset{};
    for (auto itr = rIgnoredIndices.begin(); itr != rIgnoredIndices.end(); ++itr) {
        for (IndexType i = ignored_indices_offset; i < *itr; ++i) {
            total_offset += (rContainer.begin() + i)->GetGeometry().size();
            *(data_itr++) = total_offset;
            ignored_indices_offset = (*itr + 1);
        }
    }

    for (IndexType i = ignored_indices_offset; i < rContainer.size(); ++i) {
            total_offset += (rContainer.begin() + i)->GetGeometry().size();
            *(data_itr++) = total_offset;
    }

    return p_offsets;
}

template<class TContainerType>
NDData<int>::Pointer GetConnectivities(
    const NDData<int>& rOffsets,
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap,
    const std::set<IndexType>& rIgnoredIndices)
{
    const auto offsets_span = rOffsets.ViewData();
    auto p_connectivities = Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, offsets_span.back()));
    auto connectivities_span = p_connectivities->ViewData();

    IndexPartition<IndexType>(rContainer.size()).for_each([&connectivities_span, &offsets_span, &rContainer, &rKratosVtuIndicesMap, &rIgnoredIndices](const IndexType Index) {
        if (rIgnoredIndices.find(Index) == rIgnoredIndices.end()) {
            IndexType offsets_index = Index - std::distance(rIgnoredIndices.begin(), rIgnoredIndices.upper_bound(Index));
            const auto& r_geometry = (rContainer.begin() + Index)->GetGeometry();
            auto entity_data_begin_itr = connectivities_span.begin() + offsets_span[offsets_index] - r_geometry.size();

            for (const auto& r_node : r_geometry) {
                const auto p_itr = rKratosVtuIndicesMap.find(r_node.Id());
                if (p_itr != rKratosVtuIndicesMap.end()) {
                    *(entity_data_begin_itr++) = p_itr->second;
                } else {
                    KRATOS_ERROR << "Node with id " << r_node.Id() << " not found in nodes list.";
                }
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
void AddConnectivityData(
    XmlElementsArray& rCellElement,
    std::set<IndexType>& rIgnoredIndices,
    const ModelPart::NodesContainerType& rNodes,
    VtuOutput::CellContainerPointerType pCells,
    TXmlDataElementWrapper& rXmlDataElementWrapper,
    const IndexType EchoLevel)
{
    std::visit([&rCellElement, &rXmlDataElementWrapper, &rIgnoredIndices, &rNodes, EchoLevel](auto p_container) {
        std::unordered_map<IndexType, IndexType> indices_map;
        indices_map.reserve(rNodes.size());
        IndexType vtu_index = 0;
        for (const auto& r_node : rNodes) {
            indices_map[r_node.Id()] = vtu_index++;
        }

        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2) << "------ Collecting Vtk geometry type info...\n";
        rCellElement.AddElement(rXmlDataElementWrapper.Get("types", GetGeometryTypes(rIgnoredIndices, *p_container, EchoLevel)));

        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2)
            << "------ Ignored " << rIgnoredIndices.size() << "/"
            << p_container->size() << " " << GetEntityName(p_container) << "(s).\n";

        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2) << "------ Collecting geometry offsets info...\n";
        auto p_offsets = GetOffsets(rIgnoredIndices, *p_container);

        rCellElement.AddElement(rXmlDataElementWrapper.Get("offsets", p_offsets));

        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2) << "------ Collecting geometry connectivity info...\n";
        rCellElement.AddElement(rXmlDataElementWrapper.Get("connectivity", GetConnectivities(*p_offsets, *p_container, indices_map, rIgnoredIndices)));
    }, pCells);
}

template<class TTensorAdaptorType, class TContainerPointerType, class TMapType, class TXmlDataElementWrapper, class... TArgs>
void AddFieldsFromTensorAdaptor(
    XmlElementsArray& rXmlElement,
    TContainerPointerType pContainer,
    const TMapType& rMap,
    TXmlDataElementWrapper& rXmlDataElementWrapper,
    const DataCommunicator& rDataCommunicator,
    const std::set<IndexType>& rIgnoredIndices,
    TArgs&&... rArgs)
{
    KRATOS_TRY

    for (const auto& r_pair : rMap) {
        using data_type = BareType<decltype(r_pair.second)>;
        if constexpr(std::is_same_v<data_type, Flags>) {
            // the map is of type flags
            // here we don't need to do any communication because Flags are always having a static data shape.
            TTensorAdaptorType tensor_adaptor(pContainer, *r_pair.second, rArgs...);
            tensor_adaptor.CollectData();
            rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, GetNonIgnoredNDData(rIgnoredIndices, tensor_adaptor.pGetStorage())));
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
                    } else if constexpr(std::is_same_v<TTensorAdaptorType, VariableTensorAdaptor>) {
                        ContainerIOUtils::CopyToContiguousArray<int>(
                            *pContainer, p_nd_data->ViewData(),
                            shape.data().begin(), shape.data().begin() + 1,
                            [p_variable](int& rValue, const auto& rNode) {
                                rValue = rNode.GetValue(*p_variable);
                            });
                    } else {
                        KRATOS_ERROR << "Unsupported tensor adaptor type.";
                    }

                    // since we only support Variable<int>, which is having a static data shape
                    // we don't have to do mpi communication to decide the shape on the
                    // empty ranks.
                    rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, GetNonIgnoredNDData(rIgnoredIndices, p_nd_data)));
                } else if constexpr(std::is_same_v<primitive_data_type, double>) {
                    using container_type = BareType<decltype(*pContainer)>;
                    using entity_type = typename container_type::data_type;
                    using data_type = typename BareType<decltype(*p_variable)>::Type;
                    using data_type_traits = DataTypeTraits<data_type>;

                    if constexpr(data_type_traits::IsDynamic) {
                        // this is a dynamic type such as Vector or Matrix, so
                        // we need to do communication to decide the correct size.

                        // construct the correct data_shape
                        std::vector<unsigned int> data_shape(data_type_traits::Dimension, 0);
                        data_type value;
                        if (!pContainer->empty()) {
                            if constexpr(std::is_same_v<TTensorAdaptorType, HistoricalVariableTensorAdaptor>) {
                                value = static_cast<const entity_type&>(pContainer->front()).FastGetSolutionStepValue(*p_variable);
                            } else if constexpr(std::is_same_v<TTensorAdaptorType, VariableTensorAdaptor>) {
                                value = static_cast<const entity_type&>(pContainer->front()).GetValue(*p_variable);
                            } else {
                                KRATOS_ERROR << "Unsupported tensor adaptor type.";
                            }
                            data_type_traits::Shape(value, data_shape.data(), data_shape.data() + data_type_traits::Dimension);
                        }
                        const auto max_data_shape = rDataCommunicator.MaxAll(data_shape);
                        DenseVector<unsigned int> nd_shape(max_data_shape.size() + 1);
                        std::copy(max_data_shape.begin(), max_data_shape.end(), nd_shape.begin() + 1);
                        nd_shape[0] = pContainer->size();

                        // construct the data storage
                        auto p_nd_data = Kratos::make_shared<typename TTensorAdaptorType::Storage>(nd_shape);
                        auto base_ta = TensorAdaptor<double>(pContainer, p_nd_data, false);
                        TTensorAdaptorType tensor_adaptor(base_ta, p_variable, rArgs..., false);
                        tensor_adaptor.CollectData();

                        rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, GetNonIgnoredNDData(rIgnoredIndices, tensor_adaptor.pGetStorage())));
                    } else {
                        // this is a static type such as double, array_1d<double, 3>, ...
                        // So no need of mpi communication
                        TTensorAdaptorType tensor_adaptor(pContainer, p_variable, rArgs...);
                        tensor_adaptor.CollectData();
                        rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, GetNonIgnoredNDData(rIgnoredIndices, tensor_adaptor.pGetStorage())));
                    }
                } else {
                    KRATOS_ERROR << "Unsupported variable type.";
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
    TXmlDataElementWrapper& rXmlDataElementWrapper,
    const std::set<IndexType>& rIgnoredIndices)
{
    for (const auto& r_pair : rMap) {
        std::visit([&rXmlElement, &r_pair, &rXmlDataElementWrapper, &rIgnoredIndices](auto pNDData) {
            rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, GetNonIgnoredNDData(rIgnoredIndices, pNDData)));
        }, r_pair.second);
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

void AddUnstructuredGridData(
    std::vector<VtuOutput::UnstructuredGridData>& rOutput,
    ModelPart& rModelPart,
    const IndexType EchoLevel)
{
    const bool has_elements   = rModelPart.GetCommunicator().GlobalNumberOfElements() > 0;
    const bool has_conditions = rModelPart.GetCommunicator().GlobalNumberOfConditions() > 0;

    if (has_elements) {
        // Model part has elements. Hence add a separate output
        // for elements.

        // now check if it has proper nodes
        if (rModelPart.NumberOfNodes() > 0) {
            KRATOS_INFO_IF("VtuOutput", EchoLevel > 0)
                << "Configuring output for \"" << rModelPart.FullName() << "\" elements with existing nodes container.\n";
            VtuOutput::UnstructuredGridData model_part_data{true, &rModelPart, rModelPart.pNodes(), rModelPart.pElements()};
            rOutput.push_back(model_part_data);
        } else {
            // create the nodes container.
            KRATOS_INFO_IF("VtuOutput", EchoLevel > 0)
                << "Configuring output for \"" << rModelPart.FullName() << "\" elements with new nodes container.\n";
            auto p_nodes = GetNodesContainer(rModelPart.Elements());
            VtuOutput::UnstructuredGridData model_part_data{false, &rModelPart, p_nodes, rModelPart.pElements()};
            rOutput.push_back(model_part_data);
        }
    }

    if (has_conditions) {
        // Model part has conditions. Hence add a separate output
        // for conditions.

        if (!has_elements && rModelPart.NumberOfNodes() > 0) {
            KRATOS_INFO_IF("VtuOutput", EchoLevel > 0)
                << "Configuring output for \"" << rModelPart.FullName() << "\" conditions with existing nodes container.\n";
            VtuOutput::UnstructuredGridData model_part_data{true, &rModelPart, rModelPart.pNodes(), rModelPart.pConditions()};
            rOutput.push_back(model_part_data);
        } else {
            // either this model part also contains elements, or it does not
            // contain nodes. In either case, the nodes list given by the rModelPart
            // does not reflect the actual nodes used by the conditions.
            // In order to avoid writing nodes, which are not used by the conditions,
            // this will use a new nodes list.

            KRATOS_INFO_IF("VtuOutput", EchoLevel > 0)
                << "Configuring output for \"" << rModelPart.FullName() << "\" conditions with new nodes container.\n";

            auto p_nodes = GetNodesContainer(rModelPart.Conditions());
            VtuOutput::UnstructuredGridData model_part_data{false, &rModelPart, p_nodes, rModelPart.pConditions()};
            rOutput.push_back(model_part_data);
        }
    }

    if (!has_elements && !has_conditions) {
        KRATOS_INFO_IF("VtuOutput", EchoLevel > 0)
            << "Configuring output for \"" << rModelPart.FullName() << "\" nodes.\n";
        // Model part does not have either conditions or elements.
        // Hence, only adding the nodes.
        VtuOutput::UnstructuredGridData model_part_data{true, &rModelPart, rModelPart.pNodes(), std::nullopt};
        rOutput.push_back(model_part_data);
    }
}

void AddUnstructuredGridDataRecursively(
    std::vector<VtuOutput::UnstructuredGridData>& rOutput,
    ModelPart& rModelPart,
    const IndexType EchoLevel)
{
    // add the current model part data to the output.
    AddUnstructuredGridData(rOutput, rModelPart, EchoLevel);

    // now recursively add all the sub model part data.
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        AddUnstructuredGridDataRecursively(rOutput, r_sub_model_part, EchoLevel);
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
    const bool OutputSubModelParts,
    const IndexType EchoLevel)
    : mrModelPart(rModelPart),
      mIsInitialConfiguration(IsInitialConfiguration),
      mEchoLevel(EchoLevel),
      mOutputFormat(OutputFormat),
      mPrecision(Precision)
{
    if (OutputSubModelParts) {
        // Collect all model part data recursively.
        AddUnstructuredGridDataRecursively(mListOfUnstructuredGridData, rModelPart, mEchoLevel);
    } else {
        // Only collect the data from the Passed model part
        AddUnstructuredGridData(mListOfUnstructuredGridData, rModelPart, mEchoLevel);
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
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mStepInfo.empty())
        << "Flags can be added only before the first call to the PrintOutput [ flag name = "
        << rFlagName << " ].\n";

    switch (DataLocation) {
        case Globals::DataLocation::NodeNonHistorical:
            // additionally check here in the historical containers.
            CheckDataArrayName(rFlagName, {Globals::DataLocation::NodeHistorical}, mFlags, mVariables);
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            CheckDataArrayName(rFlagName, {DataLocation}, mFlags, mVariables);
            CheckDataArrayName(rFlagName, DataLocation, mListOfUnstructuredGridData);
            mFlags[DataLocation][rFlagName] = &rFlagVariable;
            break;
        default:
            KRATOS_ERROR << "Flags can be only added to NodeNonHistorical, Condition, and Element data locations.";
            break;
    }

    KRATOS_CATCH("");
}

void VtuOutput::AddVariable(
    SupportedVariablePointerType pVariable,
    Globals::DataLocation DataLocation)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mStepInfo.empty())
        << "Variables can be added only before the first call to the PrintOutput [ variable name = "
        << GetName(pVariable) << " ].\n";

    switch (DataLocation) {
        case Globals::DataLocation::NodeHistorical:
            // additionally check here in the non historical containers.
            CheckDataArrayName(GetName(pVariable), {Globals::DataLocation::NodeNonHistorical}, mFlags, mVariables);
        case Globals::DataLocation::NodeNonHistorical:
            // additionally check here in the historical containers.
            CheckDataArrayName(GetName(pVariable), {Globals::DataLocation::NodeHistorical}, mFlags, mVariables);
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            CheckDataArrayName(GetName(pVariable), {DataLocation}, mFlags, mVariables);
            CheckDataArrayName(GetName(pVariable), DataLocation, mListOfUnstructuredGridData);
            mVariables[DataLocation][GetName(pVariable)] = pVariable;
            break;
        default:
            KRATOS_ERROR << "Variables can be only added to NodeHistorical, NodeNonHistorical, Condition, and Element data locations.";
            break;
    }

    KRATOS_CATCH("");
}

void VtuOutput::AddIntegrationPointVariable(
    SupportedVariablePointerType pVariable,
    Globals::DataLocation DataLocation)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mStepInfo.empty())
        << "Integration point variables can be added only before the first call to the PrintOutput [ integration point variable name = "
        << GetName(pVariable) << " ].\n";

    switch (DataLocation) {
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            CheckDataArrayName(GetName(pVariable), {DataLocation}, mIntegrationPointVariables);
            mIntegrationPointVariables[DataLocation][GetName(pVariable)] = pVariable;
            break;
        default:
            KRATOS_ERROR << "Integration point variables can be only added to Condition, and Element data locations.";
            break;
    }

    KRATOS_CATCH("");
}

void VtuOutput::AddContainerExpression(
    const std::string& rExpressionName,
    SupportedContainerExpressionPointerType pContainerExpression)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mStepInfo.empty())
        << "Expressions can be added only before the first call to the PrintOutput [ expression name = "
        << rExpressionName << " ].\n";

    std::visit([this, &rExpressionName](auto p_container_expression) {
        const auto& r_expression = p_container_expression->GetExpression();

        const auto number_of_data_components = r_expression.GetItemComponentCount();
        const auto& data_shape = r_expression.GetItemShape();

        auto& r_expression_container = p_container_expression->GetContainer();

        using expression_container_type = BareType<decltype(r_expression_container)>;

        for (IndexType i_model_part_data = 0; i_model_part_data < this->mListOfUnstructuredGridData.size(); ++i_model_part_data) {
            auto& r_model_part_data = this->mListOfUnstructuredGridData[i_model_part_data];

            // get the model part
            auto p_model_part = r_model_part_data.mpModelPart;

            // expressions are created from local mesh entities.
            if constexpr(std::is_same_v<expression_container_type, ModelPart::NodesContainerType>) {
                if (&r_expression_container == &p_model_part->GetCommunicator().LocalMesh().Nodes()) {
                    if (r_model_part_data.UsePointsForDataFieldOutput) {
                        // first check if there are any nodal fields present with the same name already
                        CheckDataArrayName(rExpressionName, {Globals::DataLocation::NodeHistorical, Globals::DataLocation::NodeNonHistorical}, mFlags, mVariables);
                    }

                    // now check if the rExpressionName is found in the mPointFields of this model part
                    CheckDataArrayName(rExpressionName, Globals::DataLocation::NodeNonHistorical, r_model_part_data);

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
                    if (std::visit([this, &r_model_part_data, &rExpressionName, &r_expression, &r_expression_container, &data_shape, number_of_data_components](auto p_model_part_container) {
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

                                KRATOS_ERROR_IF_NOT(r_model_part_data.mCellFields.find(rExpressionName) == r_model_part_data.mCellFields.end())
                                    << "Found an existing data array with the same expression name = \"" << rExpressionName << "\" [ "
                                    << "expression : " << r_expression << " ]\n" << *this;

                                r_model_part_data.mCellFields[rExpressionName] = p_nd_data;
                                return true;
                            }
                        }

                        return false;
                    }, r_model_part_data.mpCells.value())) {
                        return true;
                    }
                }
            }
        }

        KRATOS_ERROR
            << "The container in the ContainerExpression is not referring to any of the containers "
            << "written by this VTU output [ container expression name = " << rExpressionName
            << ", tensor_adaptor = " << *p_container_expression << " ]\n"
            << *this;

    }, pContainerExpression);

    KRATOS_CATCH("");
}

void VtuOutput::AddTensorAdaptor(
    const std::string& rTensorAdaptorName,
    SupportedTensorAdaptorPointerType pTensorAdaptor)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mStepInfo.empty())
        << "TensorAdaptors can be added only before the first call to the PrintOutput [ tensor adaptor name = "
        << rTensorAdaptorName << " ].\n";

    std::visit([this, &rTensorAdaptorName](auto p_tensor_adaptor) {
        using tensor_adaptor_type = BareType<decltype(*p_tensor_adaptor)>;
        auto shape = p_tensor_adaptor->Shape();
        auto ta_span = p_tensor_adaptor->ViewData();

        // the tensor adaptors may have different number of components in dimensions from 2 to N [Eg.
        // tensor adaptors created from dynamic variables, having zero entities. ] Hence doing
        // communication to correctly size the number of components in higher dimensions.
        std::vector<unsigned int> local_data_shape(shape.begin() + 1, shape.end());
        const auto& max_data_shape = this->mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_data_shape);
        if (shape[0] == 0) std::copy(max_data_shape.begin(), max_data_shape.end(), shape.begin() + 1);

        const auto number_of_data_components = std::accumulate(shape.begin() + 1, shape.end(), 1u, std::multiplies<unsigned int>{});

        const bool ta_added = std::visit([this, &rTensorAdaptorName, &p_tensor_adaptor, &shape, &ta_span, number_of_data_components](auto p_ta_container){
            using ta_container_type = BareType<decltype(*p_ta_container)>;

            // first check to which model part output this tenor adaptor belongs to
            for (IndexType i = 0; i < this->mListOfUnstructuredGridData.size(); ++i) {
                auto& r_model_part_data = this->mListOfUnstructuredGridData[i];
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

                    if (&*(p_ta_container) == &p_model_part->Nodes() || &*(p_ta_container) == &p_model_part->GetCommunicator().LocalMesh().Nodes()) {
                        if (r_model_part_data.UsePointsForDataFieldOutput) {
                            // first check if there are any nodal fields present with the same name already
                            CheckDataArrayName(rTensorAdaptorName, {Globals::DataLocation::NodeHistorical, Globals::DataLocation::NodeNonHistorical}, mFlags, mVariables);
                        }

                        // now check if the rTensorAdaptorName is found in the mPointFields of this model part
                        CheckDataArrayName(rTensorAdaptorName, Globals::DataLocation::NodeNonHistorical, r_model_part_data);
                    }

                    if (&*(p_ta_container) == &p_model_part->Nodes()) {
                        // now we have to check for the shapes, since TensorAdaptors may not have consistent shapes in the
                        // dimensions except for the first dimension [Eg. TensorAdaptor created from a dynamic variable, where
                        // one rank is having zero entities.] Hence needs to communicate the dimensions
                        r_model_part_data.mPointFields[rTensorAdaptorName] = Kratos::make_shared<typename tensor_adaptor_type::Storage>(ta_span.data(), shape);;
                        return true;
                    } else if (&*(p_ta_container) == &p_model_part->GetCommunicator().LocalMesh().Nodes()) {
                        // the tensor adaptor is having a container referring to the local mesh nodes.
                        // now we have to construct the nodal values for all the nodes including ghost mesh
                        // nodes.

                        // first clear the TENSOR_ADAPTOR_SYNC.
                        VariableUtils().SetNonHistoricalVariableToZero(TENSOR_ADAPTOR_SYNC, p_model_part->Nodes());

                        // secondly fill in the local nodal values to the temporary variable TENSOR_ADAPTOR_SYNC
                        IndexPartition<IndexType>(p_ta_container->size()).for_each(Vector{number_of_data_components}, [&p_ta_container, &ta_span, number_of_data_components](const auto Index, auto& rTLS) {
                            const auto data_begin_index = Index * number_of_data_components;
                            std::copy(ta_span.begin() + data_begin_index, ta_span.begin() + data_begin_index + number_of_data_components, rTLS.begin());
                            (p_ta_container->begin() + Index)->SetValue(TENSOR_ADAPTOR_SYNC, rTLS);
                        });

                        // the p_nd_data is now correctly filled. Add it to the
                        // map and then exit the for loop since, the given container expression
                        // is already found.
                        p_model_part->GetCommunicator().SynchronizeVariable(TENSOR_ADAPTOR_SYNC);
                        r_model_part_data.mPointFields[rTensorAdaptorName] = GetNodalNDData(*r_model_part_data.mpPoints, shape);
                        return true;
                    }
                } else if constexpr(std::is_same_v<ta_container_type, ModelPart::ConditionsContainerType>) {
                    // ghost, local and model part containers are the same. so we check only the model part container
                    if (r_model_part_data.mpCells.has_value() &&
                        std::holds_alternative<ModelPart::ConditionsContainerType::Pointer>(r_model_part_data.mpCells.value()) &&
                        &*p_ta_container == &p_model_part->Conditions()) {
                            // now check if the rTensorAdaptorName is found in the mPointFields of this model part
                            CheckDataArrayName(rTensorAdaptorName, Globals::DataLocation::Condition, r_model_part_data);
                            r_model_part_data.mCellFields[rTensorAdaptorName] = tensor_adaptor_type(*p_tensor_adaptor).pGetStorage();
                            return true;
                    }
                } else if constexpr(std::is_same_v<ta_container_type, ModelPart::ElementsContainerType>) {
                    // ghost, local and model part containers are the same. so we check only the model part container
                    if (r_model_part_data.mpCells.has_value() &&
                        std::holds_alternative<ModelPart::ElementsContainerType::Pointer>(r_model_part_data.mpCells.value()) &&
                        &*p_ta_container == &p_model_part->Elements()) {
                            // now check if the rTensorAdaptorName is found in the mPointFields of this model part
                            CheckDataArrayName(rTensorAdaptorName, Globals::DataLocation::Element, r_model_part_data);
                            r_model_part_data.mCellFields[rTensorAdaptorName] = tensor_adaptor_type(*p_tensor_adaptor).pGetStorage();
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

    KRATOS_CATCH("");
}

template<class TXmlDataElementWrapper>
std::pair<std::string, std::string> VtuOutput::WriteUnstructuredGridData(
    const std::string& rOutputFileNamePrefix,
    UnstructuredGridData& rUnstructuredGridData,
    TXmlDataElementWrapper& rXmlDataElementWrapper) const
{
    const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
    auto p_nodes = rUnstructuredGridData.mpPoints;

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
        std::to_string(rUnstructuredGridData.mpCells.has_value()
                           ? std::visit([](auto v) { return v->size(); },
                                        rUnstructuredGridData.mpCells.value())
                           : 0u));
    unstructured_grid_element->AddElement(piece_element);

    // create the position element
    KRATOS_INFO_IF("VtuOutput", mEchoLevel > 2) << "--- Collecting nodal " << (mIsInitialConfiguration ? "initial" : "current") << " position data...\n";
    NodePositionTensorAdaptor node_position_tensor_adaptor(
        p_nodes, mIsInitialConfiguration ? Globals::Configuration::Initial
                                         : Globals::Configuration::Current);
    node_position_tensor_adaptor.CollectData();

    // create the points element
    auto points_element = Kratos::make_shared<XmlElementsArray>("Points");
    points_element->AddElement(rXmlDataElementWrapper.Get("Position", node_position_tensor_adaptor.pGetStorage()));
    piece_element->AddElement(points_element);

    // create the cells element
    auto cells_element = Kratos::make_shared<XmlElementsArray>("Cells");
    std::set<IndexType> ignored_indices;
    if (rUnstructuredGridData.mpCells.has_value()) {
        KRATOS_INFO_IF("VtuOutput", mEchoLevel > 2) << "--- Collecting " << GetEntityName(rUnstructuredGridData.mpCells) << " connectivity data...\n";
        AddConnectivityData(*cells_element, ignored_indices, *rUnstructuredGridData.mpPoints, rUnstructuredGridData.mpCells.value(), rXmlDataElementWrapper, mEchoLevel);
    }
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
    piece_element->AddElement(point_data_element);

    if (rUnstructuredGridData.UsePointsForDataFieldOutput) {
        KRATOS_INFO_IF("VtuOutput", mEchoLevel > 2) << "--- Collecting nodal data fields...\n";
        // generate and add point field data
        std::set<IndexType> empty_ignored_indices;
        AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeNonHistorical, mFlags), rXmlDataElementWrapper, r_data_communicator, empty_ignored_indices);
        AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeNonHistorical, mVariables), rXmlDataElementWrapper, r_data_communicator, empty_ignored_indices);
        AddFieldsFromTensorAdaptor<HistoricalVariableTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeHistorical, mVariables), rXmlDataElementWrapper, r_data_communicator, empty_ignored_indices, 0);
        AddFields(*point_data_element, rUnstructuredGridData.mPointFields, rXmlDataElementWrapper, empty_ignored_indices);
    }

    // create cell data
    auto cell_data_element = Kratos::make_shared<XmlElementsArray>("CellData");
    piece_element->AddElement(cell_data_element);

    // generate and add cell field data
    if (rUnstructuredGridData.mpCells.has_value()) {
        KRATOS_INFO_IF("VtuOutput", mEchoLevel > 2) << "--- Collecting " << GetEntityName(rUnstructuredGridData.mpCells) << " data fields...\n";
        std::visit([this, &rUnstructuredGridData, &cell_data_element, &rXmlDataElementWrapper, &r_data_communicator, &ignored_indices](auto p_container){
            AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*cell_data_element, p_container, GetContainerMap(this->mVariables, rUnstructuredGridData.mpCells.value()), rXmlDataElementWrapper, r_data_communicator, ignored_indices);
            AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*cell_data_element, p_container, GetContainerMap(this->mFlags, rUnstructuredGridData.mpCells.value()), rXmlDataElementWrapper, r_data_communicator, ignored_indices);
        }, rUnstructuredGridData.mpCells.value());
        AddFields(*cell_data_element, rUnstructuredGridData.mCellFields, rXmlDataElementWrapper, ignored_indices);
    }

    std::stringstream output_vtu_file_name;
    output_vtu_file_name << rOutputFileNamePrefix << "/" << rUnstructuredGridData.mpModelPart->FullName();

    // identify suffix with the entity type.
    const std::string& suffix = "_" + GetEntityName(rUnstructuredGridData.mpCells) + "s";

    const std::string pvd_data_set_name = rUnstructuredGridData.mpModelPart->FullName() + suffix;

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
    UnstructuredGridData& rUnstructuredGridData,
    TXmlDataElementWrapper& rXmlDataElementWrapper) const
{
    if (!rUnstructuredGridData.mpCells.has_value()) {
        // nothing to do here.
        return std::make_pair("", "");
    }

    const auto& integration_point_vars = GetContainerMap(
        mIntegrationPointVariables, rUnstructuredGridData.mpCells.value());

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
    }, rUnstructuredGridData.mpCells.value());

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
    }, rUnstructuredGridData.mpCells.value());

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
        std::visit([&offsets, &point_data_element, &rXmlDataElementWrapper, &rUnstructuredGridData, &is_gauss_point_data_available, total_gauss_points](auto pVariable, auto pContainer) {
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
                        *pVariable, output, rUnstructuredGridData.mpModelPart->GetProcessInfo());

                    if (!output.empty()) {
                        // if there are available gauss point information
                        data_type_traits::Shape(output.front(), local_shape.data(),
                                                local_shape.data() + local_shape.size());
                    }
                }

                // now do the communication between ranks to get the correct size
                const auto& max_local_shape = rUnstructuredGridData.mpModelPart->GetCommunicator().GetDataCommunicator().MaxAll(local_shape);

                // now we construct the nd_data_shape
                DenseVector<unsigned int> nd_data_shape(local_shape.size() + 1);
                std::copy(max_local_shape.begin(), max_local_shape.end(), nd_data_shape.begin() + 1);
                nd_data_shape[0] = total_gauss_points;

                const auto total_number_of_components = std::accumulate(max_local_shape.begin(), max_local_shape.end(), 1u, std::multiplies<unsigned int>{});

                auto p_gauss_data = Kratos::make_shared<NDData<double>>(nd_data_shape);
                auto span = p_gauss_data->ViewData();

                if (span.size() > 0) {
                    is_gauss_point_data_available = true;
                    IndexPartition<IndexType>(pContainer->size()).for_each(output, [&span, &pContainer, &pVariable, &rUnstructuredGridData, &offsets, total_number_of_components](const auto Index, auto& rTLS) {
                        auto& r_entity = *(pContainer->begin() + Index);
                        r_entity.CalculateOnIntegrationPoints(*pVariable, rTLS, rUnstructuredGridData.mpModelPart->GetProcessInfo());
                        DataTypeTraits<std::vector<data_type>>::CopyToContiguousData(span.begin() + offsets[Index] * total_number_of_components, rTLS);
                    });
                    point_data_element->AddElement(rXmlDataElementWrapper.Get(pVariable->Name(), p_gauss_data));
                }
            }
        }, r_pair.second, rUnstructuredGridData.mpCells.value());
    }

    const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();

    if (r_data_communicator.OrReduceAll(is_gauss_point_data_available)) {
        std::stringstream output_vtu_file_name;
        output_vtu_file_name
            << rOutputFileNamePrefix << "/" << rUnstructuredGridData.mpModelPart->FullName() << "_"
            << GetEntityName(rUnstructuredGridData.mpCells)
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
                rUnstructuredGridData.mpModelPart->FullName() + "_gauss",
                WritePartitionedUnstructuredGridData(
                    *point_data_element, *Kratos::make_shared<XmlElementsArray>(""),
                    output_vtu_file_name.str(), r_data_communicator));
        }

        return std::make_pair(rUnstructuredGridData.mpModelPart->FullName() + "_gauss",
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

    for (auto& r_model_part_data : mListOfUnstructuredGridData) {
        switch (mOutputFormat) {
            case ASCII:
            {
                XmlAsciiNDDataElementWrapper data_element_wrapper{mPrecision};
                // first write the unstructured grid data
                KRATOS_INFO_IF("VtuOutput", mEchoLevel > 1)
                    << "Writing \"" << r_model_part_data.mpModelPart->FullName()
                    << "\" " << GetEntityName(r_model_part_data.mpCells) << " fields in ASCII format...\n";
                pvd_file_name_info.push_back(WriteUnstructuredGridData(
                    rOutputFileNamePrefix, r_model_part_data, data_element_wrapper));

                // now write the gauss point info
                KRATOS_INFO_IF("VtuOutput", mEchoLevel > 1)
                    << "Writing \"" << r_model_part_data.mpModelPart->FullName()
                    << "\" " << GetEntityName(r_model_part_data.mpCells) << " gauss point fields in ASCII format...\n";
                pvd_file_name_info.push_back(WriteIntegrationPointData(
                    rOutputFileNamePrefix, r_model_part_data, data_element_wrapper));
                break;
            }
            case BINARY:
            {
                XmlBase64BinaryNDDataElementWrapper data_element_wrapper;
                // first write the unstructured grid data
                KRATOS_INFO_IF("VtuOutput", mEchoLevel > 1)
                    << "Writing \"" << r_model_part_data.mpModelPart->FullName()
                    << "\" " << GetEntityName(r_model_part_data.mpCells) << " fields in binary format...\n";
                pvd_file_name_info.push_back(WriteUnstructuredGridData(
                    rOutputFileNamePrefix, r_model_part_data, data_element_wrapper));

                // now write the gauss point info
                KRATOS_INFO_IF("VtuOutput", mEchoLevel > 1)
                    << "Writing \"" << r_model_part_data.mpModelPart->FullName()
                    << "\" " << GetEntityName(r_model_part_data.mpCells) << " gauss point fields in binary format...\n";
                pvd_file_name_info.push_back(WriteIntegrationPointData(
                    rOutputFileNamePrefix, r_model_part_data, data_element_wrapper));
                break;
            }
        }
    }

    // now generate the *.pvd file
    if (mrModelPart.GetCommunicator().MyPID() == 0) {
        KRATOS_INFO_IF("VtuOutput", mEchoLevel > 1)
                    << "Writing \"" << mrModelPart.FullName()
                    << "\" PVD file...\n";
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
         << (mIsInitialConfiguration ? "yes" : "no") << ", echo level = " << mEchoLevel << " ]";
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
    for (const auto& r_model_part_data : mListOfUnstructuredGridData) {
        rOStream << "\n\tModel part: \"" << r_model_part_data.mpModelPart->FullName() << "\""
                 << " with " << GetEntityName(r_model_part_data.mpCells) << "s"
                 << ", used for point fields = " << (r_model_part_data.UsePointsForDataFieldOutput ? "yes" : "no");

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
