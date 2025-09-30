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
#include <filesystem>
#include <numeric>
#include <type_traits>

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
#include "future/utilities/xml_utilities/xml_appended_data_element_wrapper.h"
#include "future/utilities/xml_utilities/xml_elements_array.h"
#include "future/utilities/xml_utilities/xml_in_place_data_element_wrapper.h"
#include "future/utilities/xml_utilities/xml_utils.h"

// Include base h
#include "vtu_output.h"

namespace Kratos::Future {

namespace {

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

void CopyAttributes(
    const XmlElement& rSource,
    XmlElement& rDestination)
{
    for (const auto& [attribute, value] : rSource.GetAttributes()) {
        rDestination.AddAttribute(attribute, value);
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
    if (rOffsets.Size() == 0) {
        return Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, 0));
    }

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
        auto p_type_data = GetGeometryTypes(rIgnoredIndices, *p_container, EchoLevel);

        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2) << "------ Ignored " << rIgnoredIndices.size() << "/"
            << p_container->size() << " " << GetEntityName(p_container) << "(s).\n";

        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2) << "------ Collecting geometry offsets info...\n";
        auto p_offsets = GetOffsets(rIgnoredIndices, *p_container);

        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2) << "------ Collecting geometry connectivity info...\n";
        rCellElement.AddElement(rXmlDataElementWrapper.Get("connectivity", GetConnectivities(*p_offsets, *p_container, indices_map, rIgnoredIndices)));
        rCellElement.AddElement(rXmlDataElementWrapper.Get("offsets", p_offsets));
        rCellElement.AddElement(rXmlDataElementWrapper.Get("types", GetNonIgnoredNDData(rIgnoredIndices, p_type_data)));
    }, pCells);
}

template<class TTensorAdaptorType, class TContainerPointerType, class TDataType, class TXmlDataElementWrapper, class... TArgs>
void AddFieldsFromTensorAdaptorImpl(
    XmlElementsArray& rXmlElement,
    TContainerPointerType pContainer,
    const Variable<TDataType>& rVariable,
    TXmlDataElementWrapper& rXmlDataElementWrapper,
    const DataCommunicator& rDataCommunicator,
    const std::set<IndexType>& rIgnoredIndices,
    const IndexType EchoLevel,
    TArgs&&... rArgs)
{
    KRATOS_TRY

    using primitive_data_type = typename DataTypeTraits<TDataType>::PrimitiveType;

    if constexpr(std::is_same_v<primitive_data_type, int>) {
        // we only support int variable, and there are no TensorAdaptors to
        // collect data from integer variables, we do it manually here.
        auto p_nd_data = Kratos::make_shared<NDData<int>>(DenseVector<unsigned int>(1, pContainer->size()));
        const auto& shape = p_nd_data->Shape();
        if constexpr (std::is_same_v<TTensorAdaptorType, HistoricalVariableTensorAdaptor>) {
            ContainerIOUtils::CopyToContiguousArray<int>(
                *pContainer, p_nd_data->ViewData(),
                shape.data().begin(), shape.data().begin() + 1,
                [&rVariable](int& rValue, const Node& rNode) {
                    rValue = rNode.FastGetSolutionStepValue(rVariable);
                });
        } else if constexpr(std::is_same_v<TTensorAdaptorType, VariableTensorAdaptor>) {
            ContainerIOUtils::CopyToContiguousArray<int>(
                *pContainer, p_nd_data->ViewData(),
                shape.data().begin(), shape.data().begin() + 1,
                [&rVariable](int& rValue, const auto& rNode) {
                    rValue = rNode.GetValue(rVariable);
                });
        } else {
            KRATOS_ERROR << "Unsupported tensor adaptor type.";
        }

        // since we only support Variable<int>, which is having a static data shape
        // we don't have to do mpi communication to decide the shape on the
        // empty ranks.
        rXmlElement.AddElement(rXmlDataElementWrapper.Get(rVariable.Name(), GetNonIgnoredNDData(rIgnoredIndices, p_nd_data)));
    } else if constexpr(std::is_same_v<primitive_data_type, double>) {
        using container_type = BareType<decltype(*pContainer)>;
        using entity_type = typename container_type::data_type;
        using data_type_traits = DataTypeTraits<TDataType>;

        if constexpr(data_type_traits::IsDynamic) {
            // this is a dynamic type such as Vector or Matrix, so
            // we need to do communication to decide the correct size
            // because, there may be ranks with zero entities, hence these ranks will not have
            // correctly sized dynamic Matrix or Vectors. Vtu needs all the ranks to have the same number
            // of components, hence communication is a must in here.

            // construct the correct data_shape
            std::vector<unsigned int> data_shape(data_type_traits::Dimension, 0);
            if (!pContainer->empty()) {
                TDataType value;
                if constexpr(std::is_same_v<TTensorAdaptorType, HistoricalVariableTensorAdaptor>) {
                    value = static_cast<const entity_type&>(pContainer->front()).FastGetSolutionStepValue(rVariable);
                } else if constexpr(std::is_same_v<TTensorAdaptorType, VariableTensorAdaptor>) {
                    value = static_cast<const entity_type&>(pContainer->front()).GetValue(rVariable);
                } else {
                    KRATOS_ERROR << "Unsupported tensor adaptor type.";
                    value = TDataType{};
                }
                data_type_traits::Shape(value, data_shape.data(), data_shape.data() + data_type_traits::Dimension);
            }
            const auto& max_data_shape = rDataCommunicator.MaxAll(data_shape);
            DenseVector<unsigned int> nd_shape(max_data_shape.size() + 1);
            std::copy(max_data_shape.begin(), max_data_shape.end(), nd_shape.begin() + 1);
            nd_shape[0] = pContainer->size();

            // construct the data storage
            auto p_nd_data = Kratos::make_shared<typename TTensorAdaptorType::Storage>(nd_shape);
            auto base_ta = TensorAdaptor<double>(pContainer, p_nd_data, false);
            TTensorAdaptorType tensor_adaptor(base_ta, &rVariable, rArgs..., false);
            tensor_adaptor.CollectData();

            rXmlElement.AddElement(rXmlDataElementWrapper.Get(rVariable.Name(), GetNonIgnoredNDData(rIgnoredIndices, tensor_adaptor.pGetStorage())));
        } else {
            // this is a static type such as double, array_1d<double, 3>, ...
            // So no need of mpi communication
            TTensorAdaptorType tensor_adaptor(pContainer, &rVariable, rArgs...);
            tensor_adaptor.CollectData();
            rXmlElement.AddElement(rXmlDataElementWrapper.Get(rVariable.Name(), GetNonIgnoredNDData(rIgnoredIndices, tensor_adaptor.pGetStorage())));
        }
    } else {
        KRATOS_ERROR << "Unsupported variable type.";
    }

    KRATOS_CATCH("");
}

template<class TTensorAdaptorType, class TContainerPointerType, class TMapType, class TXmlDataElementWrapper, class... TArgs>
void AddFieldsFromTensorAdaptor(
    XmlElementsArray& rXmlElement,
    TContainerPointerType pContainer,
    const TMapType& rMap,
    TXmlDataElementWrapper& rXmlDataElementWrapper,
    const DataCommunicator& rDataCommunicator,
    const std::set<IndexType>& rIgnoredIndices,
    const IndexType EchoLevel,
    TArgs&&... rArgs)
{
    KRATOS_TRY

    for (const auto& r_pair : rMap) {
        using data_type = BareType<decltype(r_pair.second)>;

        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2) << "------ Collecting " << r_pair.first << " data...\n";

        if constexpr(std::is_same_v<data_type, Flags>) {
            // the map is of type flags
            // here we don't need to do any communication because Flags are always having a static data shape.
            TTensorAdaptorType tensor_adaptor(pContainer, *r_pair.second, rArgs...);
            tensor_adaptor.CollectData();
            rXmlElement.AddElement(rXmlDataElementWrapper.Get(r_pair.first, GetNonIgnoredNDData(rIgnoredIndices, tensor_adaptor.pGetStorage())));
        } else {
            std::visit([&](const auto p_variable) {
                AddFieldsFromTensorAdaptorImpl<TTensorAdaptorType>(
                    rXmlElement, pContainer, *p_variable, rXmlDataElementWrapper,
                    rDataCommunicator, rIgnoredIndices, EchoLevel, rArgs...);
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
    const std::set<IndexType>& rIgnoredIndices,
    const IndexType EchoLevel)
{
    for (const auto& r_pair : rMap) {
        KRATOS_INFO_IF("VtuOutput", EchoLevel > 2) << "------ Collecting " << r_pair.first << " data...\n";

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
    const IndexType EchoLevel,
    const bool UseSubModelParts)
{
    const std::vector<char> entity_availability{rModelPart.NumberOfNodes() > 0, rModelPart.NumberOfConditions() > 0, rModelPart.NumberOfElements() > 0};
    const auto& max_entity_availability = rModelPart.GetRootModelPart().GetCommunicator().GetDataCommunicator().MaxAll(entity_availability);
    const bool has_nodes = max_entity_availability[0];
    const bool has_conditions = max_entity_availability[1];
    const bool has_elements   = max_entity_availability[2];

    if (has_elements) {
        // Model part has elements. Hence add a separate output
        // for elements.

        // now check if it has proper nodes
        if (has_nodes) {
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

        if (!has_elements && has_nodes) {
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

    if (UseSubModelParts) {
        // now recursively add all the sub model part data.
        for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
            AddUnstructuredGridData(rOutput, r_sub_model_part, EchoLevel, UseSubModelParts);
        }
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
    const int writing_rank = 0;

    // remove the rank from the rOutputVtuFileName.
    const auto& r_base_name = rOutputVtuFileName.substr(0, rOutputVtuFileName.rfind("_"));

    const auto& p_vtu_file_name = r_base_name  + ".pvtu";

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
        p_types_element->AddAttribute("type", "UInt8");
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
        for (IndexType i_rank = 0; i_rank < rDataCommunicator.Size(); ++i_rank) {
            const auto& r_file_name = r_base_name + "_" + std::to_string(i_rank) + ".vtu";
            auto piece = Kratos::make_shared<XmlElementsArray>("Piece");
            // since we are writing to the same folder the pvtu files
            piece->AddAttribute(
                "Source", std::filesystem::relative(
                              std::filesystem::absolute(r_file_name),
                              std::filesystem::absolute(p_vtu_file_name).parent_path())
                              .generic_string());
            p_unstructured_grid_element->AddElement(piece);
        }

        // writing to file
        std::ofstream output_file;
        output_file.open(p_vtu_file_name, std::ios::out | std::ios::trunc);
        p_vtu_file_element.Write(output_file);

   }

    return p_vtu_file_name;
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

template<class TContainerType>
std::pair<Globals::DataLocation, std::vector<VtuOutput::UnstructuredGridData>::iterator> FindUnstructuredGridData(
    const TContainerType& rContainer,
    std::vector<VtuOutput::UnstructuredGridData>& rUnstructuredGridDataList)
{
    for (auto itr = rUnstructuredGridDataList.begin(); itr != rUnstructuredGridDataList.end(); ++itr) {
        auto p_model_part = itr->mpModelPart;
        if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
            if (itr->UsePointsForDataFieldOutput) {
                if (
                    &rContainer == &p_model_part->Nodes() ||
                    (!p_model_part->IsDistributed() && &rContainer == &p_model_part->GetCommunicator().LocalMesh().Nodes())
                ) {
                    return std::make_pair(Globals::DataLocation::NodeNonHistorical, itr);
                }
            }
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
            // since Kratos is doing partitioning based on elements, and there are no conditions
            // on the ghost meshes. so normal mesh and the local mesh should be having identical entities.
            if (itr->mpCells.has_value() &&
                std::holds_alternative<ModelPart::ConditionsContainerType::Pointer>(itr->mpCells.value()) &&
                (
                    &rContainer == &*std::get<ModelPart::ConditionsContainerType::Pointer>(itr->mpCells.value()) ||
                    &rContainer == &p_model_part->GetCommunicator().LocalMesh().Conditions()
                )
            ) {
                return std::make_pair(Globals::DataLocation::Condition, itr);
            }
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
            // since Kratos is doing partitioning based on elements, and there are no elements
            // on the ghost meshes. so normal mesh and the local mesh should be having identical entities.
            if (itr->mpCells.has_value() &&
                std::holds_alternative<ModelPart::ElementsContainerType::Pointer>(itr->mpCells.value()) &&
                (
                    &rContainer == &*std::get<ModelPart::ElementsContainerType::Pointer>(itr->mpCells.value()) ||
                    &rContainer == &p_model_part->GetCommunicator().LocalMesh().Elements()
                )
            ) {
                return std::make_pair(Globals::DataLocation::Element, itr);
            }
        } else {
            KRATOS_ERROR << "Unsupported container type.";
        }
    }

    return std::make_pair(Globals::DataLocation::ModelPart, rUnstructuredGridDataList.end());
}

} // namespace

VtuOutput::VtuOutput(
    ModelPart& rModelPart,
    const Globals::Configuration Configuration,
    const WriterFormat OutputFormat,
    const IndexType Precision,
    const bool OutputSubModelParts,
    const IndexType EchoLevel)
    : mIsPVDFileHeaderWritten(false),
      mrModelPart(rModelPart),
      mConfiguration(Configuration),
      mEchoLevel(EchoLevel),
      mOutputFormat(OutputFormat),
      mPrecision(Precision)
{
    AddUnstructuredGridData(mUnstructuredGridDataList, rModelPart, mEchoLevel, OutputSubModelParts);

    // sort the order of output to be consistent between different compilers
    std::sort(mUnstructuredGridDataList.begin(), mUnstructuredGridDataList.end(),
              [](const auto& rV1, const auto& rV2) {
                  return rV1.mpModelPart->FullName() < rV2.mpModelPart->FullName();
              });
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

    KRATOS_ERROR_IF(mIsPVDFileHeaderWritten)
        << "Flags can be added only before the first call to the PrintOutput [ flag name = "
        << rFlagName << " ].\n";

    switch (DataLocation) {
        case Globals::DataLocation::NodeNonHistorical:
            // additionally check here in the historical containers.
            CheckDataArrayName(rFlagName, {Globals::DataLocation::NodeHistorical}, mFlags, mVariables);
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            CheckDataArrayName(rFlagName, {DataLocation}, mFlags, mVariables); // checks in the current data location of mVariables map
            CheckDataArrayName(rFlagName, DataLocation, mUnstructuredGridDataList); // checks in the tensor adaptors list
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

    KRATOS_ERROR_IF(mIsPVDFileHeaderWritten)
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
            CheckDataArrayName(GetName(pVariable), {DataLocation}, mFlags, mVariables); // checks in the current data location of mVariables map
            CheckDataArrayName(GetName(pVariable), DataLocation, mUnstructuredGridDataList); // checks in the tensor adaptors list
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

    KRATOS_ERROR_IF(mIsPVDFileHeaderWritten)
        << "Integration point variables can be added only before the first call to the PrintOutput [ integration point variable name = "
        << GetName(pVariable) << " ].\n";

    switch (DataLocation) {
        case Globals::DataLocation::Condition:
        case Globals::DataLocation::Element:
            CheckDataArrayName(GetName(pVariable), {DataLocation}, mIntegrationPointVariables); // checks in the current data location of mIntegrationPointVariables map
            mIntegrationPointVariables[DataLocation][GetName(pVariable)] = pVariable;
            break;
        default:
            KRATOS_ERROR << "Integration point variables can be only added to Condition, and Element data locations.";
            break;
    }

    KRATOS_CATCH("");
}

void VtuOutput::AddTensorAdaptor(
    const std::string& rTensorAdaptorName,
    SupportedTensorAdaptorPointerType pTensorAdaptor)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mIsPVDFileHeaderWritten)
        << "TensorAdaptors can be added only before the first call to the PrintOutput [ tensor adaptor name = "
        << rTensorAdaptorName << " ].\n";

    std::visit([this, &rTensorAdaptorName](auto p_tensor_adaptor) {
        using tensor_adaptor_type = BareType<decltype(*p_tensor_adaptor)>;
        using storage_type = typename tensor_adaptor_type::Storage;
        auto shape = p_tensor_adaptor->Shape();
        auto ta_span = p_tensor_adaptor->ViewData();

        // the tensor adaptors may have different number of components in dimensions from 2 to N [Eg.
        // tensor adaptors created from dynamic variables, having zero entities. ] Hence doing
        // communication to correctly size the number of components in higher dimensions.
        std::vector<unsigned int> local_data_shape(shape.begin() + 1, shape.end());
        const auto& max_data_shape = this->mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_data_shape);
        if (shape[0] == 0) std::copy(max_data_shape.begin(), max_data_shape.end(), shape.begin() + 1);

        std::visit([this, &rTensorAdaptorName, &p_tensor_adaptor, &ta_span, &shape](auto pContainer){
            [[maybe_unused]] auto [mesh_type, itr] = FindUnstructuredGridData(*pContainer, this->mUnstructuredGridDataList);

            switch (mesh_type) {
                case Globals::DataLocation::Condition:
                case Globals::DataLocation::Element: {
                    auto& r_unstructured_grid_data = *itr;
                    // check for existing field names
                    CheckDataArrayName(rTensorAdaptorName, {mesh_type}, mFlags, mVariables); // checks in the current data location of mVariables map
                    CheckDataArrayName(rTensorAdaptorName, mesh_type, r_unstructured_grid_data); // checks in the tensor adaptors list
                    r_unstructured_grid_data.mCellFields[rTensorAdaptorName] = Kratos::make_shared<storage_type>(ta_span.data(), shape);
                    break;
                }
                case Globals::DataLocation::NodeNonHistorical: {
                    auto& r_unstructured_grid_data = *itr;
                    // check for existing field names
                    CheckDataArrayName(rTensorAdaptorName, {Globals::DataLocation::NodeNonHistorical, Globals::DataLocation::NodeHistorical}, mFlags, mVariables); // checks in the current data location of mVariables map
                    CheckDataArrayName(rTensorAdaptorName, Globals::DataLocation::NodeNonHistorical, r_unstructured_grid_data); // checks in the tensor adaptors list
                    r_unstructured_grid_data.mPointFields[rTensorAdaptorName] = Kratos::make_shared<storage_type>(ta_span.data(), shape);
                    break;
                }
                default:
                    KRATOS_ERROR
                        << "The container in the TensorAdaptor is not referring to any of the containers "
                        << "written by this Vtu output [ tensor adaptor name = " << rTensorAdaptorName
                        << ", tensor_adaptor = " << *p_tensor_adaptor << " ]\n"
                        << *this;
                    break;
            }
        }, p_tensor_adaptor->GetContainer());
    }, pTensorAdaptor);

    KRATOS_CATCH("");
}

void VtuOutput::UpdateTensorAdaptor(
    const std::string& rTensorAdaptorName,
    SupportedTensorAdaptorPointerType pTensorAdaptor)
{
    KRATOS_TRY

    std::visit([this, &rTensorAdaptorName](auto p_tensor_adaptor) {
        using tensor_adaptor_type = BareType<decltype(*p_tensor_adaptor)>;
        using storage_type = typename tensor_adaptor_type::Storage;
        auto shape = p_tensor_adaptor->Shape();
        auto ta_span = p_tensor_adaptor->ViewData();

        // the tensor adaptors may have different number of components in dimensions from 2 to N [Eg.
        // tensor adaptors created from dynamic variables, having zero entities. ] Hence doing
        // communication to correctly size the number of components in higher dimensions.
        std::vector<unsigned int> local_data_shape(shape.begin() + 1, shape.end());
        const auto& max_data_shape = this->mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_data_shape);
        if (shape[0] == 0) std::copy(max_data_shape.begin(), max_data_shape.end(), shape.begin() + 1);

        std::visit([this, &rTensorAdaptorName, &p_tensor_adaptor, &ta_span, &shape](auto pContainer){
            [[maybe_unused]] auto [mesh_type, itr] = FindUnstructuredGridData(*pContainer, this->mUnstructuredGridDataList);

            switch (mesh_type) {
                case Globals::DataLocation::Condition:
                case Globals::DataLocation::Element: {
                    auto data_field_itr = itr->mCellFields.find(rTensorAdaptorName);

                    KRATOS_ERROR_IF(data_field_itr == itr->mCellFields.end())
                        << "TensorAdaptor name = \""
                        << rTensorAdaptorName << "\" not found in the existing data fields. It is only allowed to update existing data fields. VtuOutput: \n"
                        << *this;

                    data_field_itr->second = Kratos::make_shared<storage_type>(ta_span.data(), shape);
                    break;
                }
                case Globals::DataLocation::NodeNonHistorical: {
                    auto data_field_itr = itr->mPointFields.find(rTensorAdaptorName);

                    KRATOS_ERROR_IF(data_field_itr == itr->mPointFields.end())
                        << "TensorAdaptor name = \""
                        << rTensorAdaptorName << "\" not found in the existing data fields. It is only allowed to update existing data fields. VtuOutput: \n"
                        << *this;

                    data_field_itr->second = Kratos::make_shared<storage_type>(ta_span.data(), shape);
                    break;
                }
                default:
                    KRATOS_ERROR
                        << "The container in the TensorAdaptor is not referring to any of the containers "
                        << "written by this Vtu output [ tensor adaptor name = " << rTensorAdaptorName
                        << ", tensor_adaptor = " << *p_tensor_adaptor << " ]\n"
                        << *this;
                    break;
            }
        }, p_tensor_adaptor->GetContainer());
    }, pTensorAdaptor);

    KRATOS_CATCH("");
}

void VtuOutput::EmplaceTensorAdaptor(
    const std::string& rTensorAdaptorName,
    SupportedTensorAdaptorPointerType pTensorAdaptor)
{
    if (!mIsPVDFileHeaderWritten) {
        AddTensorAdaptor(rTensorAdaptorName, pTensorAdaptor);
    } else {
        UpdateTensorAdaptor(rTensorAdaptorName, pTensorAdaptor);
    }
}

template<class TXmlElementDataWrapperCreateFunctor, class TXmlElementDataWrapperAppendFunctor>
std::pair<std::string, std::string> VtuOutput::WriteUnstructuredGridData(
    TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor,
    TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor,
    UnstructuredGridData& rUnstructuredGridData,
    const std::string& rOutputFileNamePrefix,
    const IndexType Step) const
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

    auto p_xml_data_element_wrapper = rElementDataWrapperCreateFunctor();
    rElementDataWrapperAppendFunctor(vtk_file_element, p_xml_data_element_wrapper);

    // create the piece element
    auto piece_element = Kratos::make_shared<XmlElementsArray>("Piece");
    // adding number of points
    piece_element->AddAttribute("NumberOfPoints", std::to_string(p_nodes->size()));
    unstructured_grid_element->AddElement(piece_element);

    // create the position element
    KRATOS_INFO_IF("VtuOutput", mEchoLevel > 2) << "--- Collecting nodal " << (mConfiguration == Globals::Configuration::Initial ? "initial" : "current") << " position data...\n";
    NodePositionTensorAdaptor node_position_tensor_adaptor(p_nodes, mConfiguration);
    node_position_tensor_adaptor.CollectData();

    // create the points element
    auto points_element = Kratos::make_shared<XmlElementsArray>("Points");
    points_element->AddElement(p_xml_data_element_wrapper->Get("Position", node_position_tensor_adaptor.pGetStorage()));
    piece_element->AddElement(points_element);

    // create the cells element
    auto cells_element = Kratos::make_shared<XmlElementsArray>("Cells");
    std::set<IndexType> ignored_indices;
    if (rUnstructuredGridData.mpCells.has_value()) {
        KRATOS_INFO_IF("VtuOutput", mEchoLevel > 2) << "--- Collecting " << GetEntityName(rUnstructuredGridData.mpCells) << " connectivity data...\n";
        AddConnectivityData(*cells_element, ignored_indices, *rUnstructuredGridData.mpPoints, rUnstructuredGridData.mpCells.value(), *p_xml_data_element_wrapper, mEchoLevel);
        // adding number of cells
        piece_element->AddAttribute(
            "NumberOfCells",
            std::to_string(std::visit([&ignored_indices](auto v) { return v->size() - ignored_indices.size(); }, rUnstructuredGridData.mpCells.value())));
    } else {
        // adding number of cells
        piece_element->AddAttribute("NumberOfCells", "0");
    }
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
    piece_element->AddElement(point_data_element);

    if (rUnstructuredGridData.UsePointsForDataFieldOutput) {
        KRATOS_INFO_IF("VtuOutput", mEchoLevel > 2) << "--- Collecting nodal data fields...\n";
        // generate and add point field data
        std::set<IndexType> empty_ignored_indices;
        AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeNonHistorical, mFlags), *p_xml_data_element_wrapper, r_data_communicator, empty_ignored_indices, mEchoLevel);
        AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeNonHistorical, mVariables), *p_xml_data_element_wrapper, r_data_communicator, empty_ignored_indices, mEchoLevel);
        AddFieldsFromTensorAdaptor<HistoricalVariableTensorAdaptor>(*point_data_element, p_nodes, GetUnorderedMapValue(Globals::DataLocation::NodeHistorical, mVariables), *p_xml_data_element_wrapper, r_data_communicator, empty_ignored_indices, mEchoLevel, 0);
        AddFields(*point_data_element, rUnstructuredGridData.mPointFields, *p_xml_data_element_wrapper, empty_ignored_indices, mEchoLevel);
    }

    // create cell data
    auto cell_data_element = Kratos::make_shared<XmlElementsArray>("CellData");
    piece_element->AddElement(cell_data_element);

    // generate and add cell field data
    if (rUnstructuredGridData.mpCells.has_value()) {
        KRATOS_INFO_IF("VtuOutput", mEchoLevel > 2) << "--- Collecting " << GetEntityName(rUnstructuredGridData.mpCells) << " data fields...\n";
        std::visit([this, &rUnstructuredGridData, &cell_data_element, &p_xml_data_element_wrapper, &r_data_communicator, &ignored_indices](auto p_container){
            AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*cell_data_element, p_container, GetContainerMap(this->mVariables, rUnstructuredGridData.mpCells.value()), *p_xml_data_element_wrapper, r_data_communicator, ignored_indices, mEchoLevel);
            AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*cell_data_element, p_container, GetContainerMap(this->mFlags, rUnstructuredGridData.mpCells.value()), *p_xml_data_element_wrapper, r_data_communicator, ignored_indices, mEchoLevel);
        }, rUnstructuredGridData.mpCells.value());
        AddFields(*cell_data_element, rUnstructuredGridData.mCellFields, *p_xml_data_element_wrapper, ignored_indices, mEchoLevel);
    }

    std::stringstream output_vtu_file_name;
    output_vtu_file_name << rOutputFileNamePrefix << "/" << rUnstructuredGridData.mpModelPart->FullName();

    // identify suffix with the entity type.
    const std::string& suffix = "_" + GetEntityName(rUnstructuredGridData.mpCells) + "s";

    const std::string pvd_data_set_name = rUnstructuredGridData.mpModelPart->FullName() + suffix;

    // append with the step value and rank and extension
    output_vtu_file_name << suffix << "_" << Step
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

template<class TXmlElementDataWrapperCreateFunctor, class TXmlElementDataWrapperAppendFunctor>
std::pair<std::string, std::string> VtuOutput::WriteIntegrationPointData(
    TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor,
    TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor,
    UnstructuredGridData& rUnstructuredGridData,
    const std::string& rOutputFileNamePrefix,
    const IndexType Step) const
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

    auto p_xml_data_element_wrapper = rElementDataWrapperCreateFunctor();
    rElementDataWrapperAppendFunctor(vtk_file_element, p_xml_data_element_wrapper);

    DenseVector<std::size_t> offsets;

    const auto total_gauss_points = std::visit([&offsets](auto p_container) {
        // resize the offsets
        offsets.resize(p_container->size(), false);

        IndexType total_number_of_gauss_points = 0;

        // now compute the offsets for each entity. This allows
        // having different number of gps in different entities.
        // which is the case if we have a model part with mixed type
        // of elements.
        for (IndexType i = 0; i < p_container->size(); ++i) {
            const auto& r_entity = *(p_container->begin() + i);
            const auto number_of_gps = r_entity.GetGeometry().IntegrationPointsNumber(r_entity.GetIntegrationMethod());
            total_number_of_gauss_points += number_of_gps;
            offsets[i] = total_number_of_gauss_points - number_of_gps;
        }

        return total_number_of_gauss_points;
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
    points_element->AddElement(p_xml_data_element_wrapper->Get("Position", gauss_point_positions));
    piece_element->AddElement(points_element);

    auto cells_element = Kratos::make_shared<XmlElementsArray>("Cells");
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
    piece_element->AddElement(point_data_element);

    const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();

    // add the gauss point data
    bool is_gauss_point_data_available = false;
    for (const auto& r_pair : integration_point_vars) {
        std::visit([&offsets, &point_data_element, &p_xml_data_element_wrapper, &rUnstructuredGridData, &is_gauss_point_data_available, &r_data_communicator, total_gauss_points](auto pVariable, auto pContainer) {
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
                const auto& max_local_shape = r_data_communicator.MaxAll(local_shape);

                // now we construct the nd_data_shape
                DenseVector<unsigned int> nd_data_shape(local_shape.size() + 1);
                std::copy(max_local_shape.begin(), max_local_shape.end(), nd_data_shape.begin() + 1);
                nd_data_shape[0] = total_gauss_points;

                const auto total_number_of_components = std::accumulate(max_local_shape.begin(), max_local_shape.end(), 1u, std::multiplies<unsigned int>{});

                auto p_gauss_data = Kratos::make_shared<NDData<double>>(nd_data_shape);
                auto span = p_gauss_data->ViewData();

                if (r_data_communicator.SumAll(static_cast<unsigned long>(span.size())) > 0) {
                    is_gauss_point_data_available = true;
                    IndexPartition<IndexType>(pContainer->size()).for_each(std::vector<data_type>{}, [&span, &pContainer, &pVariable, &rUnstructuredGridData, &offsets, total_number_of_components](const auto Index, auto& rTLS) {
                        auto& r_entity = *(pContainer->begin() + Index);
                        r_entity.CalculateOnIntegrationPoints(*pVariable, rTLS, rUnstructuredGridData.mpModelPart->GetProcessInfo());
                        DataTypeTraits<std::vector<data_type>>::CopyToContiguousData(span.begin() + offsets[Index] * total_number_of_components, rTLS);
                    });
                    point_data_element->AddElement(p_xml_data_element_wrapper->Get(pVariable->Name(), p_gauss_data));
               }
           }
        }, r_pair.second, rUnstructuredGridData.mpCells.value());
    }

    if (is_gauss_point_data_available) {
        std::stringstream output_vtu_file_name;
        output_vtu_file_name
            << rOutputFileNamePrefix << "/" << rUnstructuredGridData.mpModelPart->FullName() << "_"
            << GetEntityName(rUnstructuredGridData.mpCells)
            << "_gauss_" << Step
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
                rUnstructuredGridData.mpModelPart->FullName() + "_" + GetEntityName(rUnstructuredGridData.mpCells) + "_gauss",
                WritePartitionedUnstructuredGridData(
                    *point_data_element, *Kratos::make_shared<XmlElementsArray>(""),
                    output_vtu_file_name.str(), r_data_communicator));
        }

        return std::make_pair(rUnstructuredGridData.mpModelPart->FullName() + "_" + GetEntityName(rUnstructuredGridData.mpCells) + "_gauss",
                              output_vtu_file_name.str());
    }
    else {
        return std::make_pair("", "");
    }
}

template<class TXmlElementDataWrapperCreateFunctor, class TXmlElementDataWrapperAppendFunctor>
void VtuOutput::WriteData(
    std::vector<std::pair<std::string, std::string>>& rPVDFileNameInfo,
    TXmlElementDataWrapperCreateFunctor&& rElementDataWrapperCreateFunctor,
    TXmlElementDataWrapperAppendFunctor&& rElementDataWrapperAppendFunctor,
    UnstructuredGridData& rUnstructuredGridData,
    const std::string& rOutputPrefix,
    const IndexType Step) const
{
    KRATOS_TRY

    rPVDFileNameInfo.push_back(WriteUnstructuredGridData(
        rElementDataWrapperCreateFunctor, rElementDataWrapperAppendFunctor,
        rUnstructuredGridData, rOutputPrefix, Step));

    rPVDFileNameInfo.push_back(WriteIntegrationPointData(
        rElementDataWrapperCreateFunctor, rElementDataWrapperAppendFunctor,
        rUnstructuredGridData, rOutputPrefix, Step));

    KRATOS_CATCH("");
}

void VtuOutput::PrintOutput(const std::string& rOutputFileNamePrefix)
{
    KRATOS_TRY

    const auto& r_process_info = mrModelPart.GetProcessInfo();

    // Add the time step info.
    // check if a similar time has been already printed to vtu.
    // here we do not check whether the r_process_info has the time variable specified
    // because, in a const DataValueContainer, if the variable is not there, it returns the
    // zero value of the variable.
    const double time = r_process_info[TIME];

    const IndexType step = r_process_info[STEP];

    std::filesystem::create_directories(rOutputFileNamePrefix);

    std::vector<std::pair<std::string, std::string>> pvd_file_name_info;

    for (auto& r_unstructured_grid_data : mUnstructuredGridDataList) {
        switch (mOutputFormat) {
            case ASCII:
            {
                WriteData(
                    pvd_file_name_info,
                    [this](){ return Kratos::make_shared<XmlInPlaceDataElementWrapper>(XmlInPlaceDataElementWrapper::ASCII, this->mPrecision); },
                    [](auto& rVtkFileElement, auto pXmlDataElementWrapper) {},
                    r_unstructured_grid_data,
                    rOutputFileNamePrefix,
                    step);
                break;
            }
            case BINARY:
            {
                WriteData(
                    pvd_file_name_info,
                    [this](){ return Kratos::make_shared<XmlInPlaceDataElementWrapper>(XmlInPlaceDataElementWrapper::BINARY, this->mPrecision); },
                    [](auto& rVtkFileElement, auto pXmlDataElementWrapper) { rVtkFileElement.AddAttribute("header_type", "UInt64"); },
                    r_unstructured_grid_data,
                    rOutputFileNamePrefix,
                    step);
                break;
            }
            case RAW:
            {
                WriteData(
                    pvd_file_name_info,
                    [](){ return Kratos::make_shared<XmlAppendedDataElementWrapper>(XmlAppendedDataElementWrapper::RAW); },
                    [](auto& rVtkFileElement, auto pXmlDataElementWrapper) { rVtkFileElement.AddElement(pXmlDataElementWrapper); rVtkFileElement.AddAttribute("header_type", "UInt64"); },
                    r_unstructured_grid_data,
                    rOutputFileNamePrefix,
                    step);
                break;
            }
            case COMPRESSED_RAW:
            {
                WriteData(
                    pvd_file_name_info,
                    [](){ return Kratos::make_shared<XmlAppendedDataElementWrapper>(XmlAppendedDataElementWrapper::RAW_COMPRESSED); },
                    [](auto& rVtkFileElement, auto pXmlDataElementWrapper) { rVtkFileElement.AddElement(pXmlDataElementWrapper); rVtkFileElement.AddAttribute("header_type", "UInt64"); rVtkFileElement.AddAttribute("compressor", "vtkZLibDataCompressor"); },
                    r_unstructured_grid_data,
                    rOutputFileNamePrefix,
                    step);
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

        if (!mIsPVDFileHeaderWritten) {
            mIsPVDFileHeaderWritten = true;

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
            IndexType local_index = 0;
            for (IndexType i = 0; i < pvd_file_name_info.size(); ++i) {
                if (pvd_file_name_info[i].second != "") {
                    auto current_element = Kratos::make_shared<XmlElementsArray>("DataSet");

                    // write the time with the specified precision.
                    std::stringstream str_time;
                    str_time << std::scientific << std::setprecision(mPrecision) << time;

                    current_element->AddAttribute("timestep", str_time.str());
                    current_element->AddAttribute("name", pvd_file_name_info[i].first);
                    current_element->AddAttribute("part", std::to_string(local_index++));
                    current_element->AddAttribute(
                        "file", std::filesystem::relative(
                                    std::filesystem::absolute(pvd_file_name_info[i].second),
                                    std::filesystem::absolute(rOutputFileNamePrefix).parent_path())
                                    .generic_string());
                    collection_element->AddElement(current_element);
                }
            }

            std::ofstream output_file;
            output_file.open(rOutputFileNamePrefix + ".pvd", std::ios::out | std::ios::trunc);
            pvd_file_element.Write(output_file);
            output_file.close();
        } else {
            std::ofstream output_file;
            output_file.open(rOutputFileNamePrefix + ".pvd", std::ios::in | std::ios::out);
            output_file.seekp(-28, std::ios::end);

            // now iterate through all the time steps and correctly write
            // the file names for each time step.
            IndexType local_index = 0;
            for (IndexType i = 0; i < pvd_file_name_info.size(); ++i) {
                if (pvd_file_name_info[i].second != "") {
                    auto current_element = Kratos::make_shared<XmlElementsArray>("DataSet");

                    // write the time with the specified precision.
                    std::stringstream str_time;
                    str_time << std::scientific << std::setprecision(mPrecision) << time;

                    current_element->AddAttribute("timestep", str_time.str());
                    current_element->AddAttribute("name", pvd_file_name_info[i].first);
                    current_element->AddAttribute("part", std::to_string(local_index++));
                    current_element->AddAttribute(
                        "file", std::filesystem::relative(
                                    std::filesystem::absolute(pvd_file_name_info[i].second),
                                    std::filesystem::absolute(rOutputFileNamePrefix).parent_path())
                                    .generic_string());

                    current_element->Write(output_file, 2);
                }
            }

            output_file << "   </Collection>\n</VTKFile>";
            output_file.close();
        }

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
        case RAW:
            info << "RAW";
            break;
        case COMPRESSED_RAW:
            info << "COMPRESSED_RAW";
            break;
    }

    info << ", precision = " << mPrecision << ", configuration = ";

    switch (mConfiguration){
        case Globals::Configuration::Initial:
            info << "initial";
            break;
        case Globals::Configuration::Current:
            info << "current";
            break;
    }

    info << ", echo level = " << mEchoLevel << " ]";

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
    for (const auto& r_model_part_data : mUnstructuredGridDataList) {
        rOStream << "\n\tModel part: \"" << r_model_part_data.mpModelPart->FullName() << "\""
                 << " with " << GetEntityName(r_model_part_data.mpCells) << "s"
                 << ", used for point fields = " << (r_model_part_data.UsePointsForDataFieldOutput ? "yes" : "no");

        if (r_model_part_data.UsePointsForDataFieldOutput) {
            rOStream << "\n\t\t" << "Point fields:";
            for (const auto& r_pair : r_model_part_data.mPointFields) {
                std::visit([&rOStream, &r_pair](auto pNDData){
                    rOStream << "\n\t\t\t" << r_pair.first << ": " << *pNDData;
                }, r_pair.second);
            }
        }

        rOStream << "\n\t\t" << "Cell fields:";
        for (const auto& r_pair : r_model_part_data.mCellFields) {
            std::visit([&rOStream, &r_pair](auto pNDData){
                rOStream << "\n\t\t\t" << r_pair.first << ": " << *pNDData;
            }, r_pair.second);
        }
    }
}

} // namespace Kratos
