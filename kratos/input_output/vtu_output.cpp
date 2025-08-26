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
#include <vector>
#include <iomanip>
#include <numeric>
#include <type_traits>

// External includes

// Project includes
#include "expression/container_data_io.h"
#include "expression/container_expression.h"
#include "expression/literal_flat_expression.h"
#include "expression/variable_expression_io.h"
#include "includes/data_communicator.h"
#include "input_output/base_64_encoded_output.h"
#include "input_output/vtk_definitions.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/string_utilities.h"
#include "utilities/xml_utilities/xml_expression_element.h"
#include "utilities/xml_utilities/xml_ostream_ascii_writer.h"
#include "utilities/xml_utilities/xml_ostream_base64_binary_writer.h"

#include "utilities/container_io_utils.h"
#include "tensor_adaptors/node_position_tensor_adaptor.h"
#include "tensor_adaptors/variable_tensor_adaptor.h"
#include "tensor_adaptors/historical_variable_tensor_adaptor.h"
#include "tensor_adaptors/flags_tensor_adaptor.h"
#include "tensor_adaptors/gauss_point_variable_tensor_adaptor.h"

#include "utilities/data_type_traits.h"

#include "utilities/xml_utilities/xml_elements_array.h"
#include "utilities/xml_utilities/xml_ascii_nd_data_element.h"
#include "utilities/xml_utilities/xml_base64_binary_nd_data_element.h"

// Include base h
#include "vtu_output.h"

namespace Kratos {

    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, NODES,   1);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, CONDITIONS,  2);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, ELEMENTS, 3);

namespace {

class Xml
{
public:
    Xml()
    {

    }


private:

    template <class TDataType, class... TArgs>
    XmlElement::Pointer CreateDataArrayElement(
        const VtuOutput::WriterFormat& WriterFormat,
        TArgs&&... rArgs)
    {
        switch (WriterFormat) {
            case VtuOutput::WriterFormat::ASCII:
                return Kratos::make_shared<XmlAsciiNDDataElement<TDataType>>(rArgs...);
            case VtuOutput::WriterFormat::BINARY:
                return Kratos::make_shared<XmlBase64BinaryNDDataElement<TDataType>>(rArgs...);
            default:
                KRATOS_ERROR << "Unsupported writer format.";
                return nullptr;
        }
    }


};

template <class TDataType, class... TArgs>
XmlElement::Pointer CreateDataArrayElement(
    const VtuOutput::WriterFormat& WriterFormat,
    TArgs&&... rArgs)
{
    switch (WriterFormat) {
        case VtuOutput::WriterFormat::ASCII:
            return Kratos::make_shared<XmlAsciiNDDataElement<TDataType>>(rArgs...);
        case VtuOutput::WriterFormat::BINARY:
            return Kratos::make_shared<XmlBase64BinaryNDDataElement<TDataType>>(rArgs...);
        default:
            KRATOS_ERROR << "Unsupported writer format.";
            return nullptr;
    }
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

template<class TContainerType>
XmlElement::Pointer CreateCellsXmlElement(
    const VtuOutput::WriterFormat& rWriterFormat,
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    auto p_cells_xml_element = Kratos::make_shared<XmlElementsArray>("Cells");

    auto p_offsets = GetOffsets(rContainer);

    p_cells_xml_element->AddElement(CreateDataArrayElement<int>(rWriterFormat, "connectivity", GetConnectivities(*p_offsets, rContainer, rKratosVtuIndicesMap)));
    p_cells_xml_element->AddElement(CreateDataArrayElement<int>(rWriterFormat, "offsets", p_offsets));
    p_cells_xml_element->AddElement(CreateDataArrayElement<unsigned char>(rWriterFormat, "types", GetGeometryTypes(rContainer)));

    return p_cells_xml_element;
}

template<class TTensorAdaptorType, class TContainerPointerType, class TMapType, class... TArgs>
void AddFieldsFromTensorAdaptor(
    const VtuOutput::WriterFormat& rWriterFormat,
    XmlElementsArray& rXmlElement,
    TContainerPointerType pContainer,
    const TMapType& rMap,
    TArgs&&... rArgs)
{
    KRATOS_TRY

    for (const auto& r_pair : rMap) {
        using data_type = BareType<decltype(r_pair.second)>;
        if constexpr(std::is_same_v<data_type, Flags>) {
            // the map is of type flags
            TTensorAdaptorType tensor_adaptor(pContainer, *r_pair.second, rArgs...);
            tensor_adaptor.Check();
            tensor_adaptor.CollectData();
            rXmlElement.AddElement(CreateDataArrayElement<int>(rWriterFormat, r_pair.first, tensor_adaptor.pGetStorage()));
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
                    rXmlElement.AddElement(CreateDataArrayElement<int>(rWriterFormat, r_pair.first, p_nd_data));
                } else {
                    TTensorAdaptorType tensor_adaptor(pContainer, p_variable, rArgs...);
                    tensor_adaptor.Check();
                    tensor_adaptor.CollectData();
                    rXmlElement.AddElement(CreateDataArrayElement<double>(rWriterFormat, r_pair.first, tensor_adaptor.pGetStorage()));
                }
            }, r_pair.second);
        }
    }

    KRATOS_CATCH("");
}

void AddFields(
    const VtuOutput::WriterFormat& rWriterFormat,
    XmlElementsArray& rXmlElement,
    const std::map<std::string, VtuOutput::FieldPointerType>& rMap)
{
    for (const auto& r_pair : rMap) {
        std::visit([&rXmlElement, &r_pair, &rWriterFormat](auto pNDData) {
            using data_type = typename BareType<decltype(*pNDData)>::DataType;
            rXmlElement.AddElement(CreateDataArrayElement<data_type>(rWriterFormat, r_pair.first, pNDData));
        }, r_pair.second);
    }
}

// template<template <class> class TXmlDataArrayType, class TCellsContainerPointerType, class... TArgs>
// XmlElement::Pointer GetVtkXmlElement(
//     ModelPart::NodesContainerType::Pointer pNodes,
//     TCellsContainerPointerType pCells,
//     const Globals::Configuration& rConfiguration,
//     const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap,
//     const std::unordered_map<std::string, const Flags*>& rNodalFlagsMap,
//     const std::unordered_map<std::string, VtuOutput::SupportedVariables>& rHistoricalVariablesMap,
//     const std::unordered_map<std::string, VtuOutput::SupportedVariables>& rNodalVariablesMap,
//     const std::unordered_map<std::string, const Flags*>& rCellFlagsMap,
//     const std::unordered_map<std::string, VtuOutput::SupportedVariables>& rCellVariablesMap,
//     const std::unordered_map<std::string, VtuOutput::SupportedVariables>& rCellGaussPointVariablesMap,
//     TArgs&&... rArgs)
// {
//     KRATOS_TRY

//     // create the vtk file
//     auto vtk_file_element = Kratos::make_shared<XmlElementsArray>("VTKFile");
//     vtk_file_element->AddAttribute("type", "UnstructuredGrid");
//     vtk_file_element->AddAttribute("version", "0.1");
//     vtk_file_element->AddAttribute("byte_order", VtuOutputHelperUtilities::GetEndianness());

//     // create the unstructured grid
//     auto unstructured_grid_element = Kratos::make_shared<XmlElementsArray>("UnstructuredGrid");
//     vtk_file_element->AddElement(unstructured_grid_element);

//     // create the piece element
//     auto piece_element = Kratos::make_shared<XmlElementsArray>("Piece");
//     piece_element->AddAttribute("NumberOfPoints", std::to_string(pNodes.size()));
//     piece_element->AddAttribute("NumberOfCells", std::to_string(pCells.size()));
//     unstructured_grid_element->AddElement(piece_element);

//     // create the position element
//     NodePositionTensorAdaptor node_position_tensor_adaptor(pNodes, rConfiguration);
//     node_position_tensor_adaptor.Check();
//     node_position_tensor_adaptor.CollectData();

//     // create the points element
//     auto points_element = Kratos::make_shared<XmlElementsArray>("Points");
//     points_element->AddElement(Kratos::make_shared<TXmlDataArrayType<double>>("Position", node_position_tensor_adaptor.pGetStorage(), rArgs...));
//     piece_element->AddElement(points_element);

//     // create the cells element
//     auto cells_element = CreateCellsXmlElement(*pCells, rKratosVtuIndicesMap);
//     piece_element->AddElement(cells_element);

//     // create the point data
//     auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
//     piece_element->AddElement(point_data_element);

//     // generate and add point field data
//     AddFields(*point_data_element, mPointFieldsMap);
//     AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*point_data_element, mrModelPart.pNodes(), mNodalFlagsMap);
//     AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*point_data_element, mrModelPart.pNodes(), mNonHistoricalNodalVariablesMap);
//     AddFieldsFromTensorAdaptor<HistoricalVariableTensorAdaptor>(*point_data_element, mrModelPart.pNodes(), mHistoricalVariablesMap);

//     // create cell data
//     auto cell_data_element = Kratos::make_shared<XmlElementsArray>("CellData");
//     piece_element->AddElement(cell_data_element);

//     // generate and add cell field data
//     AddFields(*cell_data_element, mCellFieldsMap);
//     AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*point_data_element, pCells, mCellFlagsMap);
//     AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*point_data_element, pCells, mNonHistoricalCellVariablesMap);
//     AddFieldsFromTensorAdaptor<GaussPointVariableTensorAdaptor>(*point_data_element, pCells, mGaussPointVariablesMap, mrModelPart.pGetProcessInfo());

//     return vtk_file_element;

//     KRATOS_CATCH("");
// }

} // namespace

namespace VtuOutputHelperUtilities {

Expression::ConstPointer CreatePositionsExpression(
    const ModelPart::NodesContainerType& rNodes,
    const bool IsInitialConfiguration)
{
    auto p_position_expression = LiteralFlatExpression<double>::Create(rNodes.size(), {3});
    auto& r_position_expression = *p_position_expression;

    if (IsInitialConfiguration) {
        IndexPartition<IndexType>(rNodes.size()).for_each([&r_position_expression, &rNodes](const IndexType Index) {
            const auto& r_coordinates = (rNodes.begin() + Index)->GetInitialPosition().Coordinates();
            const IndexType start_index = Index * 3;
            r_position_expression.SetData(start_index, 0, r_coordinates[0]);
            r_position_expression.SetData(start_index, 1, r_coordinates[1]);
            r_position_expression.SetData(start_index, 2, r_coordinates[2]);
        });
    } else {
        IndexPartition<IndexType>(rNodes.size()).for_each([&r_position_expression, &rNodes](const IndexType Index) {
            const auto& r_coordinates = (rNodes.begin() + Index)->Coordinates();
            const IndexType start_index = Index * 3;
            r_position_expression.SetData(start_index, 0, r_coordinates[0]);
            r_position_expression.SetData(start_index, 1, r_coordinates[1]);
            r_position_expression.SetData(start_index, 2, r_coordinates[2]);
        });
    }

    return p_position_expression;
}

XmlExpressionElement::Pointer CreateDataArrayElement(
    const std::string& rDataArrayName,
    const std::vector<const Expression*>& rExpressions)
{
    std::vector<Expression::ConstPointer> expressions;
    for (const auto& p_expression : rExpressions) {
        if (p_expression) {
            expressions.push_back(p_expression);
        }
    }
    if (expressions.size() > 0) {
        return Kratos::make_shared<XmlExpressionElement>(rDataArrayName, expressions);
    } else {
        return nullptr;
    }
}

template<class TContainerType, MeshType TMeshType>
const Expression* pGetExpression(const ContainerExpression<TContainerType, TMeshType>& rContainerExpression)
{
    if (rContainerExpression.HasExpression()) {
        return &*rContainerExpression.pGetExpression();
    } else {
        return nullptr;
    }
}

void AddDataArrayElement(
    XmlExpressionElement::Pointer pParentElement,
    XmlExpressionElement::Pointer pDataArrayElement)
{
    if (pDataArrayElement) {
        pParentElement->AddElement(pDataArrayElement);
    }
}

XmlExpressionElement::Pointer CreatePointsXmlElement(
    const ModelPart& rModelPart,
    const bool IsInitialConfiguration)
{
    auto p_points_xml_element = Kratos::make_shared<XmlExpressionElement>("Points");

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_local_nodes = r_communicator.LocalMesh().Nodes();
    const auto& r_ghost_nodes = r_communicator.GhostMesh().Nodes();

    auto local_position = CreatePositionsExpression(r_local_nodes, IsInitialConfiguration);
    auto ghost_position = CreatePositionsExpression(r_ghost_nodes, IsInitialConfiguration);

    AddDataArrayElement(
        p_points_xml_element,
        CreateDataArrayElement("Position", {&*local_position, &*ghost_position}));

    return p_points_xml_element;
}

template<class TContainerType>
LiteralFlatExpression<int>::Pointer CreateOffsetsExpression(const TContainerType& rContainer)
{
    auto p_offsets_expression = LiteralFlatExpression<int>::Create(rContainer.size(), {});
    auto data_itr = p_offsets_expression->begin();

    IndexType total_offset = 0;
    for (const auto& r_entity : rContainer) {
        total_offset += r_entity.GetGeometry().size();
        *(data_itr++) = total_offset;
    }
    return p_offsets_expression;
}

template<class TContainerType>
Expression::Pointer CreateGeometryTypesExpression(const TContainerType& rContainer)
{
    auto p_geometry_expression = LiteralFlatExpression<char>::Create(rContainer.size(), {});
    auto data_itr = p_geometry_expression->begin();

    IndexPartition<IndexType>(rContainer.size()).for_each([data_itr, &rContainer](const IndexType Index) {
        const auto p_itr = VtkDefinitions::KratosVtkGeometryTypes.find((rContainer.begin() + Index)->GetGeometry().GetGeometryType());
        if (p_itr != VtkDefinitions::KratosVtkGeometryTypes.end()) {
            *(data_itr + Index) = p_itr->second;
        } else {
            KRATOS_ERROR << "Element with id " << (rContainer.begin() + Index)->Id() << " has unsupported geometry.";
        }
    });
    return p_geometry_expression;
}

template<class TContainerType>
Expression::Pointer CreateConnectivityExpression(
    const LiteralFlatExpression<int>::Pointer pOffsetsExpression,
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    auto offset_data_itr = pOffsetsExpression->begin();

    auto p_connectivity_expression = LiteralFlatExpression<int>::Create(*(offset_data_itr + rContainer.size() - 1), {});
    auto data_itr = p_connectivity_expression->begin();

    IndexPartition<IndexType>(rContainer.size()).for_each([data_itr, offset_data_itr, &rContainer, &rKratosVtuIndicesMap](const IndexType Index) {
        const auto& r_geometry = (rContainer.begin() + Index)->GetGeometry();
        auto entity_data_begin_itr = data_itr + *(offset_data_itr + Index) - r_geometry.size();

        for (const auto& r_node : r_geometry) {
            const auto p_itr = rKratosVtuIndicesMap.find(r_node.Id());
            if (p_itr != rKratosVtuIndicesMap.end()) {
                *(entity_data_begin_itr++) = p_itr->second;
            } else {
                KRATOS_ERROR << "Node with id " << r_node.Id() << " not found in nodes list.";
            }
        }
    });
    return p_connectivity_expression;
}

template<class TContainerType>
XmlExpressionElement::Pointer CreateCellsXmlElement(
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    auto p_cells_xml_element = Kratos::make_shared<XmlExpressionElement>("Cells");

    auto p_offsets_expression = CreateOffsetsExpression(rContainer);
    auto p_connectivity_expression = CreateConnectivityExpression(p_offsets_expression, rContainer, rKratosVtuIndicesMap);
    auto p_geometry_type_expression = CreateGeometryTypesExpression(rContainer);

    AddDataArrayElement(
        p_cells_xml_element,
        CreateDataArrayElement("connectivity", {&*p_connectivity_expression}));
    AddDataArrayElement(p_cells_xml_element,
                        CreateDataArrayElement("offsets", {&*p_offsets_expression}));
    AddDataArrayElement(p_cells_xml_element,
                        CreateDataArrayElement("types", {&*p_geometry_type_expression}));
    return p_cells_xml_element;
}

template<class TDataType, class TContainerType, class TContainerDataIOTag>
XmlExpressionElement::Pointer CreateVariableDataXmlElement(
    const Variable<TDataType>& rVariable,
    ModelPart& rModelPart)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        ContainerExpression<TContainerType> local_container(rModelPart);
        if constexpr(std::is_same_v<TContainerDataIOTag, ContainerDataIOTags::Historical>) {
            VariableExpressionIO::Read(local_container, &rVariable, true);
        } else {
            VariableExpressionIO::Read(local_container, &rVariable, false);
        }

        if (rModelPart.GetCommunicator().GhostMesh().NumberOfNodes() > 0) {
            ContainerExpression<TContainerType, MeshType::Ghost> ghost_container(rModelPart);
            if constexpr(std::is_same_v<TContainerDataIOTag, ContainerDataIOTags::Historical>) {
                VariableExpressionIO::Read(ghost_container, &rVariable, true);
            } else {
                VariableExpressionIO::Read(ghost_container, &rVariable, false);
            }

            return CreateDataArrayElement(
                rVariable.Name(),
                {pGetExpression(local_container), pGetExpression(ghost_container)});
        } else {
            return CreateDataArrayElement(
                rVariable.Name(),
                {pGetExpression(local_container)});
        }

    } else {
        ContainerExpression<TContainerType> local_container(rModelPart);
        VariableExpressionIO::Read(local_container, &rVariable);

        return CreateDataArrayElement(rVariable.Name(), {pGetExpression(local_container)});
    }
}

template<class TContainerType>
Expression::Pointer CreateContainerFlagExpression(
    const TContainerType& rContainer,
    const Flags& rFlag)
{
    auto p_flag_expression = LiteralFlatExpression<int>::Create(rContainer.size(), {});
    auto data_itr = p_flag_expression->begin();

    IndexPartition<IndexType>(rContainer.size()).for_each([data_itr, &rContainer, &rFlag](const IndexType Index) {
        *(data_itr + Index) = (rContainer.begin() + Index)->Is(rFlag);
    });

    return p_flag_expression;
}

template<class TContainerType>
XmlExpressionElement::Pointer CreateFlagDataXmlElement(
    const std::string& rFlagName,
    const Flags& rFlags,
    const ModelPart& rModelPart)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        const auto& r_communicator = rModelPart.GetCommunicator();
        const auto& r_local_nodes = r_communicator.LocalMesh().Nodes();
        const auto& r_ghost_nodes = r_communicator.GhostMesh().Nodes();

        return CreateDataArrayElement(
            rFlagName,
            {&*CreateContainerFlagExpression(r_local_nodes, rFlags),
             &*CreateContainerFlagExpression(r_ghost_nodes, rFlags)});
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return CreateDataArrayElement(
            rFlagName,
            {&*CreateContainerFlagExpression(rModelPart.GetCommunicator().LocalMesh().Conditions(), rFlags)});
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return CreateDataArrayElement(
            rFlagName,
            {&*CreateContainerFlagExpression(rModelPart.GetCommunicator().LocalMesh().Elements(), rFlags)});
    }
}

Expression::Pointer CreateGhostNodeExpression(
    const DataCommunicator& rDataCommunicator,
    const Expression& rLocalNodesExpression,
    const ModelPart::NodesContainerType& rLocalNodes,
    const ModelPart::NodesContainerType& rGhostNodes,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    const IndexType number_of_ghost_nodes = rGhostNodes.size();

    std::vector<int> ghost_indices(number_of_ghost_nodes);
    std::transform(rGhostNodes.begin(), rGhostNodes.end(), ghost_indices.begin(), [](const auto& rNode) { return rNode.Id(); });
    auto gp_list = GlobalPointerUtilities::RetrieveGlobalIndexedPointers(rLocalNodes, ghost_indices, rDataCommunicator);

    GlobalPointerCommunicator<ModelPart::NodeType> pointer_comm(rDataCommunicator, gp_list.ptr_begin(), gp_list.ptr_end());

    const IndexType number_of_components = rLocalNodesExpression.GetItemComponentCount();

    auto values_proxy = pointer_comm.Apply(
        [&rLocalNodesExpression, number_of_components, &rKratosVtuIndicesMap](GlobalPointer<ModelPart::NodeType>& rGP) -> std::vector<double> {
            std::vector<double> values(number_of_components);
            const auto p_itr = rKratosVtuIndicesMap.find(rGP->Id());
            if (p_itr != rKratosVtuIndicesMap.end()) {
                const IndexType enitity_data_begin_index = p_itr->second * number_of_components;
                for (IndexType i = 0; i < number_of_components; ++i) {
                    values[i] = rLocalNodesExpression.Evaluate(p_itr->second, enitity_data_begin_index, i);
                }
            } else {
                KRATOS_ERROR << "The node with id " << rGP->Id() << " not found in the owning rank local expression.";
            }
            return values;
        }
    );

    auto p_ghost_nodes_expression = LiteralFlatExpression<double>::Create(number_of_ghost_nodes, rLocalNodesExpression.GetItemShape());
    auto data_itr = p_ghost_nodes_expression->begin();

    for(IndexType i = 0; i < number_of_ghost_nodes; ++i) {
        const auto& r_gp_value = values_proxy.Get(gp_list(i));
        for (IndexType j = 0; j < number_of_components; ++j) {
            *(data_itr++) = r_gp_value[j];
        }
    }
    return p_ghost_nodes_expression;
}

template<class TContainerType>
XmlExpressionElement::Pointer CreateContainerExpressionXmlElement(
    const std::string& rExpressionName,
    const ContainerExpression<TContainerType>& rContainerExpression,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        const auto& r_communicator = rContainerExpression.GetModelPart().GetCommunicator();
        const auto& r_local_nodes = r_communicator.LocalMesh().Nodes();
        const auto& r_ghost_nodes = r_communicator.GhostMesh().Nodes();

        auto ghost_node_expression = CreateGhostNodeExpression(
            r_communicator.GetDataCommunicator(), rContainerExpression.GetExpression(),
            r_local_nodes, r_ghost_nodes, rKratosVtuIndicesMap);

        return CreateDataArrayElement(
            rExpressionName, {pGetExpression(rContainerExpression), &*ghost_node_expression});
    } else {
        return CreateDataArrayElement(rExpressionName,
                                      {pGetExpression(rContainerExpression)});
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

template<class TContaierType, class TContainerIOTag>
void AddListOfVariables(
    XmlExpressionElement::Pointer pXmlDataElement,
    ModelPart& rModelPart,
    const std::map<std::string, VtuOutput::SupportedVariables>& rVariablesMap)
{
    for (const auto& p_variable_variant : rVariablesMap){
        std::visit([&rModelPart, &pXmlDataElement](auto pVariable) {
                using variable_type =
                    typename std::remove_const_t<std::remove_reference_t<decltype(*pVariable)>>::Type;
                AddDataArrayElement(
                    pXmlDataElement,
                    CreateVariableDataXmlElement<variable_type, TContaierType, TContainerIOTag>(
                        *pVariable, rModelPart));
            },
            p_variable_variant.second);
    }
}

template<class TContaierType>
void AddListOfFlags(
    XmlExpressionElement::Pointer pXmlDataElement,
    ModelPart& rModelPart,
    const std::map<std::string, const Flags*>& rFlagsMap)
{
    for (const auto& r_flag_pair : rFlagsMap){
        AddDataArrayElement(pXmlDataElement, CreateFlagDataXmlElement<TContaierType>(r_flag_pair.first, *r_flag_pair.second, rModelPart));
    }
}

template<class TContainerExpressionType>
void AddListOfContainerExpressions(
    XmlExpressionElement::Pointer pXmlDataElement,
    ModelPart& rModelPart,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap,
    const std::map<std::string, TContainerExpressionType>& rContainerExpressionsMap)
{
    for (const auto& r_container_expression_pair : rContainerExpressionsMap){
        const std::string& r_name = r_container_expression_pair.first;
        if constexpr(std::is_same_v<TContainerExpressionType, ContainerExpression<ModelPart::NodesContainerType>::Pointer>) {
            AddDataArrayElement(
                pXmlDataElement,
                CreateContainerExpressionXmlElement(
                    r_name,
                    *r_container_expression_pair.second, rKratosVtuIndicesMap));
        } else {
            std::visit([&r_name, &pXmlDataElement, &rKratosVtuIndicesMap](auto pContainerExpression) {
                    AddDataArrayElement(
                        pXmlDataElement,
                        CreateContainerExpressionXmlElement(
                            r_name, *pContainerExpression, rKratosVtuIndicesMap));
                },
                r_container_expression_pair.second);
        }
    }
}

void CopyAttributes(
    XmlExpressionElement& rOutputElement,
    const XmlExpressionElement& rInputElement)
{
    for (const auto& r_attribute_data : rInputElement.GetAttributes()) {
        rOutputElement.AddAttribute(r_attribute_data.first, r_attribute_data.second);
    }
}

void CreatePDataArrays(
    XmlExpressionElement& rOutputElement,
    const XmlExpressionElement& rInputElement)
{
    for (const auto& p_element : rInputElement.GetElements()) {
        auto new_pdata_array = Kratos::make_shared<XmlExpressionElement>("PDataArray");
        CopyAttributes(*new_pdata_array, *p_element);
        rOutputElement.AddElement(new_pdata_array);
    }
}

std::string GetEndianness()
{
    int i = 0x0001;

    if (*reinterpret_cast<char*>(&i) != 0) {
        return "LittleEndian";
    } else {
        return "BigEndian";
    }
}

void WritePvtuFile(
    XmlExpressionElement::Pointer pVtkFileElement,
    const ModelPart& rModelPart,
    const std::string& rOutputFileNamePrefix,
    const std::string& rOutputFileName)
{
    // get list of file names
    std::stringstream list_of_file_names;

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_data_communicator = r_communicator.GetDataCommunicator();

    // TODO: May be we want to check a rank which has some entities (not empty ranks)
    //       Then write on that rank.
    const int writing_rank = 0;

    list_of_file_names << rOutputFileName << "\n";
    if (r_data_communicator.Rank() == writing_rank) {
        for (int rank = 0; rank < r_data_communicator.Size(); ++rank) {
            if (rank != writing_rank) {
                std::string msg;
                r_data_communicator.Recv(msg, rank);
                list_of_file_names << msg;
            }
        }
    } else {
        r_data_communicator.Send(list_of_file_names.str(), writing_rank);
    }
    r_data_communicator.Barrier();

    if (r_data_communicator.Rank() == writing_rank) {
        // create the vtk file
        auto vtk_file_element = Kratos::make_shared<XmlExpressionElement>("VTKFile");
        vtk_file_element->AddAttribute("type", "PUnstructuredGrid");
        vtk_file_element->AddAttribute("version", "0.1");
        vtk_file_element->AddAttribute("byte_order", GetEndianness());

        // create the unstructured grid
        auto unstructured_grid_element = Kratos::make_shared<XmlExpressionElement>("PUnstructuredGrid");
        unstructured_grid_element->AddAttribute("GhostLevel", "0");
        vtk_file_element->AddElement(unstructured_grid_element);

        // get the points xml element
        auto piece = pVtkFileElement->GetElements("UnstructuredGrid")[0]->GetElements("Piece")[0];

        auto points = piece->GetElements("Points")[0];
        auto p_points = Kratos::make_shared<XmlExpressionElement>("PPoints");
        VtuOutputHelperUtilities::CreatePDataArrays(*p_points, *points);
        unstructured_grid_element->AddElement(p_points);

        auto cells = piece->GetElements("Cells")[0];
        auto p_cells = Kratos::make_shared<XmlExpressionElement>("PCells");
        VtuOutputHelperUtilities::CreatePDataArrays(*p_cells, *cells);
        unstructured_grid_element->AddElement(p_cells);

        auto point_data = piece->GetElements("PointData")[0];
        auto p_point_data = Kratos::make_shared<XmlExpressionElement>("PPointData");
        VtuOutputHelperUtilities::CreatePDataArrays(*p_point_data, *point_data);
        unstructured_grid_element->AddElement(p_point_data);

        auto cell_data = piece->GetElements("CellData")[0];
        auto p_cell_data = Kratos::make_shared<XmlExpressionElement>("PCellData");
        VtuOutputHelperUtilities::CreatePDataArrays(*p_cell_data, *cell_data);
        unstructured_grid_element->AddElement(p_cell_data);

        // now write all the pieces
        const auto& r_file_names = StringUtilities::SplitStringByDelimiter(list_of_file_names.str(), '\n');
        for (const auto& r_file_name : r_file_names) {
            auto piece = Kratos::make_shared<XmlExpressionElement>("Piece");
            piece->AddAttribute("Source", r_file_name);
            unstructured_grid_element->AddElement(piece);
        }

        // writing to file
        std::stringstream output_pvtu_file_name;
        output_pvtu_file_name << rOutputFileNamePrefix << ".pvtu";
        std::ofstream output_file;
        output_file.open(output_pvtu_file_name.str(), std::ios::out | std::ios::trunc);
        XmlOStreamAsciiWriter writer(output_file, 1);
        writer.WriteElement(*vtk_file_element);
        output_file.close();
    }
}

}; // namespace VtuOutputHelperUtilities

VtuOutput::VtuOutput(
    ModelPart& rModelPart,
    const bool IsInitialConfiguration,
    const WriterFormat OutputFormat,
    const IndexType Precision)
    : mrModelPart(rModelPart),
      mIsInitialConfiguration(IsInitialConfiguration),
      mOutputFormat(OutputFormat),
      mPrecision(Precision)
{
    const auto& r_communicator = rModelPart.GetCommunicator();

    mIsConditionsConsidered = r_communicator.GlobalNumberOfConditions() > 0;
    mIsElementsConsidered = r_communicator.GlobalNumberOfElements() > 0;

    KRATOS_WARNING_IF("VtuOutput", mIsElementsConsidered && mIsConditionsConsidered)
        << "Conditions and Elements vtu output chosen for " << mrModelPart.FullName()
        << " which is not supported. Giving priority to elements.\n";

    mIsConditionsConsidered = mIsElementsConsidered ? false : mIsConditionsConsidered;

    if (mIsConditionsConsidered || mIsElementsConsidered) {
        // we first always use the local mesh
        IndexType vtu_index = 0;
        for (const auto& r_node : mrModelPart.GetCommunicator().LocalMesh().Nodes()) {
            mKratosVtuIndicesMap[r_node.Id()] = vtu_index++;
        }

        // then we add the ghost mesh
        for (const auto& r_node : mrModelPart.GetCommunicator().GhostMesh().Nodes()) {
            mKratosVtuIndicesMap[r_node.Id()] = vtu_index++;
        }
    }
}

template<class TDataType>
void VtuOutput::AddHistoricalVariable(const Variable<TDataType>& rVariable)
{
    VtuOutputHelperUtilities::CheckDataArrayName(
        rVariable.Name(), mNonHistoricalNodalVariablesMap, mNodalFlagsMap, mPointContainerExpressionsMap);
    mHistoricalVariablesMap[rVariable.Name()] = &rVariable;
}

template<class TDataType>
void VtuOutput::AddNonHistoricalVariable(
    const Variable<TDataType>& rVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        VtuOutputHelperUtilities::CheckDataArrayName(
            rVariable.Name(), mHistoricalVariablesMap, mNodalFlagsMap,
            mPointContainerExpressionsMap);
        mNonHistoricalNodalVariablesMap[rVariable.Name()] = &rVariable;
    } else {
        KRATOS_ERROR_IF_NOT(!rEntityFlags.Is(CONDITIONS) || (rEntityFlags.Is(CONDITIONS) && mIsConditionsConsidered))
            << "Condition variable \""
            << rVariable.Name() << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";

        KRATOS_ERROR_IF_NOT(!rEntityFlags.Is(ELEMENTS) || (rEntityFlags.Is(ELEMENTS) && mIsElementsConsidered))
            << "Element variable \""
            << rVariable.Name() << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";

        VtuOutputHelperUtilities::CheckDataArrayName(
            rVariable.Name(), mCellFlagsMap, mCellContainerExpressionsMap);
        mNonHistoricalCellVariablesMap[rVariable.Name()] = &rVariable;
    }
}

void VtuOutput::AddFlagVariable(
    const std::string& rFlagName,
    const Flags& rFlagVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        VtuOutputHelperUtilities::CheckDataArrayName(
            rFlagName, mHistoricalVariablesMap, mNonHistoricalNodalVariablesMap,
            mPointContainerExpressionsMap);
        mNodalFlagsMap[rFlagName] = &rFlagVariable;
    } else {
        KRATOS_ERROR_IF_NOT(!rEntityFlags.Is(CONDITIONS) || (rEntityFlags.Is(CONDITIONS) && mIsConditionsConsidered))
            << "Condition flag \""
            << rFlagName << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";

        KRATOS_ERROR_IF_NOT(!rEntityFlags.Is(ELEMENTS) || (rEntityFlags.Is(ELEMENTS) && mIsElementsConsidered))
            << "Element flag \""
            << rFlagName << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";

        VtuOutputHelperUtilities::CheckDataArrayName(
            rFlagName, mNonHistoricalCellVariablesMap, mCellContainerExpressionsMap);
        mCellFlagsMap[rFlagName] = &rFlagVariable;
    }
}

template <class TContainerType>
void VtuOutput::AddContainerExpression(
    const std::string& rExpressionName,
    const typename ContainerExpression<TContainerType>::Pointer pContainerExpression)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        KRATOS_ERROR_IF_NOT(mIsConditionsConsidered)
            << "Conditions container expression \"" << rExpressionName
            << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        KRATOS_ERROR_IF_NOT(mIsElementsConsidered)
            << "Elements container expression \"" << rExpressionName
            << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";
    }

    KRATOS_ERROR_IF_NOT(&pContainerExpression->GetModelPart() == &mrModelPart)
        << "Model part mismatch in container expression addition. [ Vtu output model part name = \""
        << mrModelPart.FullName() << "\", container expression model part name = \""
        << pContainerExpression->GetModelPart().FullName() << "\" ].\n";

    if constexpr (std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        VtuOutputHelperUtilities::CheckDataArrayName(
            rExpressionName, mHistoricalVariablesMap,
            mNonHistoricalNodalVariablesMap, mNodalFlagsMap);
        mPointContainerExpressionsMap[rExpressionName] = pContainerExpression;
    } else {
        VtuOutputHelperUtilities::CheckDataArrayName(
            rExpressionName, mNonHistoricalCellVariablesMap, mCellFlagsMap);
        mCellContainerExpressionsMap[rExpressionName] = pContainerExpression;
    }
}

void VtuOutput::PrintModelPart(
    const std::string& rOutputFileNamePrefix,
    ModelPart& rModelPart) const
{
    // create the vtk file
    XmlElementsArray vtk_file_element("VTKFile");
    vtk_file_element.AddAttribute("type", "UnstructuredGrid");
    vtk_file_element.AddAttribute("version", "0.1");
    vtk_file_element.AddAttribute("byte_order", VtuOutputHelperUtilities::GetEndianness());

    // create the unstructured grid
    auto unstructured_grid_element = Kratos::make_shared<XmlElementsArray>("UnstructuredGrid");
    vtk_file_element.AddElement(unstructured_grid_element);

    // create the piece element
    auto piece_element = Kratos::make_shared<XmlElementsArray>("Piece");
    piece_element->AddAttribute("NumberOfPoints", std::to_string(rModelPart.NumberOfNodes()));
    piece_element->AddAttribute("NumberOfCells",
        std::to_string(mIsElementsConsidered ? rModelPart.NumberOfElements()
                                             : rModelPart.NumberOfConditions()));
    unstructured_grid_element->AddElement(piece_element);

    // create the position element
    NodePositionTensorAdaptor node_position_tensor_adaptor(rModelPart.pNodes(), mIsInitialConfiguration ? Globals::Configuration::Initial : Globals::Configuration::Current);
    node_position_tensor_adaptor.Check();
    node_position_tensor_adaptor.CollectData();

    // create the points element
    auto points_element = Kratos::make_shared<XmlElementsArray>("Points");
    points_element->AddElement(CreateDataArrayElement<double>(mOutputFormat, "Position", node_position_tensor_adaptor.pGetStorage()));
    piece_element->AddElement(points_element);

    // create the cells element
    auto cells_element = mIsElementsConsidered ? CreateCellsXmlElement(mOutputFormat, rModelPart.Elements(), mKratosVtuIndicesMap) : CreateCellsXmlElement(mOutputFormat, rModelPart.Conditions(), mKratosVtuIndicesMap);
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
    piece_element->AddElement(point_data_element);

    // generate and add point field data
    AddFields(mOutputFormat, *point_data_element, mPointFieldsMap);
    AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pNodes(), mNodalFlagsMap);
    AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pNodes(), mNonHistoricalNodalVariablesMap);
    AddFieldsFromTensorAdaptor<HistoricalVariableTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pNodes(), mHistoricalVariablesMap);

    // create cell data
    auto cell_data_element = Kratos::make_shared<XmlElementsArray>("CellData");
    piece_element->AddElement(cell_data_element);

    // generate and add point field data
    AddFields(mOutputFormat, *cell_data_element, mCellFieldsMap);
    if (mIsElementsConsidered) {
        AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pElements(), mCellFlagsMap);
        AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pElements(), mNonHistoricalCellVariablesMap);
        AddFieldsFromTensorAdaptor<GaussPointVariableTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pElements(), mGaussPointVariablesMap, mrModelPart.pGetProcessInfo());
    } else if (mIsConditionsConsidered) {
        AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pConditions(), mCellFlagsMap);
        AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pConditions(), mNonHistoricalCellVariablesMap);
        AddFieldsFromTensorAdaptor<GaussPointVariableTensorAdaptor>(mOutputFormat, *point_data_element, mrModelPart.pConditions(), mGaussPointVariablesMap, mrModelPart.pGetProcessInfo());
    }

    const auto& r_communicator = rModelPart.GetCommunicator();

    std::stringstream output_vtu_file_name;
    output_vtu_file_name << rOutputFileNamePrefix;
    if (r_communicator.IsDistributed()) {
        output_vtu_file_name << "_" << r_communicator.MyPID();
    }
    output_vtu_file_name << ".vtu";

    std::ofstream output_file;
    output_file.open(output_vtu_file_name.str(), std::ios::out | std::ios::trunc);

    vtk_file_element.Write(output_file);
    output_file.close();

    // if (r_communicator.IsDistributed()) {
    //     VtuOutputHelperUtilities::WritePvtuFile(vtk_file_element, rModelPart, rOutputFileNamePrefix,
    //                                             output_vtu_file_name.str());
    // }
}

void VtuOutput::ClearHistoricalVariables()
{
    mHistoricalVariablesMap.clear();
}

void VtuOutput::ClearNodalNonHistoricalVariables()
{
    mNonHistoricalNodalVariablesMap.clear();
}

void VtuOutput::ClearCellNonHistoricalVariables()
{
    mNonHistoricalCellVariablesMap.clear();
}

void VtuOutput::ClearNodalFlags()
{
    mNodalFlagsMap.clear();
}

void VtuOutput::ClearCellFlags()
{
    mCellFlagsMap.clear();
}

void VtuOutput::ClearNodalContainerExpressions()
{
    mPointContainerExpressionsMap.clear();
}

void VtuOutput::ClearCellContainerExpressions()
{
    mCellContainerExpressionsMap.clear();
}

const ModelPart& VtuOutput::GetModelPart() const
{
    return mrModelPart;
}

void VtuOutput::PrintOutput(const std::string& rOutputFilenamePrefix)
{
    PrintModelPart(rOutputFilenamePrefix, mrModelPart);
}

// template instantiations
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<int>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<double>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 3>>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 4>>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 6>>&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddHistoricalVariable(const Variable<array_1d<double, 9>>&);

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<int>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<double>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 3>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 4>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 6>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddNonHistoricalVariable(const Variable<array_1d<double, 9>>&, const Flags&);

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::NodesContainerType>(const std::string&, const typename ContainerExpression<ModelPart::NodesContainerType>::Pointer);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ConditionsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ConditionsContainerType>::Pointer);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ElementsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ElementsContainerType>::Pointer);

} // namespace Kratos
