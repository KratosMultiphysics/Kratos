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
#include <numeric>
#include <type_traits>
#include <filesystem>

// External includes

// Project includes
#include "includes/data_communicator.h"
#include "input_output/base_64_encoded_output.h"
#include "input_output/vtk_definitions.h"
#include "tensor_adaptors/flags_tensor_adaptor.h"
#include "tensor_adaptors/gauss_point_variable_tensor_adaptor.h"
#include "tensor_adaptors/historical_variable_tensor_adaptor.h"
#include "tensor_adaptors/node_position_tensor_adaptor.h"
#include "tensor_adaptors/variable_tensor_adaptor.h"
#include "utilities/container_io_utils.h"
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/xml_utilities/xml_ascii_nd_data_element.h"
#include "utilities/xml_utilities/xml_base64_binary_nd_data_element.h"
#include "utilities/xml_utilities/xml_elements_array.h"

// Include base h
#include "vtu_output.h"
namespace Kratos
{

    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, NODES,   1);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, CONDITIONS,  2);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, ELEMENTS, 3);

namespace
{
struct XmlAsciiNDDataElementWrapper
{
    IndexType mPrecision;

    template<class NDDataTypePointer>
    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataTypePointer pNDData)
    {
        using data_type = typename BareType<decltype(*pNDData)>::DataType;
        return Kratos::make_shared<XmlAsciiNDDataElement<data_type>>(rDataArrayName, pNDData, mPrecision);
    }

    template<class NDDataTypePointer>
    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataTypePointer pNDData,
        const DenseVector<unsigned int>& rShape)
    {
        using data_type = typename BareType<decltype(*pNDData)>::DataType;
        return Kratos::make_shared<XmlAsciiNDDataElement<data_type>>(rDataArrayName, pNDData, rShape, mPrecision);
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
        return Kratos::make_shared<XmlBase64BinaryNDDataElement<data_type>>(rDataArrayName, pNDData);
    }

    template<class NDDataTypePointer>
    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataTypePointer pNDData,
        const DenseVector<unsigned int>& rShape)
    {
        using data_type = typename BareType<decltype(*pNDData)>::DataType;
        return Kratos::make_shared<XmlBase64BinaryNDDataElement<data_type>>(rDataArrayName, pNDData, rShape);
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

std::size_t GetNumberOfCells(const std::tuple<ModelPart*, VtuOutput::CellType, std::unordered_map<IndexType, IndexType>>& rModelPartCellData)
{
    const auto& r_model_part = *std::get<0>(rModelPartCellData);
    const auto cell_type = std::get<1>(rModelPartCellData);

    switch (cell_type) {
        case VtuOutput::CellType::Elements:
            return r_model_part.NumberOfElements();
        case VtuOutput::CellType::Conditions:
            return r_model_part.NumberOfConditions();
        default:
            return 0;
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

template<class TCellsPointerType, class TXmlDataElementWrapper>
XmlElement::Pointer CreateCellsXmlElement(
    TCellsPointerType pCells,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap,
    TXmlDataElementWrapper& rXmlDataElementWrapper)
{
    auto p_cells_xml_element = Kratos::make_shared<XmlElementsArray>("Cells");

    auto p_offsets = GetOffsets(*pCells);
    p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("connectivity", GetConnectivities(*p_offsets, *pCells, rKratosVtuIndicesMap)));
    p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("offsets", p_offsets));
    p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("types", GetGeometryTypes(*pCells)));

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

void FillCellData(
    std::vector<std::tuple<ModelPart*, VtuOutput::CellType, std::unordered_map<IndexType, IndexType>*>>& rOutput,
    ModelPart& rModelPart)
{
    const bool has_elements = rModelPart.GetCommunicator().GlobalNumberOfElements() > 0;
    const bool has_conditions = rModelPart.GetCommunicator().GlobalNumberOfConditions() > 0;

    auto p_indices_map = new std::unordered_map<IndexType, IndexType>{};

    if (has_elements) {
        rOutput.push_back(std::make_tuple(&rModelPart, VtuOutput::CellType::Elements, p_indices_map));
    }

    if (has_conditions) {
        rOutput.push_back(std::make_tuple(&rModelPart, VtuOutput::CellType::Conditions, p_indices_map));
    }

    if (!has_elements && !has_conditions) {
        rOutput.push_back(std::make_tuple(&rModelPart, VtuOutput::CellType::None, p_indices_map));
    }

    if (has_elements || has_conditions) {
        IndexType vtu_index = 0;
        for (const auto& r_node : rModelPart.Nodes()) {
            (*p_indices_map)[r_node.Id()] = vtu_index++;
        }
    }
}

void FillCellDataRecursively(
    std::vector<std::tuple<ModelPart*, VtuOutput::CellType, std::unordered_map<IndexType, IndexType>*>>& rOutput,
    ModelPart& rModelPart)
{
    FillCellData(rOutput, rModelPart);
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        FillCellDataRecursively(rOutput, r_sub_model_part);
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
        FillCellDataRecursively(mModelPartCellData, rModelPart);
    } else {
        FillCellData(mModelPartCellData, rModelPart);
    }
}

VtuOutput::~VtuOutput()
{
    std::vector<std::unordered_map<IndexType, IndexType>*> list_of_maps;
    for ([[maybe_unused]] auto [p_model_part, cell_type, p_kratos_vtu_indices_map] : mModelPartCellData) {
        list_of_maps.push_back(p_kratos_vtu_indices_map);
    }

    std::sort(list_of_maps.begin(), list_of_maps.end());
    list_of_maps.erase(std::unique(list_of_maps.begin(), list_of_maps.end()), list_of_maps.end() );

    // now delete
    for (auto p_map : list_of_maps) {
        delete p_map;
    }
}

template<class TDataType>
void VtuOutput::AddHistoricalVariable(const Variable<TDataType>& rVariable)
{
    CheckDataArrayName(
        rVariable.Name(), mNonHistoricalNodalVariablesMap, mNodalFlagsMap, mPointFieldsMap);
    mHistoricalVariablesMap[rVariable.Name()] = &rVariable;
}

template<class TDataType>
void VtuOutput::AddNonHistoricalVariable(
    const Variable<TDataType>& rVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        CheckDataArrayName(
            rVariable.Name(), mHistoricalVariablesMap, mNodalFlagsMap,
            mPointFieldsMap);
        mNonHistoricalNodalVariablesMap[rVariable.Name()] = &rVariable;
    } else {
        CheckDataArrayName(
            rVariable.Name(), mCellFlagsMap, mCellFieldsMap);
        mNonHistoricalCellVariablesMap[rVariable.Name()] = &rVariable;
    }
}

template<class TDataType>
void VtuOutput::AddIntegrationPointVariable(
    const Variable<TDataType>& rVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        KRATOS_ERROR << "Nodes";
    } else {
        CheckDataArrayName(
            rVariable.Name(), mCellFlagsMap, mCellFieldsMap);
        mGaussPointCellVariablesMap[rVariable.Name()] = &rVariable;
    }
}

void VtuOutput::AddFlagVariable(
    const std::string& rFlagName,
    const Flags& rFlagVariable,
    const Flags& rEntityFlags)
{
    if (rEntityFlags.Is(NODES)) {
        CheckDataArrayName(
            rFlagName, mHistoricalVariablesMap, mNonHistoricalNodalVariablesMap,
            mPointFieldsMap);
        mNodalFlagsMap[rFlagName] = &rFlagVariable;
    } else {
        CheckDataArrayName(
            rFlagName, mNonHistoricalCellVariablesMap, mCellFieldsMap);
        mCellFlagsMap[rFlagName] = &rFlagVariable;
    }
}

template <class TContainerType>
void VtuOutput::AddContainerExpression(
    const std::string& rExpressionName,
    const typename ContainerExpression<TContainerType>::Pointer pContainerExpression)
{
    KRATOS_ERROR_IF_NOT(&pContainerExpression->GetModelPart() == &mrModelPart)
        << "Model part mismatch in container expression addition. [ Vtu output model part name = \""
        << mrModelPart.FullName() << "\", container expression model part name = \""
        << pContainerExpression->GetModelPart().FullName() << "\" ].\n";

    // generate the nd data
    const auto& r_expression = pContainerExpression->GetExpression();
    const auto& r_data_shape = r_expression.GetItemShape();
    const auto number_of_components = r_expression.GetItemComponentCount();

    DenseVector<unsigned int> shape(r_data_shape.size() + 1);
    std::copy(r_data_shape.begin(), r_data_shape.end(), shape.begin() + 1);
    shape[0] = r_expression.NumberOfEntities();

    auto p_nd_data = Kratos::make_shared<NDData<double>>(shape);
    const auto span = p_nd_data->ViewData();

    IndexPartition<IndexType>(r_expression.NumberOfEntities()).for_each([&span, &r_expression, number_of_components](const auto EntityIndex) {
        const auto data_begin_index = EntityIndex * number_of_components;
        for (IndexType i = 0; i < number_of_components; ++i) {
            span[data_begin_index + i] = r_expression.Evaluate(EntityIndex, data_begin_index, i);
        }
    });

    if constexpr (std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        CheckDataArrayName(
            rExpressionName, mHistoricalVariablesMap,
            mNonHistoricalNodalVariablesMap, mNodalFlagsMap);
        mPointFieldsMap[rExpressionName] = p_nd_data;
    } else {
        CheckDataArrayName(
            rExpressionName, mNonHistoricalCellVariablesMap, mCellFlagsMap);
        mCellFieldsMap[rExpressionName] = p_nd_data;
    }
}

template<class TCellsPointerType, class TXmlDataElementWrapper>
std::string VtuOutput::PrintModelPart(
    const std::string& rOutputFileNamePrefix,
    ModelPart::NodesContainerType::Pointer pNodes,
    TCellsPointerType pCells,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap,
    const DataCommunicator& rDataCommunicator,
    TXmlDataElementWrapper& rXmlDataElementWrapper) const
{
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
    piece_element->AddAttribute("NumberOfPoints", std::to_string(pNodes->size()));
    piece_element->AddAttribute("NumberOfCells", std::to_string(pCells->size()));
    unstructured_grid_element->AddElement(piece_element);

    // create the position element
    NodePositionTensorAdaptor node_position_tensor_adaptor(pNodes, mIsInitialConfiguration ? Globals::Configuration::Initial : Globals::Configuration::Current);
    node_position_tensor_adaptor.Check();
    node_position_tensor_adaptor.CollectData();

    // create the points element
    auto points_element = Kratos::make_shared<XmlElementsArray>("Points");
    points_element->AddElement(rXmlDataElementWrapper.Get("Position", node_position_tensor_adaptor.pGetStorage()));
    piece_element->AddElement(points_element);

    // create the cells element
    auto cells_element = CreateCellsXmlElement(pCells, rKratosVtuIndicesMap, rXmlDataElementWrapper);
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
    piece_element->AddElement(point_data_element);

    // generate and add point field data
    AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*point_data_element, pNodes, mNodalFlagsMap, rXmlDataElementWrapper);
    AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*point_data_element, pNodes, mNonHistoricalNodalVariablesMap, rXmlDataElementWrapper);
    AddFieldsFromTensorAdaptor<HistoricalVariableTensorAdaptor>(*point_data_element, pNodes, mHistoricalVariablesMap, rXmlDataElementWrapper);
    AddFields(*point_data_element, mPointFieldsMap, rXmlDataElementWrapper);

    // create cell data
    auto cell_data_element = Kratos::make_shared<XmlElementsArray>("CellData");
    piece_element->AddElement(cell_data_element);

    // generate and add cell field data
    AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*cell_data_element, pCells, mCellFlagsMap, rXmlDataElementWrapper);
    AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*cell_data_element, pCells, mNonHistoricalCellVariablesMap, rXmlDataElementWrapper);
    AddFields(*cell_data_element, mCellFieldsMap, rXmlDataElementWrapper);

    std::stringstream output_vtu_file_name;
    output_vtu_file_name << rOutputFileNamePrefix;
    if (rDataCommunicator.IsDistributed()) {
        output_vtu_file_name << "_" << rDataCommunicator.Rank();
    }
    output_vtu_file_name << ".vtu";

    std::ofstream output_file;
    output_file.open(output_vtu_file_name.str(), std::ios::out | std::ios::trunc);
    vtk_file_element.Write(output_file);
    output_file.close();

    return output_vtu_file_name.str();
}

template<class TCellsPointerType, class TXmlDataElementWrapper>
std::string VtuOutput::PrintGaussPointFields(
    const std::string& rOutputFileNamePrefix,
    TCellsPointerType pCells,
    const DataCommunicator& rDataCommunicator,
    ProcessInfo::Pointer pProcessInfo,
    TXmlDataElementWrapper& rXmlDataElementWrapper) const
{
    // create the vtk file
    XmlElementsArray vtk_file_element("VTKFile");
    vtk_file_element.AddAttribute("type", "UnstructuredGrid");
    vtk_file_element.AddAttribute("version", "0.1");
    vtk_file_element.AddAttribute("byte_order", GetEndianness());

    // create the unstructured grid
    auto unstructured_grid_element = Kratos::make_shared<XmlElementsArray>("UnstructuredGrid");
    vtk_file_element.AddElement(unstructured_grid_element);

    const auto total_gauss_points = block_for_each<SumReduction<IndexType>>(*pCells, [](const auto& rEntity) {
            return rEntity.GetGeometry().IntegrationPointsNumber(rEntity.GetIntegrationMethod());
        });

    // create the piece element
    auto piece_element = Kratos::make_shared<XmlElementsArray>("Piece");
    piece_element->AddAttribute("NumberOfPoints", std::to_string(total_gauss_points));
    piece_element->AddAttribute("NumberOfCells", "0");
    unstructured_grid_element->AddElement(piece_element);

    DenseVector<unsigned int> gauss_point_nd_data_shape(2);
    gauss_point_nd_data_shape[0] = total_gauss_points;
    gauss_point_nd_data_shape[1] = 3;
    auto gauss_point_positions = Kratos::make_shared<NDData<double>>(gauss_point_nd_data_shape);
    const auto span = gauss_point_positions->ViewData();
    auto gauss_point_itr = span.begin();
    array_1d<double, 3> coordinates;

    for (const auto& r_entity : (*pCells)) {
        const auto number_of_gauss_points = r_entity.GetGeometry().IntegrationPointsNumber(r_entity.GetIntegrationMethod());
        for (IndexType i = 0; i < number_of_gauss_points; ++i) {
            r_entity.GetGeometry().GlobalCoordinates(coordinates, i);
            *(gauss_point_itr++) = coordinates[0];
            *(gauss_point_itr++) = coordinates[1];
            *(gauss_point_itr++) = coordinates[2];
        }
    }

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
    for (const auto& r_pair : mGaussPointCellVariablesMap) {
        std::visit([&pCells, &point_data_element, &rXmlDataElementWrapper, &pProcessInfo, &is_gauss_point_data_available](auto pVariable) {
            using data_type = typename DataTypeTraits<typename BareType<decltype(*pVariable)>::Type>::PrimitiveType;
            if constexpr(std::is_same_v<data_type, double>) {
                GaussPointVariableTensorAdaptor gp_ta(pCells, pVariable,  pProcessInfo);
                gp_ta.CollectData();
                const auto& gp_ta_shape = gp_ta.Shape();

                KRATOS_ERROR_IF(gp_ta_shape.size() < 2)
                    << "Gauss points requires at least 2 dimensional tensor data";

                DenseVector<unsigned int> shape(gp_ta_shape.size() - 1);
                std::copy(gp_ta_shape.begin() + 2, gp_ta_shape.end(), shape.begin());
                shape[0] = gp_ta_shape[0] * gp_ta_shape[1];

                if (gp_ta.Size() > 0) {
                    is_gauss_point_data_available = true;
                    point_data_element->AddElement(rXmlDataElementWrapper.Get(pVariable->Name(), gp_ta.pGetStorage(), shape));
                }
            }
        }, r_pair.second);
    }

    if (rDataCommunicator.OrReduceAll(is_gauss_point_data_available)) {
        std::stringstream output_vtu_file_name;
        output_vtu_file_name << rOutputFileNamePrefix << "_GAUSS";
        if (rDataCommunicator.IsDistributed()) {
            output_vtu_file_name << "_" << rDataCommunicator.Rank();
        }
        output_vtu_file_name << ".vtu";

        std::ofstream output_file;
        output_file.open(output_vtu_file_name.str(), std::ios::out | std::ios::trunc);
        vtk_file_element.Write(output_file);
        output_file.close();

        return output_vtu_file_name.str();
    } else {
        return "";
    }
}

template<class TCellsPointerType>
std::string VtuOutput::PrintModelPartDispatcher(
    const std::string& rOutputFileNamePrefix,
    ModelPart::NodesContainerType::Pointer pNodes,
    TCellsPointerType pCells,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap,
    const DataCommunicator& rDataCommunicator) const
{
    switch (mOutputFormat)  {
        case ASCII: {
                XmlAsciiNDDataElementWrapper data_element_wrapper{mPrecision};
                return PrintModelPart(rOutputFileNamePrefix, pNodes, pCells, rKratosVtuIndicesMap, rDataCommunicator, data_element_wrapper);
            }
        case BINARY: {
                XmlBase64BinaryNDDataElementWrapper data_element_wrapper;
                return PrintModelPart(rOutputFileNamePrefix, pNodes, pCells, rKratosVtuIndicesMap, rDataCommunicator, data_element_wrapper);
            }
    }
}

template<class TCellsPointerType>
std::string VtuOutput::PrintGaussPointFieldsDispatcher(
    const std::string& rOutputFileNamePrefix,
    TCellsPointerType pCells,
    const DataCommunicator& rDataCommunicator,
    ProcessInfo::Pointer pProcessInfo) const
{
    switch (mOutputFormat)  {
        case ASCII: {
                XmlAsciiNDDataElementWrapper data_element_wrapper{mPrecision};
                return PrintGaussPointFields(rOutputFileNamePrefix, pCells, rDataCommunicator, pProcessInfo, data_element_wrapper);
            }
        case BINARY: {
                XmlBase64BinaryNDDataElementWrapper data_element_wrapper;
                return PrintGaussPointFields(rOutputFileNamePrefix, pCells, rDataCommunicator, pProcessInfo, data_element_wrapper);
            }
    }
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
    mPointFieldsMap.clear();
}

void VtuOutput::ClearCellContainerExpressions()
{
    mCellFieldsMap.clear();
}

const ModelPart& VtuOutput::GetModelPart() const
{
    return mrModelPart;
}

void VtuOutput::PrintOutput(const std::string& rOutputFolderName)
{
    std::filesystem::create_directories(rOutputFolderName + "/" + mrModelPart.FullName());

    const auto& r_process_info = mrModelPart.GetProcessInfo();

    const std::string suffix = (r_process_info.Has(STEP) ? "_" + std::to_string(r_process_info[STEP]) : "");

    std::vector<std::pair<std::string, std::string>> vtm_info;
    for ([[maybe_unused]] auto [p_model_part, cell_type, p_kratos_vtu_indices_map] : mModelPartCellData) {
        const auto& r_data_communicator = p_model_part->GetCommunicator().GetDataCommunicator();
        auto p_process_info = p_model_part->pGetProcessInfo();

        std::string cell_type_suffix;
        switch (cell_type) {
            case Elements:
                cell_type_suffix = "_ELEMENTS";
                break;
            case Conditions:
                cell_type_suffix = "_CONDITIONS";
                break;
            default:
                cell_type_suffix = "_NODES";
                break;
        }

        const auto& current_mp_name = rOutputFolderName + "/" + mrModelPart.FullName() + "/" + p_model_part->FullName() + cell_type_suffix + suffix;
        const auto& r_model_part_name = p_model_part->FullName() + cell_type_suffix;
        const auto& r_kratos_vtu_indices_map = *p_kratos_vtu_indices_map;
        switch (cell_type) {
            case Elements:
                vtm_info.push_back(std::make_pair(r_model_part_name, PrintModelPartDispatcher(current_mp_name, p_model_part->pNodes(), p_model_part->pElements(), r_kratos_vtu_indices_map, r_data_communicator)));
                if (!mGaussPointCellVariablesMap.empty()) {
                    vtm_info.push_back(std::make_pair(r_model_part_name + "_GAUSS", PrintGaussPointFieldsDispatcher(current_mp_name, p_model_part->pElements(), r_data_communicator, p_process_info)));
                }
                break;
            case Conditions:
                vtm_info.push_back(std::make_pair(r_model_part_name, PrintModelPartDispatcher(current_mp_name, p_model_part->pNodes(), p_model_part->pConditions(), r_kratos_vtu_indices_map, r_data_communicator)));
                if (!mGaussPointCellVariablesMap.empty()) {
                    vtm_info.push_back(std::make_pair(r_model_part_name + "_GAUSS", PrintGaussPointFieldsDispatcher(current_mp_name, p_model_part->pConditions(), r_data_communicator, p_process_info)));
                }
                break;
            default:
                vtm_info.push_back(std::make_pair(r_model_part_name, PrintModelPartDispatcher(current_mp_name, p_model_part->pNodes(), p_model_part->pElements(), r_kratos_vtu_indices_map, r_data_communicator)));
                break;
        }
    }

    if (mrModelPart.GetCommunicator().MyPID() == 0) {
        // now generate the vtm file
        XmlElementsArray vtm_file_element("VTKFile");
        vtm_file_element.AddAttribute("type", "vtkMultiBlockDataSet");
        vtm_file_element.AddAttribute("version", "1.0");
        vtm_file_element.AddAttribute("byte_order", GetEndianness());

        auto multi_block_set_element = Kratos::make_shared<XmlElementsArray>("vtkMultiBlockDataSet");
        vtm_file_element.AddElement(multi_block_set_element);

        IndexType local_index = 0;
        for (IndexType i = 0; i < vtm_info.size(); ++i) {
            if (vtm_info[i].second != "") {
                auto current_element = Kratos::make_shared<XmlElementsArray>("DataSet");
                current_element->AddAttribute("index", std::to_string(local_index++));
                current_element->AddAttribute("name", vtm_info[i].first);
                current_element->AddAttribute("file", vtm_info[i].second.substr(rOutputFolderName.size() + 1));
                multi_block_set_element->AddElement(current_element);
            }
        }

        std::ofstream output_file;
        output_file.open(rOutputFolderName + "/" + mrModelPart.FullName() + suffix + ".vtm", std::ios::out | std::ios::trunc);
        vtm_file_element.Write(output_file);
        output_file.close();

        // now generate the pvd file
        if (!suffix.empty()) {
            // get the step
            const auto step = r_process_info[STEP];

            // only generate if the STEP variable is found
            const double time = (r_process_info.Has(TIME) ? r_process_info[TIME] : static_cast<double>(step));

            mStepInfo.push_back(std::make_pair(step, time));

            XmlElementsArray pvd_file_element("VTKFile");
            pvd_file_element.AddAttribute("type", "Collection");
            pvd_file_element.AddAttribute("version", "0.1");

            auto collection_element = Kratos::make_shared<XmlElementsArray>("Collection");
            pvd_file_element.AddElement(collection_element);

            for ([[maybe_unused]] auto [current_step, current_time] : mStepInfo) {
                std::stringstream str_time;
                str_time << std::scientific << std::setprecision(mPrecision) << current_time;

                auto p_step_element = Kratos::make_shared<XmlElementsArray>("DataSet");
                p_step_element->AddAttribute("timestep", str_time.str());
                p_step_element->AddAttribute("file", mrModelPart.FullName() + "_" + std::to_string(current_step) + ".vtm");
                collection_element->AddElement(p_step_element);
            }

            std::ofstream output_file;
            output_file.open(rOutputFolderName + "/" + mrModelPart.FullName() + ".pvd", std::ios::out | std::ios::trunc);
            pvd_file_element.Write(output_file);
            output_file.close();
        }
    }
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

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddIntegrationPointVariable(const Variable<int>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddIntegrationPointVariable(const Variable<double>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddIntegrationPointVariable(const Variable<array_1d<double, 3>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddIntegrationPointVariable(const Variable<array_1d<double, 4>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddIntegrationPointVariable(const Variable<array_1d<double, 6>>&, const Flags&);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddIntegrationPointVariable(const Variable<array_1d<double, 9>>&, const Flags&);

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::NodesContainerType>(const std::string&, const typename ContainerExpression<ModelPart::NodesContainerType>::Pointer);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ConditionsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ConditionsContainerType>::Pointer);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ElementsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ElementsContainerType>::Pointer);

} // namespace Kratos
