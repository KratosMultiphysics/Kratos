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
};

struct XmlBase64BinaryNDDataElementWrapper
{
    IndexType mPrecision;

    template<class NDDataTypePointer>
    XmlElement::Pointer Get(
        const std::string& rDataArrayName,
        NDDataTypePointer pNDData)
    {
        using data_type = typename BareType<decltype(*pNDData)>::DataType;
        return Kratos::make_shared<XmlBase64BinaryNDDataElement<data_type>>(rDataArrayName, pNDData);
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

template<class TContainerType, class TXmlDataElementWrapper>
XmlElement::Pointer CreateCellsXmlElement(
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap,
    TXmlDataElementWrapper& rXmlDataElementWrapper)
{
    auto p_cells_xml_element = Kratos::make_shared<XmlElementsArray>("Cells");

    auto p_offsets = GetOffsets(rContainer);

    p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("connectivity", GetConnectivities(*p_offsets, rContainer, rKratosVtuIndicesMap)));
    p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("offsets", p_offsets));
    p_cells_xml_element->AddElement(rXmlDataElementWrapper.Get("types", GetGeometryTypes(rContainer)));

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
            tensor_adaptor.Check();
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
                    tensor_adaptor.Check();
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

} // namespace

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
        KRATOS_ERROR_IF_NOT(!rEntityFlags.Is(CONDITIONS) || (rEntityFlags.Is(CONDITIONS) && mIsConditionsConsidered))
            << "Condition variable \""
            << rVariable.Name() << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";

        KRATOS_ERROR_IF_NOT(!rEntityFlags.Is(ELEMENTS) || (rEntityFlags.Is(ELEMENTS) && mIsElementsConsidered))
            << "Element variable \""
            << rVariable.Name() << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";

        CheckDataArrayName(
            rVariable.Name(), mCellFlagsMap, mCellFieldsMap);
        mNonHistoricalCellVariablesMap[rVariable.Name()] = &rVariable;
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
        KRATOS_ERROR_IF_NOT(!rEntityFlags.Is(CONDITIONS) || (rEntityFlags.Is(CONDITIONS) && mIsConditionsConsidered))
            << "Condition flag \""
            << rFlagName << "\" cannot be written for a model part with elements [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";

        KRATOS_ERROR_IF_NOT(!rEntityFlags.Is(ELEMENTS) || (rEntityFlags.Is(ELEMENTS) && mIsElementsConsidered))
            << "Element flag \""
            << rFlagName << "\" cannot be written for a model part with only conditions [ model part name = \""
            << mrModelPart.FullName() << "\" ].\n";

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

template<class TXmlDataElementWrapper>
void VtuOutput::PrintModelPart(
    const std::string& rOutputFileNamePrefix,
    ModelPart& rModelPart,
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
    points_element->AddElement(rXmlDataElementWrapper.Get("Position", node_position_tensor_adaptor.pGetStorage()));
    piece_element->AddElement(points_element);

    // create the cells element
    auto cells_element = mIsElementsConsidered
                            ? CreateCellsXmlElement(rModelPart.Elements(), mKratosVtuIndicesMap, rXmlDataElementWrapper)
                            : CreateCellsXmlElement(rModelPart.Conditions(), mKratosVtuIndicesMap, rXmlDataElementWrapper);
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<XmlElementsArray>("PointData");
    piece_element->AddElement(point_data_element);

    // generate and add point field data
    AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*point_data_element, rModelPart.pNodes(), mNodalFlagsMap, rXmlDataElementWrapper);
    AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*point_data_element, rModelPart.pNodes(), mNonHistoricalNodalVariablesMap, rXmlDataElementWrapper);
    AddFieldsFromTensorAdaptor<HistoricalVariableTensorAdaptor>(*point_data_element, rModelPart.pNodes(), mHistoricalVariablesMap, rXmlDataElementWrapper);
    AddFields(*point_data_element, mPointFieldsMap, rXmlDataElementWrapper);

    // create cell data
    auto cell_data_element = Kratos::make_shared<XmlElementsArray>("CellData");
    piece_element->AddElement(cell_data_element);

    // generate and add cell field data
    if (mIsElementsConsidered) {
        AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*cell_data_element, rModelPart.pElements(), mCellFlagsMap, rXmlDataElementWrapper);
        AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*cell_data_element, rModelPart.pElements(), mNonHistoricalCellVariablesMap, rXmlDataElementWrapper);
    } else if (mIsConditionsConsidered) {
        AddFieldsFromTensorAdaptor<FlagsTensorAdaptor>(*cell_data_element, rModelPart.pConditions(), mCellFlagsMap, rXmlDataElementWrapper);
        AddFieldsFromTensorAdaptor<VariableTensorAdaptor>(*cell_data_element, rModelPart.pConditions(), mNonHistoricalCellVariablesMap, rXmlDataElementWrapper);
    }
    AddFields(*cell_data_element, mCellFieldsMap, rXmlDataElementWrapper);

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

void VtuOutput::PrintOutput(const std::string& rOutputFilenamePrefix)
{
    switch (mOutputFormat) {
        case ASCII:
            {
                XmlAsciiNDDataElementWrapper data_element_wrapper{mPrecision};
                PrintModelPart(rOutputFilenamePrefix, mrModelPart, data_element_wrapper);
                break;
            }
        case BINARY:
            {
                XmlBase64BinaryNDDataElementWrapper data_element_wrapper;
                PrintModelPart(rOutputFilenamePrefix, mrModelPart, data_element_wrapper);
                break;
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

template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::NodesContainerType>(const std::string&, const typename ContainerExpression<ModelPart::NodesContainerType>::Pointer);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ConditionsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ConditionsContainerType>::Pointer);
template KRATOS_API(KRATOS_CORE) void VtuOutput::AddContainerExpression<ModelPart::ElementsContainerType>(const std::string&, const typename ContainerExpression<ModelPart::ElementsContainerType>::Pointer);

} // namespace Kratos
