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
#include <string>
#include <vector>
#include <iomanip>
#include <numeric>

// External includes

// Project includes
#include "containers/container_expression/expressions/expression.h"
#include "containers/container_expression/expressions/expression_iterator.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "containers/container_expression/specialized_container_expression.h"
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/io.h"
#include "utilities/global_pointer_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/reduction_utilities.h"
#include "utilities/string_utilities.h"

// Include base h
#include "vtu_output.h"

namespace Kratos {

    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, NODES,   1);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, CONDITIONS,  2);
    KRATOS_CREATE_LOCAL_FLAG(VtuOutput, ELEMENTS, 3);

namespace VtuOutputHelperUtilities {

class KRATOS_API(KRATOS_CORE) XmlOStreamWriter
{
public:
    ///@name Life cycle
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    XmlOStreamWriter(
        std::ostream& rOStream,
        const VtuOutput::WriterFormat OutputFormat,
        const IndexType Precision)
        : mrOStream(rOStream),
        mOutputFormat(OutputFormat)
    {
        mrOStream << std::scientific << std::setprecision(Precision);
    }

    ///@}
    ///@name Public operations
    ///@{

    void WriteElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level,
        const bool IsEmptyElement)
    {
        WriteAttributes(rTagName, rAttributes, Level);

        if (IsEmptyElement) {
            mrOStream << "/>\n";
        } else {
            mrOStream << ">\n";
        }

    }

    void CloseElement(
        const std::string& rTagName,
        const IndexType Level)
    {
        const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
        mrOStream << tabbing << "</" << rTagName << ">\n";
    }

    void WriteDataElement(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const std::vector<Expression::Pointer>& rExpressions,
        const IndexType Level)
    {

        switch (mOutputFormat){
            case VtuOutput::WriterFormat::ASCII:
                if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<char>*>(&*pExpression); })) {
                    WriteDataElementAscii<LiteralFlatExpression<char>>(rTagName, rAttributes, Level, rExpressions);
                } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<int>*>(&*pExpression); })) {
                    WriteDataElementAscii<LiteralFlatExpression<int>>(rTagName, rAttributes, Level, rExpressions);
                } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<double>*>(&*pExpression); })) {
                    WriteDataElementAscii<LiteralFlatExpression<double>>(rTagName, rAttributes, Level, rExpressions);
                } else {
                    WriteDataElementAscii<Expression>(rTagName, rAttributes, Level, rExpressions);
                }
                break;
            case VtuOutput::WriterFormat::BINARY:
                if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<char>*>(&*pExpression); })) {
                    WriteDataElementBinary<LiteralFlatExpression<char>>(rTagName, rAttributes, Level, rExpressions);
                } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<int>*>(&*pExpression); })) {
                    WriteDataElementBinary<LiteralFlatExpression<int>>(rTagName, rAttributes, Level, rExpressions);
                } else if (std::all_of(rExpressions.begin(), rExpressions.end(), [](const auto& pExpression) { return dynamic_cast<LiteralFlatExpression<double>*>(&*pExpression); })) {
                    WriteDataElementBinary<LiteralFlatExpression<double>>(rTagName, rAttributes, Level, rExpressions);
                } else {
                    WriteDataElementBinary<Expression>(rTagName, rAttributes, Level, rExpressions);
                }
                break;
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::ostream& mrOStream;

    const VtuOutput::WriterFormat mOutputFormat;

    ///@}
    ///@name Private operations
    ///@{

    void WriteAttributes(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level)
    {
        const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
        mrOStream << tabbing << "<" << rTagName;
        if (rAttributes.size() > 0) {
            for (const auto& r_pair : rAttributes) {
                mrOStream << " " << r_pair.first << "=\"" << r_pair.second << "\"";
            }
        }
    }

    template<class TExpressionType>
    void WriteDataElementAscii(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level,
        const std::vector<Expression::Pointer>& rExpressions)
    {
        using exp_itr_type = ExpressionIterator<TExpressionType>;

        using data_itr_type = typename exp_itr_type::ConstIteratorType;

        WriteAttributes(rTagName, rAttributes, Level);
        // add format
        const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
        mrOStream << " format=\"ascii\">\n" << tabbing;

        std::vector<TExpressionType*> transformed_expressions(rExpressions.size());
        std::transform(rExpressions.begin(), rExpressions.end(),
                    transformed_expressions.begin(), [](auto& pExpression) {
                        return dynamic_cast<TExpressionType*>(&*(pExpression));
                    });

        for (const auto& p_expression : transformed_expressions) {
            exp_itr_type expression_iterator(p_expression);
            auto data_begin = expression_iterator.cbegin();
            auto data_end   = expression_iterator.cend();
            for (data_itr_type itr = data_begin; itr != data_end; ++itr) {
                if constexpr(std::is_same_v<typename exp_itr_type::value_type, char>) {
                    mrOStream << "  " << static_cast<int>(*(itr));
                } else {
                    mrOStream << "  " << *itr;
                }
            }
        }

        mrOStream << "\n";
        CloseElement(rTagName, Level);
    }

    template<class TExpressionType>
    void WriteDataElementBinary(
        const std::string& rTagName,
        const std::vector<std::pair<const std::string, const std::string>>& rAttributes,
        const IndexType Level,
        const std::vector<Expression::Pointer>& rExpressions)
    {
        using exp_itr_type = ExpressionIterator<TExpressionType>;

        using value_type = typename exp_itr_type::value_type;

        using data_itr_type = typename exp_itr_type::ConstIteratorType;

        WriteAttributes(rTagName, rAttributes, Level);

        std::vector<TExpressionType*> transformed_expressions(rExpressions.size());
        std::transform(rExpressions.begin(), rExpressions.end(),
                    transformed_expressions.begin(), [](auto& pExpression) {
                        return dynamic_cast<TExpressionType*>(&*(pExpression));
                    });

        if (rExpressions.size() == 0) {
            mrOStream << " format=\"binary\"/>\n";
            return;
        } else {
            const std::string& tabbing = XmlOStreamWriter::GetTabbing(Level);
            mrOStream << " format=\"binary\">\n" << tabbing << "  ";
        }

        constexpr char base64_map[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

        const IndexType total_entities = std::accumulate(rExpressions.begin(), rExpressions.end(), 0U, [](const IndexType LHS, const auto& pExpression) { return LHS + pExpression->NumberOfEntities();});

        using writing_data_type = value_type;
        constexpr IndexType data_type_size = sizeof(writing_data_type);
        constexpr IndexType size_type_size = sizeof(unsigned int);
        const IndexType total_number_of_values = total_entities * rExpressions[0]->GetFlattenedShapeSize();
        const IndexType size_type_triplets = size_type_size / 3;
        const unsigned int total_data_size = total_number_of_values * data_type_size;
        const IndexType total_number_of_triplets = (total_data_size + size_type_size)  / 3;

        IndexType byte_index = 0;
        auto p_expression = transformed_expressions.data();
        data_itr_type data_itr = exp_itr_type(*p_expression).cbegin();
        data_itr_type data_end = exp_itr_type(*p_expression).cend();
        writing_data_type current_value =  *data_itr;

        auto get_next_byte = [&]() -> char {
            if (byte_index == data_type_size) {
                byte_index = 0;
                ++data_itr;

                if (data_itr == data_end) {
                    ++p_expression;
                    data_itr = exp_itr_type(*p_expression).cbegin();
                    data_end = exp_itr_type(*p_expression).cend();
                }

                current_value = *data_itr;
            }

            const char byte = *(reinterpret_cast<const char*>(&current_value) + byte_index++);

            return byte;
        };

        auto write_encoded_triplet = [&](const std::array<char, 3>& bytes, size_t padding) {
            char tmp[5] = {
                base64_map[(bytes[0] & 0xfc) >> 2],
                base64_map[((bytes[0] & 0x03) << 4) + ((bytes[1] & 0xf0) >> 4)],
                base64_map[((bytes[1] & 0x0f) << 2) + ((bytes[2] & 0xc0) >> 6)],
                base64_map[bytes[2] & 0x3f], '\0'};

            std::fill(tmp + 4 - padding, tmp + 4, '=');

            mrOStream << tmp;
        };

        // first write the total number of bytes in the array this will be of 8 bytes
        auto total_data_size_begin = reinterpret_cast<const char*>(&total_data_size);
        IndexType local_index = 0;
        IndexType initial_number_of_triplets = 0;
        for (IndexType i = 0; i < size_type_triplets; ++i) {
            write_encoded_triplet({*(total_data_size_begin + local_index++),
                                   *(total_data_size_begin + local_index++),
                                   *(total_data_size_begin + local_index++)},
                                  0);
            ++initial_number_of_triplets;
        }

        std::array<char, 3> data;
        switch (size_type_size % 3) {
            case 1:
                data[0] = *(total_data_size_begin + local_index++);
                data[1] = get_next_byte();
                data[2] = get_next_byte();
                ++initial_number_of_triplets;
                write_encoded_triplet(data, 0);
                break;
            case 2:
                data[0] = *(total_data_size_begin + local_index++);
                data[1] = *(total_data_size_begin + local_index++);
                data[2] = get_next_byte();
                ++initial_number_of_triplets;
                write_encoded_triplet(data, 0);
                break;
            default:
                break;
        }

        for (IndexType i = initial_number_of_triplets; i < total_number_of_triplets; ++i) {
            write_encoded_triplet({get_next_byte(), get_next_byte(), get_next_byte()}, 0);
        }

        const IndexType number_of_bytes_remaining = (total_data_size + size_type_size) % 3;
        if (number_of_bytes_remaining != 0) {
            std::array<char, 3> bytes{'\0', '\0', '\0'};

            for (IndexType i = 0; i < number_of_bytes_remaining; ++i) {
                bytes[i] = get_next_byte();
            }
            write_encoded_triplet(bytes, 3 - number_of_bytes_remaining);
        }

        mrOStream << "\n";
        CloseElement(rTagName, Level);
    }

    static const std::string GetTabbing(const IndexType Level)
    {
        std::stringstream ss_tabbing;
        for (IndexType i = 0; i < Level; ++i) {
            ss_tabbing << "   ";
        }
        return ss_tabbing.str();
    }

    ///@}
};

class XmlElement {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(XmlElement);

    ///@}
    ///@name Life cycle
    ///@{

    XmlElement(const std::string& rTagName)
        : mTagName(rTagName)
    {
    }

    XmlElement(
        const std::string& rDataName,
        const std::vector<Expression::Pointer>& rExpressions)
        : mTagName("DataArray"),
        mExpressions(rExpressions)
    {
        KRATOS_ERROR_IF(rExpressions.size() == 0)
            << "Empty expression lists are not allowed.";

        IndexType number_of_components = 0;

        for (const auto& p_expression : rExpressions) {
            if (number_of_components == 0) {
                number_of_components = p_expression->GetFlattenedShapeSize();
            }
            KRATOS_ERROR_IF(number_of_components != p_expression->GetFlattenedShapeSize())
                << "Found expressions with mismatching shapes.";
        }

        if (std::all_of(mExpressions.begin(), mExpressions.end(), [](const auto& pExpression) {
                return dynamic_cast<LiteralFlatExpression<char>*>(&*pExpression);
            })) {
            AddAttribute("type", "UInt8");
        } else if (std::all_of(mExpressions.begin(), mExpressions.end(), [](const auto& pExpression) {
                return dynamic_cast<LiteralFlatExpression<int>*>(&*pExpression);
            })) {
            AddAttribute("type", "Int32");
        } else {
            AddAttribute("type", "Float64");
        }

        AddAttribute("Name", rDataName);

        if (number_of_components > 0) {
            AddAttribute("NumberOfComponents", std::to_string(number_of_components));
        }
    }

    ///@}
    ///@name Public operations
    ///@{

    const std::string GetTagName() const { return mTagName; };

    void AddAttribute(
        const std::string& rName,
        const std::string& rValue)
    {
        for (const auto& r_attribute : mAttributes) {
            KRATOS_ERROR_IF(r_attribute.first == rName)
                << "There exists an attribute named \"" << rName
                << "\" in the xml element with value = \""
                << r_attribute.second << "\" [ given new value = \""
                << rValue << "\" ].\n";
        }
        mAttributes.push_back(std::make_pair(rName, rValue));
    }

    const std::vector<std::pair<const std::string, const std::string>>& GetAttributes() const { return mAttributes; }

    void ClearAttributes() { mAttributes.clear(); }

    void AddElement(const XmlElement::Pointer pXmlElement)
    {
        if (mExpressions.size() == 0) {
            for (const auto& p_element : mXmlElements) {
                KRATOS_ERROR_IF(&*(p_element) == &*(pXmlElement))
                    << "The xml element is already aded.";
            }
            mXmlElements.push_back(pXmlElement);
        } else {
            KRATOS_ERROR << "Cannot add element to an Xml element which has "
                            "data [ current xml tag = \""
                        << GetTagName() << "\", new element tag name = \""
                        << pXmlElement->GetTagName() << "\" ].\n";
        }
    }

    std::vector<XmlElement::Pointer> GetElements(const std::string& rTagName) const
    {
        std::vector<XmlElement::Pointer> results;
        for (const auto& p_element : mXmlElements) {
            if (p_element->GetTagName() == rTagName) {
                results.push_back(p_element);
            }
        }
        return results;
    }

    const std::vector<XmlElement::Pointer>& GetElements() const { return mXmlElements; }

    void ClearElements() { mXmlElements.clear(); }

    void Write(
        XmlOStreamWriter& rWriter,
        const IndexType Level = 0) const
    {
        if (mXmlElements.size() > 0) {
            rWriter.WriteElement(GetTagName(), GetAttributes(), Level, false);
            for (const auto& p_element : mXmlElements) {
                p_element->Write(rWriter, Level + 1);
            }
            rWriter.CloseElement(GetTagName(), Level);
        } else if (mExpressions.size() > 0) {
            rWriter.WriteDataElement(GetTagName(), GetAttributes(), mExpressions, Level);
        } else {
            rWriter.WriteElement(GetTagName(), GetAttributes(), Level, true);
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::string mTagName;

    std::vector<std::pair<const std::string, const std::string>> mAttributes;

    std::vector<XmlElement::Pointer> mXmlElements;

    const std::vector<Expression::Pointer> mExpressions;

    ///@}
};

const std::map<GeometryData::KratosGeometryType, char> KratosVtuGeometryTypes = {
    { GeometryData::KratosGeometryType::Kratos_Point2D,          1 },
    { GeometryData::KratosGeometryType::Kratos_Point3D,          1 },
    { GeometryData::KratosGeometryType::Kratos_Line2D2,          3 },
    { GeometryData::KratosGeometryType::Kratos_Line3D2,          3 },
    { GeometryData::KratosGeometryType::Kratos_Triangle2D3,      5 },
    { GeometryData::KratosGeometryType::Kratos_Triangle3D3,      5 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, 9 },
    { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    10 },
    { GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     12 },
    { GeometryData::KratosGeometryType::Kratos_Prism3D6,         13 },
    { GeometryData::KratosGeometryType::Kratos_Pyramid3D5,       14 },
    { GeometryData::KratosGeometryType::Kratos_Line2D3,          21 },
    { GeometryData::KratosGeometryType::Kratos_Line3D3,          21 },
    { GeometryData::KratosGeometryType::Kratos_Triangle2D6,      22 },
    { GeometryData::KratosGeometryType::Kratos_Triangle3D6,      22 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, 23 },
    { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, 23 },
    { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,   24 },
    { GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,    25 },
    { GeometryData::KratosGeometryType::Kratos_Prism3D15,        26 },
    { GeometryData::KratosGeometryType::Kratos_Pyramid3D13,      27 }
};

Expression::Pointer CreatePositionsExpression(
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

XmlElement::Pointer CreateDataArrayElement(
    const std::string& rDataArrayName,
    const std::vector<Expression*>& rExpressions)
{
    std::vector<Expression::Pointer> expressions;
    for (const auto& p_expression : rExpressions) {
        if (p_expression) {
            expressions.push_back(p_expression);
        }
    }
    if (expressions.size() > 0) {
        return Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>(rDataArrayName, expressions);
    } else {
        return nullptr;
    }
}

template<class TContainerType, class TMeshType>
Expression* pGetExpression(const ContainerExpression<TContainerType, TMeshType>& rContainerExpression)
{
    if (rContainerExpression.HasExpression()) {
        return &*rContainerExpression.pGetExpression();
    } else {
        return nullptr;
    }
}

void AddDataArrayElement(
    VtuOutputHelperUtilities::XmlElement::Pointer pParentElement,
    VtuOutputHelperUtilities::XmlElement::Pointer pDataArrayElement)
{
    if (pDataArrayElement) {
        pParentElement->AddElement(pDataArrayElement);
    }
}

VtuOutputHelperUtilities::XmlElement::Pointer CreatePointsXmlElement(
    const ModelPart& rModelPart,
    const bool IsInitialConfiguration)
{
    auto p_points_xml_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("Points");

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
        const auto p_itr = KratosVtuGeometryTypes.find((rContainer.begin() + Index)->GetGeometry().GetGeometryType());
        if (p_itr != KratosVtuGeometryTypes.end()) {
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
VtuOutputHelperUtilities::XmlElement::Pointer CreateCellsXmlElement(
    const TContainerType& rContainer,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap)
{
    auto p_cells_xml_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("Cells");

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
VtuOutputHelperUtilities::XmlElement::Pointer CreateVariableDataXmlElement(
    const Variable<TDataType>& rVariable,
    ModelPart& rModelPart)
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        SpecializedContainerExpression<TContainerType, ContainerDataIO<TContainerDataIOTag>, MeshType::Local> local_container(rModelPart);
        local_container.Read(rVariable);
        SpecializedContainerExpression<TContainerType, ContainerDataIO<TContainerDataIOTag>, MeshType::Ghost> ghost_container(rModelPart);
        ghost_container.Read(rVariable);

        return CreateDataArrayElement(
            rVariable.Name(),
            {pGetExpression(local_container), pGetExpression(ghost_container)});
    } else {
        SpecializedContainerExpression<TContainerType, ContainerDataIO<TContainerDataIOTag>> local_container(rModelPart);
        local_container.Read(rVariable);

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
VtuOutputHelperUtilities::XmlElement::Pointer CreateFlagDataXmlElement(
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

    const IndexType number_of_components = rLocalNodesExpression.GetFlattenedShapeSize();

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

    auto p_ghost_nodes_expression = LiteralFlatExpression<double>::Create(number_of_ghost_nodes, rLocalNodesExpression.GetShape());
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
VtuOutputHelperUtilities::XmlElement::Pointer CreateContainerExpressionXmlElement(
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
    VtuOutputHelperUtilities::XmlElement::Pointer pXmlDataElement,
    ModelPart& rModelPart,
    const std::unordered_map<std::string, VtuOutput::SupportedVariables>& rVariablesMap)
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
    VtuOutputHelperUtilities::XmlElement::Pointer pXmlDataElement,
    ModelPart& rModelPart,
    const std::unordered_map<std::string, const Flags*>& rFlagsMap)
{
    for (const auto& r_flag_pair : rFlagsMap){
        AddDataArrayElement(pXmlDataElement, CreateFlagDataXmlElement<TContaierType>(r_flag_pair.first, *r_flag_pair.second, rModelPart));
    }
}

template<class TContainerExpressionType>
void AddListOfContainerExpressions(
    VtuOutputHelperUtilities::XmlElement::Pointer pXmlDataElement,
    ModelPart& rModelPart,
    const std::unordered_map<IndexType, IndexType>& rKratosVtuIndicesMap,
    const std::unordered_map<std::string, TContainerExpressionType>& rContainerExpressionsMap)
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
    VtuOutputHelperUtilities::XmlElement& rOutputElement,
    const VtuOutputHelperUtilities::XmlElement& rInputElement)
{
    for (const auto& r_attribute_data : rInputElement.GetAttributes()) {
        rOutputElement.AddAttribute(r_attribute_data.first, r_attribute_data.second);
    }
}

void CreatePDataArrays(
    VtuOutputHelperUtilities::XmlElement& rOutputElement,
    const VtuOutputHelperUtilities::XmlElement& rInputElement)
{
    for (const auto& p_element : rInputElement.GetElements()) {
        auto new_pdata_array = Kratos::make_shared<XmlElement>("PDataArray");
        CopyAttributes(*new_pdata_array, *p_element);
        rOutputElement.AddElement(new_pdata_array);
    }
}

std::string GetEndianess()
{
    int i = 0x0001;

    if (*reinterpret_cast<char*>(&i) != 0) {
        return "LittleEndian";
    } else {
        return "BigEndian";
    }
}

void WritePvtuFile(
    VtuOutputHelperUtilities::XmlElement::Pointer pVtkFileElement,
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
        auto vtk_file_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("VTKFile");
        vtk_file_element->AddAttribute("type", "PUnstructuredGrid");
        vtk_file_element->AddAttribute("version", "0.1");
        vtk_file_element->AddAttribute("byte_order", GetEndianess());

        // create the unstructured grid
        auto unstructured_grid_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("PUnstructuredGrid");
        unstructured_grid_element->AddAttribute("GhostLevel", "0");
        vtk_file_element->AddElement(unstructured_grid_element);

        // get the points xml element
        auto piece = pVtkFileElement->GetElements("UnstructuredGrid")[0]->GetElements("Piece")[0];

        auto points = piece->GetElements("Points")[0];
        auto p_points = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("PPoints");
        VtuOutputHelperUtilities::CreatePDataArrays(*p_points, *points);
        unstructured_grid_element->AddElement(p_points);

        auto cells = piece->GetElements("Cells")[0];
        auto p_cells = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("PCells");
        VtuOutputHelperUtilities::CreatePDataArrays(*p_cells, *cells);
        unstructured_grid_element->AddElement(p_cells);

        auto point_data = piece->GetElements("PointData")[0];
        auto p_point_data = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("PPointData");
        VtuOutputHelperUtilities::CreatePDataArrays(*p_point_data, *point_data);
        unstructured_grid_element->AddElement(p_point_data);

        auto cell_data = piece->GetElements("CellData")[0];
        auto p_cell_data = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("PCellData");
        VtuOutputHelperUtilities::CreatePDataArrays(*p_cell_data, *cell_data);
        unstructured_grid_element->AddElement(p_cell_data);

        // now write all the pieces
        const auto& r_file_names = StringUtilities::SplitStringByDelimiter(list_of_file_names.str(), '\n');
        for (const auto& r_file_name : r_file_names) {
            auto piece = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("Piece");
            piece->AddAttribute("Source", r_file_name);
            unstructured_grid_element->AddElement(piece);
        }

        // writing to file
        std::stringstream output_pvtu_file_name;
        output_pvtu_file_name << rOutputFileNamePrefix << ".pvtu";
        std::ofstream output_file;
        output_file.open(output_pvtu_file_name.str(), std::ios::out | std::ios::trunc);
        VtuOutputHelperUtilities::XmlOStreamWriter writer(output_file, VtuOutput::WriterFormat::ASCII, 1);
        vtk_file_element->Write(writer);
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
        << " which is not supported. Giving priority to elements.";

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

void VtuOutput::WriteModelPart(
    const std::string& rOutputFileNamePrefix,
    ModelPart& rModelPart) const
{
    // create the vtk file
    auto vtk_file_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("VTKFile");
    vtk_file_element->AddAttribute("type", "UnstructuredGrid");
    vtk_file_element->AddAttribute("version", "0.1");
    vtk_file_element->AddAttribute("byte_order", VtuOutputHelperUtilities::GetEndianess());

    // create the unstructured grid
    auto unstructured_grid_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("UnstructuredGrid");
    vtk_file_element->AddElement(unstructured_grid_element);

    // create the piece element
    auto piece_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("Piece");
    piece_element->AddAttribute("NumberOfPoints",
                                std::to_string(rModelPart.NumberOfNodes()));
    piece_element->AddAttribute(
        "NumberOfCells",
        std::to_string(mIsElementsConsidered ? rModelPart.NumberOfElements()
                                             : rModelPart.NumberOfConditions()));
    unstructured_grid_element->AddElement(piece_element);

    // create the points element
    auto points_element = VtuOutputHelperUtilities::CreatePointsXmlElement(
        rModelPart, mIsConditionsConsidered);
    piece_element->AddElement(points_element);

    // create the cells element
    VtuOutputHelperUtilities::XmlElement::Pointer cells_element;
    if (mIsElementsConsidered) {
        cells_element = VtuOutputHelperUtilities::CreateCellsXmlElement(
            rModelPart.GetCommunicator().LocalMesh().Elements(), mKratosVtuIndicesMap);
    } else {
        cells_element = VtuOutputHelperUtilities::CreateCellsXmlElement(
            rModelPart.GetCommunicator().LocalMesh().Conditions(), mKratosVtuIndicesMap);
    }
    piece_element->AddElement(cells_element);

    // create the point data
    auto point_data_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("PointData");
    piece_element->AddElement(point_data_element);
    VtuOutputHelperUtilities::AddListOfVariables<ModelPart::NodesContainerType, ContainerDataIOTags::Historical>(
        point_data_element, rModelPart, mHistoricalVariablesMap);
    VtuOutputHelperUtilities::AddListOfVariables<ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical>(
        point_data_element, rModelPart, mNonHistoricalNodalVariablesMap);
    VtuOutputHelperUtilities::AddListOfFlags<ModelPart::NodesContainerType>(
        point_data_element, rModelPart, mNodalFlagsMap);
    VtuOutputHelperUtilities::AddListOfContainerExpressions(
        point_data_element, rModelPart, mKratosVtuIndicesMap, mPointContainerExpressionsMap);

    // create cell data
    auto cell_data_element = Kratos::make_shared<VtuOutputHelperUtilities::XmlElement>("CellData");
    piece_element->AddElement(cell_data_element);
    if (mIsElementsConsidered) {
        VtuOutputHelperUtilities::AddListOfVariables<ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical>(
            cell_data_element, rModelPart, mNonHistoricalCellVariablesMap);
        VtuOutputHelperUtilities::AddListOfFlags<ModelPart::ElementsContainerType>(
            cell_data_element, rModelPart, mCellFlagsMap);
    } else if (mIsConditionsConsidered) {
        VtuOutputHelperUtilities::AddListOfVariables<ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical>(
            cell_data_element, rModelPart, mNonHistoricalCellVariablesMap);
        VtuOutputHelperUtilities::AddListOfFlags<ModelPart::ConditionsContainerType>(
            cell_data_element, rModelPart, mCellFlagsMap);
    }

    VtuOutputHelperUtilities::AddListOfContainerExpressions(
        cell_data_element, rModelPart, mKratosVtuIndicesMap, mCellContainerExpressionsMap);

    const auto& r_communiator = rModelPart.GetCommunicator();

    std::stringstream output_vtu_file_name;
    output_vtu_file_name << rOutputFileNamePrefix;
    if (r_communiator.IsDistributed()) {
        output_vtu_file_name << "_" << r_communiator.MyPID();
    }
    output_vtu_file_name << ".vtu";

    std::ofstream output_file;
    output_file.open(output_vtu_file_name.str(), std::ios::out | std::ios::trunc);
    VtuOutputHelperUtilities::XmlOStreamWriter writer(output_file, mOutputFormat, mPrecision);
    vtk_file_element->Write(writer);
    output_file.close();

    if (r_communiator.IsDistributed()) {
        VtuOutputHelperUtilities::WritePvtuFile(vtk_file_element, rModelPart, rOutputFileNamePrefix,
                                                output_vtu_file_name.str());
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
    mPointContainerExpressionsMap.clear();
}

void VtuOutput::ClearCellContainerExpressions()
{
    mCellContainerExpressionsMap.clear();
}

void VtuOutput::PrintOutput(const std::string& rOutputFilenamePrefix)
{
    WriteModelPart(rOutputFilenamePrefix, mrModelPart);
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