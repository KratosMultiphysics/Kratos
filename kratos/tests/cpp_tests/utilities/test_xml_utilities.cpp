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
#include <sstream>

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/xml_utilities/xml_elements_array.h"
#include "utilities/xml_utilities/xml_ascii_nd_data_element.h"
#include "utilities/xml_utilities/xml_base64_binary_nd_data_element.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(XmlElementGetTagName, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    KRATOS_EXPECT_EQ(element.GetTagName(), "TestElement");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementAddAndGetAttributes, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    element.AddAttribute("Attribute1", "Value1");
    element.AddAttribute("Attribute2", "Value2");

    auto attributes = element.GetAttributes();

    auto itr_1 = attributes.find("Attribute1");
    KRATOS_EXPECT_FALSE(itr_1 == attributes.end());
    KRATOS_EXPECT_EQ(itr_1->second, "Value1");

    auto itr_2 = attributes.find("Attribute2");
    KRATOS_EXPECT_FALSE(itr_2 == attributes.end());
    KRATOS_EXPECT_EQ(itr_2->second, "Value2");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementClearAttributes, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    element.AddAttribute("Attribute1", "Value1");
    element.ClearAttributes();

    auto attributes = element.GetAttributes();
    KRATOS_EXPECT_EQ(attributes.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementAddElement, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    auto childElement = Kratos::make_shared<XmlElementsArray>("ChildElement");
    element.AddElement(childElement);

    auto children = element.GetElements();
    KRATOS_EXPECT_EQ(children[0]->GetTagName(), "ChildElement");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementClearElements, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    auto childElement = Kratos::make_shared<XmlElementsArray>("ChildElement");
    element.AddElement(childElement);
    element.ClearElements();

    auto children = element.GetElements();
    KRATOS_EXPECT_EQ(children.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWrite, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    element.AddAttribute("Attribute1", "Value1");
    auto childElement = Kratos::make_shared<XmlElementsArray>("ChildElement");
    element.AddElement(childElement);

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
                       "   <TestElement Attribute1=\"Value1\">\n"
                       "      <ChildElement/>\n"
                       "   </TestElement>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementAsciiChar, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2);
    shape[0] = 5;
    shape[1] = 3;
    auto p_data = Kratos::make_shared<NDData<unsigned char>>(shape);
    auto span = p_data->ViewData();
    std::size_t local_index = 0;
    for (IndexType i = 0; i < 6; ++i) {
        span[i]= local_index++;
    }
    local_index = 0;
    for (IndexType i = 6; i < 15; ++i) {
        span[i]= local_index++;
    }

    XmlAsciiNDDataElement<unsigned char> element("data_1", p_data, 4);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\" type=\"UInt8\">\n"
            "     0  1  2  3  4  5  0  1  2  3  4  5  6  7  8\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementAsciiInt, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2);
    shape[0] = 5;
    shape[1] = 3;
    auto p_data = Kratos::make_shared<NDData<int>>(shape);
    auto span = p_data->ViewData();
    std::size_t local_index = 0;
    for (IndexType i = 0; i < 6; ++i) {
        span[i]= local_index++;
    }
    local_index = 0;
    for (IndexType i = 6; i < 15; ++i) {
        span[i]= local_index++;
    }

    XmlAsciiNDDataElement<int> element("data_1", p_data, 4);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\" type=\"Int32\">\n"
            "     0  1  2  3  4  5  0  1  2  3  4  5  6  7  8\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementAsciiDouble, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2);
    shape[0] = 5;
    shape[1] = 3;
    auto p_data = Kratos::make_shared<NDData<double>>(shape);
    auto span = p_data->ViewData();
    std::size_t local_index = 0;
    for (IndexType i = 0; i < 6; ++i) {
        span[i]= local_index++;
    }
    local_index = 0;
    for (IndexType i = 6; i < 15; ++i) {
        span[i]= local_index++;
    }

    XmlAsciiNDDataElement<double> element("data_1", p_data, 1);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\" type=\"Float64\">\n"
            "     0.0e+00  1.0e+00  2.0e+00  3.0e+00  4.0e+00  5.0e+00  0.0e+00  1.0e+00  2.0e+00  3.0e+00  4.0e+00  5.0e+00  6.0e+00  7.0e+00  8.0e+00\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementBinaryChar, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2);
    shape[0] = 5;
    shape[1] = 3;
    auto p_data = Kratos::make_shared<NDData<unsigned char>>(shape);
    auto span = p_data->ViewData();
    std::size_t local_index = 0;
    for (IndexType i = 0; i < 6; ++i) {
        span[i]= local_index++;
    }
    local_index = 0;
    for (IndexType i = 6; i < 15; ++i) {
        span[i]= local_index++;
    }

    XmlBase64BinaryNDDataElement<unsigned char> element("data_1", p_data);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

        KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\" type=\"UInt8\">\n"
            "     DwAAAAABAgMEBQABAgMEBQYHCA==\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementBinaryInt, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2);
    shape[0] = 5;
    shape[1] = 3;
    auto p_data = Kratos::make_shared<NDData<int>>(shape);
    auto span = p_data->ViewData();
    std::size_t local_index = 0;
    for (IndexType i = 0; i < 6; ++i) {
        span[i]= local_index++;
    }
    local_index = 0;
    for (IndexType i = 6; i < 15; ++i) {
        span[i]= local_index++;
    }

    XmlBase64BinaryNDDataElement<int> element("data_1", p_data);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\" type=\"Int32\">\n"
            "     PAAAAAAAAAABAAAAAgAAAAMAAAAEAAAABQAAAAAAAAABAAAAAgAAAAMAAAAEAAAABQAAAAYAAAAHAAAACAAAAA==\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlOStreamWriterWriteDataElementBinaryDouble, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2);
    shape[0] = 5;
    shape[1] = 3;
    auto p_data = Kratos::make_shared<NDData<double>>(shape);
    auto span = p_data->ViewData();
    std::size_t local_index = 0;
    for (IndexType i = 0; i < 6; ++i) {
        span[i]= local_index++;
    }
    local_index = 0;
    for (IndexType i = 6; i < 15; ++i) {
        span[i]= local_index++;
    }
    XmlBase64BinaryNDDataElement<double> element("data_1", p_data);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\" type=\"Float64\">\n"
            "     eAAAAAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAQAAAAAAAAAhAAAAAAAAAEEAAAAAAAAAUQAAAAAAAAAAAAAAAAAAA8D8AAAAAAAAAQAAAAAAAAAhAAAAAAAAAEEAAAAAAAAAUQAAAAAAAABhAAAAAAAAAHEAAAAAAAAAgQA==\n"
            "   </DataArray>\n");
}

} // namespace Kratos::Testing