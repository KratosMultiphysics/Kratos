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

KRATOS_TEST_CASE_IN_SUITE(XmlElementsArrayGetTagName, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    KRATOS_EXPECT_EQ(element.GetTagName(), "TestElement");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementsArrayAddAndGetAttributes, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    element.AddAttribute("Attribute1", "Value1");
    element.AddAttribute("Attribute2", "Value2");

    auto attributes = element.GetAttributes();

    KRATOS_EXPECT_EQ(attributes["Attribute1"], "Value1");
    KRATOS_EXPECT_EQ(attributes["Attribute2"], "Value2");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementsArrayClearAttributes, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    element.AddAttribute("Attribute1", "Value1");
    element.ClearAttributes();

    auto attributes = element.GetAttributes();
    KRATOS_EXPECT_EQ(attributes.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementsArrayAddElement, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    auto childElement = Kratos::make_shared<XmlElementsArray>("ChildElement");
    element.AddElement(childElement);

    auto children = element.GetElements();
    KRATOS_EXPECT_EQ(children[0]->GetTagName(), "ChildElement");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementsArrayGetElements, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    auto childElement1 = Kratos::make_shared<XmlElementsArray>("ChildElement");
    auto childElement2 = Kratos::make_shared<XmlElementsArray>("ChildElement2");
    element.AddElement(childElement1);
    element.AddElement(childElement2);

    auto children = element.GetElements();
    KRATOS_EXPECT_EQ(children.size(), 2);
    KRATOS_EXPECT_EQ(children[0]->GetTagName(), "ChildElement");
    KRATOS_EXPECT_EQ(children[1]->GetTagName(), "ChildElement2");
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementsArrayClearElements, KratosCoreFastSuite)
{
    XmlElementsArray element("TestElement");
    auto childElement = Kratos::make_shared<XmlElementsArray>("ChildElement");
    element.AddElement(childElement);
    element.ClearElements();

    auto children = element.GetElements();
    KRATOS_EXPECT_EQ(children.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(XmlElementsArrayWrite, KratosCoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(XmlAsciiNDDataElementChar, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2, 0);
    shape[0] = 5; shape[1] = 3;
    auto char_nd_data_1 = Kratos::make_shared<NDData<unsigned char>>(shape);
    auto char_nd_data_1_view = char_nd_data_1->ViewData();
    std::iota(char_nd_data_1_view.begin(), char_nd_data_1_view.begin() + 6, 0U);
    std::iota(char_nd_data_1_view.begin() + 6, char_nd_data_1_view.end(), 0U);

    XmlAsciiNDDataElement<unsigned char> element("data_1", char_nd_data_1, 9);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\" type=\"UInt8\">\n"
            "     0  1  2  3  4  5  0  1  2  3  4  5  6  7  8\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlAsciiNDDataElementInt, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2, 0);
    shape[0] = 5; shape[1] = 3;
    auto char_nd_data_1 = Kratos::make_shared<NDData<int>>(shape);
    auto char_nd_data_1_view = char_nd_data_1->ViewData();
    std::iota(char_nd_data_1_view.begin(), char_nd_data_1_view.begin() + 6, 0U);
    std::iota(char_nd_data_1_view.begin() + 6, char_nd_data_1_view.end(), 0U);

    XmlAsciiNDDataElement<int> element("data_1", char_nd_data_1, 9);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\" type=\"Int32\">\n"
            "     0  1  2  3  4  5  0  1  2  3  4  5  6  7  8\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlAsciiNDDataElementDouble, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2, 0);
    shape[0] = 5; shape[1] = 3;
    auto char_nd_data_1 = Kratos::make_shared<NDData<double>>(shape);
    auto char_nd_data_1_view = char_nd_data_1->ViewData();
    std::iota(char_nd_data_1_view.begin(), char_nd_data_1_view.begin() + 6, 0.0);
    std::iota(char_nd_data_1_view.begin() + 6, char_nd_data_1_view.end(), 0.0);

    XmlAsciiNDDataElement<double> element("data_1", char_nd_data_1, 1);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"ascii\" type=\"Float64\">\n"
            "     0.0e+00  1.0e+00  2.0e+00  3.0e+00  4.0e+00  5.0e+00  0.0e+00  1.0e+00  2.0e+00  3.0e+00  4.0e+00  5.0e+00  6.0e+00  7.0e+00  8.0e+00\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlBase64BinaryNDDataElementChar, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2, 0);
    shape[0] = 5; shape[1] = 3;
    auto char_nd_data_1 = Kratos::make_shared<NDData<unsigned char>>(shape);
    auto char_nd_data_1_view = char_nd_data_1->ViewData();
    std::iota(char_nd_data_1_view.begin(), char_nd_data_1_view.begin() + 6, 0U);
    std::iota(char_nd_data_1_view.begin() + 6, char_nd_data_1_view.end(), 0U);

    XmlBase64BinaryNDDataElement<unsigned char> element("data_1", char_nd_data_1);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

        KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\" type=\"UInt8\">\n"
            "     DwAAAAAAAAAAAQIDBAUAAQIDBAUGBwg=\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlBase64BinaryNDDataElementInt, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2, 0);
    shape[0] = 5; shape[1] = 3;
    auto char_nd_data_1 = Kratos::make_shared<NDData<int>>(shape);
    auto char_nd_data_1_view = char_nd_data_1->ViewData();
    std::iota(char_nd_data_1_view.begin(), char_nd_data_1_view.begin() + 6, 0U);
    std::iota(char_nd_data_1_view.begin() + 6, char_nd_data_1_view.end(), 0U);

    XmlBase64BinaryNDDataElement<int> element("data_1", char_nd_data_1);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\" type=\"Int32\">\n"
            "     PAAAAAAAAAAAAAAAAQAAAAIAAAADAAAABAAAAAUAAAAAAAAAAQAAAAIAAAADAAAABAAAAAUAAAAGAAAABwAAAAgAAAA=\n"
            "   </DataArray>\n");
}

KRATOS_TEST_CASE_IN_SUITE(XmlBase64BinaryNDDataElementDouble, KratosCoreFastSuite)
{
    DenseVector<unsigned int> shape(2, 0);
    shape[0] = 5; shape[1] = 3;
    auto char_nd_data_1 = Kratos::make_shared<NDData<double>>(shape);
    auto char_nd_data_1_view = char_nd_data_1->ViewData();
    std::iota(char_nd_data_1_view.begin(), char_nd_data_1_view.begin() + 6, 0U);
    std::iota(char_nd_data_1_view.begin() + 6, char_nd_data_1_view.end(), 0U);

    XmlBase64BinaryNDDataElement<double> element("data_1", char_nd_data_1);
    element.AddAttribute("attribute1", "value1");

    std::stringstream ss;
    element.Write(ss, 1);

    KRATOS_EXPECT_EQ(ss.str(),
            "   <DataArray Name=\"data_1\" NumberOfComponents=\"3\" attribute1=\"value1\" format=\"binary\" type=\"Float64\">\n"
            "     eAAAAAAAAAAAAAAAAAAAAAAAAAAAAPA/AAAAAAAAAEAAAAAAAAAIQAAAAAAAABBAAAAAAAAAFEAAAAAAAAAAAAAAAAAAAPA/AAAAAAAAAEAAAAAAAAAIQAAAAAAAABBAAAAAAAAAFEAAAAAAAAAYQAAAAAAAABxAAAAAAAAAIEA=\n"
            "   </DataArray>\n");
}

} // namespace Kratos::Testing