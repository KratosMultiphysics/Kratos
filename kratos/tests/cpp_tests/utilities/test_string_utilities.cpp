//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/string_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ConvertCammelCaseToSnakeCase, KratosCoreFastSuite)
{
    const std::string CammelCase = "TestInCammelCase";
    const std::string snake_case = StringUtilities::ConvertCammelCaseToSnakeCase(CammelCase);
    KRATOS_CHECK_STRING_EQUAL(snake_case, "test_in_cammel_case");
}

KRATOS_TEST_CASE_IN_SUITE(ErasePartialString, KratosCoreFastSuite)
{
    const std::string text_with_stuff = "textwithstuff";
    std::string text = StringUtilities::ErasePartialString(text_with_stuff, "withstuff");
    KRATOS_CHECK_STRING_EQUAL(text, "text");
    text = StringUtilities::ErasePartialString(text_with_stuff, "with");
    KRATOS_CHECK_STRING_EQUAL(text, "textstuff");
}

KRATOS_TEST_CASE_IN_SUITE(ContainsPartialString, KratosCoreFastSuite)
{
    const std::string text_with_stuff = "textwithstuff";
    bool found = StringUtilities::ContainsPartialString(text_with_stuff, "withstuff");
    KRATOS_CHECK_EQUAL(found, true);
    found = StringUtilities::ContainsPartialString(text_with_stuff, "daleatucuerpoalegriamacarena");
    KRATOS_CHECK_EQUAL(found, false);
}

KRATOS_TEST_CASE_IN_SUITE(RemoveWhiteSpaces, KratosCoreFastSuite)
{
    std::string text_with_spaces = "textwithspaces                ";
    std::string text = StringUtilities::RemoveWhiteSpaces(text_with_spaces);
    KRATOS_CHECK_STRING_EQUAL(text, "textwithspaces");
    text_with_spaces = "textwithspaces                jojobizarreadventure";
    text = StringUtilities::RemoveWhiteSpaces(text_with_spaces);
    KRATOS_CHECK_STRING_EQUAL(text, "textwithspacesjojobizarreadventure");
}

KRATOS_TEST_CASE_IN_SUITE(SplitStringByDelimiter, KratosCoreFastSuite)
{
    const std::string string_to_split = "abc, eeee,;esdsdsd,   ";
    const std::vector<std::string> splitted_string = StringUtilities::SplitStringByDelimiter(string_to_split, ',');
    KRATOS_CHECK_EQUAL(splitted_string.size(), 4);

    KRATOS_CHECK_STRING_EQUAL(splitted_string[0], "abc");
    KRATOS_CHECK_STRING_EQUAL(splitted_string[1], " eeee");
    KRATOS_CHECK_STRING_EQUAL(splitted_string[2], ";esdsdsd");
    KRATOS_CHECK_STRING_EQUAL(splitted_string[3], "   ");
}

KRATOS_TEST_CASE_IN_SUITE(ReplaceAllSubstrings, KratosCoreFastSuite)
{
    const std::string string_to_replace = "Pokemon is an awesome show. Pokemon is the best!";
    const std::string correct_string = StringUtilities::ReplaceAllSubstrings(string_to_replace, "Pokemon", "JoJo Bizarre Adventure");

    KRATOS_CHECK_STRING_EQUAL(correct_string, "JoJo Bizarre Adventure is an awesome show. JoJo Bizarre Adventure is the best!");
}

}   // namespace Testing
}  // namespace Kratos.
