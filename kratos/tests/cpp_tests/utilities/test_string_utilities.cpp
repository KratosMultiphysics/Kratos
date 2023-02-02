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

TEST(ConvertCammelCaseToSnakeCase, KratosCoreFastSuite)
{
    const std::string CammelCase = "TestInCammelCase";
    const std::string snake_case = StringUtilities::ConvertCammelCaseToSnakeCase(CammelCase);
    KRATOS_EXPECT_EQ(snake_case, "test_in_cammel_case");
}

TEST(ErasePartialString, KratosCoreFastSuite)
{
    const std::string text_with_stuff = "textwithstuff";
    std::string text = StringUtilities::ErasePartialString(text_with_stuff, "withstuff");
    KRATOS_EXPECT_EQ(text, "text");
    text = StringUtilities::ErasePartialString(text_with_stuff, "with");
    KRATOS_EXPECT_EQ(text, "textstuff");
}

TEST(ContainsPartialString, KratosCoreFastSuite)
{
    const std::string text_with_stuff = "textwithstuff";
    bool found = StringUtilities::ContainsPartialString(text_with_stuff, "withstuff");
    KRATOS_EXPECT_EQ(found, true);
    found = StringUtilities::ContainsPartialString(text_with_stuff, "daleatucuerpoalegriamacarena");
    KRATOS_EXPECT_EQ(found, false);
}

TEST(RemoveWhiteSpaces, KratosCoreFastSuite)
{
    std::string text_with_spaces = "textwithspaces                ";
    std::string text = StringUtilities::RemoveWhiteSpaces(text_with_spaces);
    KRATOS_EXPECT_EQ(text, "textwithspaces");
    text_with_spaces = "textwithspaces                jojobizarreadventure";
    text = StringUtilities::RemoveWhiteSpaces(text_with_spaces);
    KRATOS_EXPECT_EQ(text, "textwithspacesjojobizarreadventure");
}

TEST(SplitStringByDelimiter, KratosCoreFastSuite)
{
    const std::string string_to_split = "abc, eeee,;esdsdsd,   ";
    const std::vector<std::string> splitted_string = StringUtilities::SplitStringByDelimiter(string_to_split, ',');
    KRATOS_EXPECT_EQ(splitted_string.size(), 4);

    KRATOS_EXPECT_EQ(splitted_string[0], "abc");
    KRATOS_EXPECT_EQ(splitted_string[1], " eeee");
    KRATOS_EXPECT_EQ(splitted_string[2], ";esdsdsd");
    KRATOS_EXPECT_EQ(splitted_string[3], "   ");
}

TEST(ReplaceAllSubstrings, KratosCoreFastSuite)
{
    const std::string string_to_replace = "Pokemon is an awesome show. Pokemon is the best!";
    const std::string correct_string = StringUtilities::ReplaceAllSubstrings(string_to_replace, "Pokemon", "JoJo Bizarre Adventure");

    KRATOS_EXPECT_EQ(correct_string, "JoJo Bizarre Adventure is an awesome show. JoJo Bizarre Adventure is the best!");
}

}   // namespace Testing
}  // namespace Kratos.
