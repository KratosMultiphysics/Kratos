//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @namespace StringUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for computing string operations
 * @author Vicente Mataix Ferrandiz
 */
namespace StringUtilities
{
    /**
     * @brief This method converts CamelCase to snake_case
     * @param rString The string to be transformed into snake_case
     * @return The string in snake_case
     */
    [[nodiscard]] std::string KRATOS_API(KRATOS_CORE) ConvertCamelCaseToSnakeCase(const std::string& rString);

    /**
     *  @brief Convert snake_case to CamelCase.
     *  @param rString String to convert.
     *  @throws If the input string
     *          - contains capital letters                              [A-Z]
     *          - contains special characters other than underscores    (?![a-z0-9_])
     *          - contains repeated underscores                         __+
     */
    [[nodiscard]] std::string KRATOS_API(KRATOS_CORE) ConvertSnakeCaseToCamelCase(const std::string& rString);

    /**
     * @brief Erase first occurrence of given  substring from main string.
     * @param rMainString The string to be transformed
     * @param rToErase The string to remove
     * @return The string without the part to remove
     */
    [[nodiscard]] std::string KRATOS_API(KRATOS_CORE) ErasePartialString(
        const std::string& rMainString,
        const std::string& rToErase
        );

    /**
     * @brief Checks the existence of a substring from main string.
     * @param rMainString The string to be transformed
     * @param rToCheck The string to search
     * @return True if the substring is found and false otherwise
     */
    [[nodiscard]] bool KRATOS_API(KRATOS_CORE) ContainsPartialString(
        const std::string& rMainString,
        const std::string& rToCheck
        );

    /**
     * @brief This method removes whitespaces
     * @param rString The string to be transformed
     * @return The string without white spaces
     */
    [[nodiscard]] std::string KRATOS_API(KRATOS_CORE) RemoveWhiteSpaces(const std::string& rString);

    /**
     * @brief This method splits a string by a delimiter
     * @param rString The string to be splitted
     * @param Delimiter The delimiter by which the string is to be splitted
     * @return a vector containing the splitted string
     */
    [[nodiscard]] std::vector<std::string> KRATOS_API(KRATOS_CORE) SplitStringByDelimiter(
        const std::string& rString,
        const char Delimiter
        );

    /**
     * @brief This function replaces from a string all times a certain substring is repeated
     * @param rInputString The input string to replace the substring
     * @param rStringToBeReplaced The original string to be replaced
     * @param rStringToReplace The string which replaces the substring
     * @return The string updated with the new substring
     */
    [[nodiscard]] std::string KRATOS_API(KRATOS_CORE) ReplaceAllSubstrings(
        const std::string& rInputString,
        const std::string& rStringToBeReplaced,
        const std::string& rStringToReplace
        );

    /**
     * @brief This function trims a string by removing whitespaces, tabs etc from left and right. Same as "strip" in Python
     * @param rInputString The input string to trim
     * @param RemoveNullChar Whether or not null-characters ('\0') should be removed
     * @return The trimmed string
     */
    [[nodiscard]] std::string KRATOS_API(KRATOS_CORE) Trim(
        const std::string& rInputString,
        const bool RemoveNullChar = false);

    /**
     * @brief This function trims a string by removing whitespaces, tabs etc from left. Same as "lstrip" in Python
     * @param rInputString The input string to trim
     * @param RemoveNullChar Whether or not null-characters ('\0') should be removed
     * @return The trimmed string
     */
    [[nodiscard]] std::string KRATOS_API(KRATOS_CORE) TrimLeft(
        const std::string& rInputString,
        const bool RemoveNullChar = false);

    /**
     * @brief This function trims a string by removing whitespaces, tabs etc and right. Same as "rstrip" in Python
     * @param rInputString The input string to trim
     * @param RemoveNullChar Whether or not null-characters ('\0') should be removed
     * @return The trimmed string
     */
    [[nodiscard]] std::string KRATOS_API(KRATOS_CORE) TrimRight(
        const std::string& rInputString,
        const bool RemoveNullChar = false);


    /**
     * @brief Prints the data of an object of type TClass to the given output stream with indentation.
     * @param rOStream The output stream where the data will be printed.
     * @param rThisClass The object of type TClass whose data will be printed.
     * @param Identation (optional) The string used for the indentation. Default is four spaces.
     */
    template<class TClass>
    static void PrintDataWithIdentation(
        std::ostream& rOStream,
        const TClass& rThisClass,
        const std::string Identation = "\t"
        )
    {
        // Auxiliary stream and line
        std::stringstream ss;
        std::string line;
        rThisClass.PrintData(ss);

        // Get the output string from the stringstream.
        const std::string& r_output = ss.str();

        // Split the output string into lines.
        std::istringstream iss(r_output);
        while (std::getline(iss, line)) {
            // Here, 'line' is a single line from the output.
            rOStream << Identation << line << "\n";
        }
    }

}; // namespace StringUtilities
}  // namespace Kratos