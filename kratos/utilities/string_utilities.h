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
//                   Máté Kelemen
//

#pragma once

// Core includes
#include "includes/kratos_export_api.h"
#include "includes/model_part.h"

// System includes
#include <string>
#include <vector>
#include <regex>
#include <map>
#include <filesystem>
#include <optional>

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
    std::string KRATOS_API(KRATOS_CORE) ConvertCamelCaseToSnakeCase(const std::string& rString);

    /**
     *  @brief Convert snake_case to CamelCase.
     *  @param rString String to convert.
     *  @throws If the input string
     *          - contains capital letters                              [A-Z]
     *          - contains special characters other than underscores    (?![a-z0-9_])
     *          - contains repeated underscores                         __+
     */
    std::string KRATOS_API(KRATOS_CORE) ConvertSnakeCaseToCamelCase(const std::string& rString);

    /**
     * @brief Erase first occurrence of given  substring from main string.
     * @param rMainString The string to be transformed
     * @param rToErase The string to remove
     * @return The string without the part to remove
     */
    std::string KRATOS_API(KRATOS_CORE) ErasePartialString(
        const std::string& rMainString,
        const std::string& rToErase
        );

    /**
     * @brief Checks the existence of a substring from main string.
     * @param rMainString The string to be transformed
     * @param rToCheck The string to search
     * @return True if the substring is found and false otherwise
     */
    bool KRATOS_API(KRATOS_CORE) ContainsPartialString(
        const std::string& rMainString,
        const std::string& rToCheck
        );

    /**
     * @brief This method removes whitespaces
     * @param rString The string to be transformed
     * @return The string without white spaces
     */
    std::string KRATOS_API(KRATOS_CORE) RemoveWhiteSpaces(const std::string& rString);

    /**
     * @brief This method splits a string by a delimiter
     * @param rString The string to be splitted
     * @param Delimiter The delimiter by which the string is to be splitted
     * @return a vector containing the splitted string
     */
    std::vector<std::string> KRATOS_API(KRATOS_CORE) SplitStringByDelimiter(
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
    std::string KRATOS_API(KRATOS_CORE) ReplaceAllSubstrings(
        const std::string& rInputString,
        const std::string& rStringToBeReplaced,
        const std::string& rStringToReplace
        );

}; // namespace StringUtilities


/// @addtogroup KratosCore
/// @{


/** @brief Static-only utility class for common regexes.
 *
 *  @details Static access to pairs of regex patterns (as @c std::string)
 *           and their associated regexes (as @c std::regex).
 *  @note No instances are allowed to be constructed from this class.
 */
class RegexUtility
{
public:
    ///@name Inquiry
    ///@{

    /** @brief Get a string-regex pair representing an integer.
     *
     *  @details Matches the following pattern:
     *           - optional leading '-'
     *           - 0 or [1-9]+[0-9]*
     *
     *  @return A pair of
     *          - @c std::string of a regex pattern representing an @a integer.
     *          - @c std::regex of that string.
     *  @note This function can be called safely from other static functions.
     */
    static std::pair<std::string,std::regex> Integer();

    /** @brief Get a string-regex pair representing an unsigned integer.
     *
     *  @details Matches the following pattern:
     *           - 0 or [1-9]+[0-9]*
     *
     *  @return A pair of
     *          - @c std::string of a regex pattern representing an <em>unsigned integer</em>.
     *          - @c std::regex of that string.
     *  @note This function can be called safely from other static functions.
     */
    static std::pair<std::string,std::regex> UnsignedInteger();

    /** @brief Get a string-regex pair representing an integer.
     *
     *  @details Matches the following pattern:
     *           - optional leading '-'
     *           - whole part identical to @ref RegexUtility::Integer
     *           - optional decimal point, optionally followed by the fractional part [0-9]*
     *           - optional scientific notation suffix [eE][\+-][0-9]+
     *
     *  @return A pair of
     *          - @c std::string of a regex pattern representing a <em>floating point</em>
     *            number in decimal or scientific notation.
     *          - @c std::regex of that string.
     *
     *  @note This function can be called safely from other static functions.
     */
    static std::pair<std::string,std::regex> FloatingPoint();

    ///@}

private:
    ///@name Life Cycle
    ///@{

    /// @brief Deleted default constructor (no instances of this class are allowed).
    RegexUtility() = delete;

    /// @brief Deleted move constructor (no instances of this class are allowed).
    RegexUtility(RegexUtility&& rOther) = delete;

    /// @brief Deleted copy constructor (no instances of this class are allowed).
    RegexUtility(const RegexUtility& rOther) = delete;

    ///@}
}; // class RegexUtility


/** @brief A class for interfacing placeholders and regular expressions.
 *
 *  @note placeholders should be separated by literals, otherwise
 *        regex will probably not capture them as you'd expect.
 *        <b>BAD</b> example:
 *            @code
 *            pattern:    "<placeholder_1><placeholder_2>"
 *            string:     "abcdefg"
 *            result:
 *                "<placeholder_1>" : "abcdef"
 *                "<placeholder_2>" : "g"
 *            @endcode
 *        <b>CORRECT</b> example:
 *            @code
 *            pattern:    "<placeholder_1>d<placeholder_2>"
 *            string:     "abcdefg"
 *            result:
 *                "<placeholder_1>" : "abc"
 *                "<placeholder_2>" : "efg"
 *            @endcode
 */
class KRATOS_API(KRATOS_CORE) PlaceholderPattern
{
public:
    ///@name Type Definitions
    ///@{

    using PlaceholderMap = std::map<std::string,std::string>;

    using PlaceholderGroupMap = std::map<std::string,std::optional<std::vector<std::size_t>>>;

    using MatchType = std::map<std::string,std::vector<std::string>>;

    using PathType = std::filesystem::path;

    KRATOS_CLASS_POINTER_DEFINITION(PlaceholderPattern);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    PlaceholderPattern() = default;

    /** @brief Construct from a placeholder pattern and its associated map.
     *  @param rPattern Pattern string with placeholders.
     *  @param rPlaceholderMap Pairs of placeholders and their corresponding regex strings.
     *                         Example: @code {{"<name>", ".+"}, {"<identifier>", "[0-9]+"}} @endcode
     *
     *  @warning The corresponding regexes must be bare, not containing groups (checked)
     *           or position constraints such as line begin or end modifiers (not checked).
     */
    PlaceholderPattern(const std::string& rPattern,
                       const PlaceholderMap& rPlaceholderMap);

    /// Move constructor
    PlaceholderPattern(PlaceholderPattern&& rOther) = default;

    /// Copy constructor
    PlaceholderPattern(const PlaceholderPattern& rOther) = default;

    virtual ~PlaceholderPattern() = default;

    ///@}
    ///@name Operators
    ///@{

    /// @brief Move assignment operator
    PlaceholderPattern& operator=(PlaceholderPattern&& rOther) = default;

    /// @brief Copy assignment operator
    PlaceholderPattern& operator=(const PlaceholderPattern& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Check whether a string satisfies the pattern
    bool IsAMatch(const std::string& rString) const;

    /** @brief Find all placeholders' values in the input string.
     *
     *  @param rString String matching the input pattern.
     *  @return Map associating a vector of strings, i.e. the values
     *          of placeholders in the input string, to the placeholders.
     *
     *  @note The returned placeholder values appear in the same order
     *        they do in the input pattern.
     */
    MatchType Match(const std::string& rString) const;

    /** @brief Substitute values in the stored pattern.
     *
     *  @details Return a copy of the pattern that has its placeholders replaced
     *           with the corresponding values specified in the input map.
     *  @param rPlaceholderValueMap string - string map associating values to placeholders
     *                              @code {"palceholder" : "placeholder_value"} @endcode
     */
    std::string Apply(const PlaceholderMap& rPlaceholderValueMap) const;

    /** @brief Collect all file/directory paths that match the pattern.
     *  @tparam TOutputIterator: output iterator with value type constructible from @ref PathType.
     *  @note the search begins from the filesystem root if the pattern is an absolute path,
     *        otherwise it begins from @c cwd.
     */
    template <class TOutputIterator>
    void Glob(TOutputIterator it) const;

    ///@}
    ///@name Inquiry
    ///@{

    /// @brief Return true if the input pattern contains no placeholders.
    bool IsConst() const;

    /// @brief Get the regex for the input pattern.
    const std::regex& GetRegex() const;

    /// @brief Get the string representation of the regex.
    const std::string& GetRegexString() const;

    /// @brief Get the pattern with placeholders.
    const std::string& GetPatternString() const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mPattern;

    // Placeholders and their assigned group indices in the pattern
    PlaceholderGroupMap mPlaceholderGroupMap;

    std::string mRegexString;

    std::regex mRegex;

    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Escape characters in the input that interfere with regex syntax.
    static std::string FormatRegexLiteral(const std::string& rLiteral);

    ///@}
}; // class PlaceholderPattern


/** @brief A class for working with formatted strings related to @ref ModelParts.
 *
 *  @details Operations on strings with the following placeholders are supported:
 *           - <model_part_name>
 *           - <step>
 *           - <time>
 *           - <rank>
 *           See @ref PlaceholderPattern for supported functionalities. Other
 *           placeholders can be added at compile time by tweaking the construction
 *           of the static member @ref ModelPartPattern::mModelpartPlaceholderMap.
 */
class KRATOS_API(KRATOS_CORE) ModelPartPattern : public PlaceholderPattern
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ModelPartPattern);

    ///@}
    ///@name Life Cycle
    ///@{

    ModelPartPattern() = default;

    ModelPartPattern(const std::string& rPattern);

    ModelPartPattern(ModelPartPattern&& rOther) = default;

    ModelPartPattern(const ModelPartPattern& rOther) = default;

    ///@}
    ///@name Operators
    ///@{

    using PlaceholderPattern::operator=;

    ///@}
    ///@name Operations
    ///@{

    using PlaceholderPattern::Apply;

    /** @brief Substitute values from the specified @ref ModelPart in the stored pattern.
     *
     *  @param rModelPart: Model part to extract the values of placeholders from.
     *  @note @p rModelPart must store @ref STEP and @ref TIME.
     *  @todo Add support for string formatting. Options are:
     *        - @c snprintf from the C STL (security issues due to allowing the user to freely specify the format)
     *        - @c fmt::format (requires fmtlib)
     *        - @c std::format (==fmtlib adopted in the standard, requires C++20)
     *        - @c boost::format (requires boost)
     */
    std::string Apply(const ModelPart& rModelPart) const;

    ///@}

protected:
    ///@name Protected Methods
    ///@{

    /// @brief Forwarding constructor for derived classes.
    ModelPartPattern(const std::string& rPattern, const PlaceholderMap& rPlaceholderMap);

    /// @brief Populate a key-value map of registered placeholders from a @ref ModelPart.
    virtual void PopulatePlaceholderMap(PlaceholderMap& rMap, const ModelPart& rModelPart) const;

    static PlaceholderMap GetPlaceholderMap();

    ///@}
}; // class ModelPartPattern


///@}


}  // namespace Kratos

#include "utilities/string_utilities_impl.h"
