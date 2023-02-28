---
title: Documenting C++ Code
keywords: documentation, doxygen, c++
tags: [documenting_cpp.md]
sidebar: documentation_guide
summary: 
---
## Intro

This page is meant to provide a quick primer on how <a>Kratos</a> code should be documented in-source, using doxygen.

Doxygen is the most popular tool for generating documentation from annotated source files of major languages (such as C++ or Python). It performs static code analysis on the repository to deduce a lot useful information from namespace contents to class composition and inheritance hierarchies. Although Doxygen can do this automatically without any comments or annotations, the proper description of classes and functions is essential for a useful documentation, which in turn is the cornerstone of attracting new users and developers to the project (and remembering what you were thinking a month ago when you wrote that class template).

## Documenting C++ Code

Doxygen is compatible with a vast variety of languages and the basic structure of annotating the code does not change, though minor syntax variations can occur between them. This section introduces essential tags and commands used throughout <a>Kratos</a>, with examples in C++.

Doxygen annotations are distinguished from regular comments by appending an extra comment character (```/``` or ```*```) to the comment signal (```//``` or ```/*``` respectively). Consecutive lines with annotations are considered to be part of the same annotation block, but regular comments without the exta characters are not included in the generated documentation.
```cpp
// A regular comment ignored by doxygen.

/// Single-line doxygen annotation.

/** Scoped doxygen annotation */

/// Another single-line annotation.
/** Followed by a multiline
 *  scoped doxygen annotation,
 *  in the same block.
 */
```

In C++, annotations should **directly precede** the constructs they document (**variable declarations**, **function declarations** or **class definitions**), and as such, should mostly be in header files. In case of forward-declared classes, annotations should only be included once, in front of the class definition.
```cpp
// A forward declaration of Node. Don't annotate this.
struct Node;

/** @brief Print data stored in a @ref Node to an output stream.
 *
 *  @param rStream: output stream to print the data to.
 *  @param rNode: node to extract data from.
 */
void PrintNode(std::ostream& rStream, const Node& rNode);

// Definition of Node => annotate it!
/// Simple 2D node, storing only its position and ID.
struct Node
{
    ///@name Data members
    ///@{

    /// Position on the x-axis.
    double m_x;

    /// Position on the y-axis.
    double m_y;

    ///@}

    ///@name Miscellaneous members
    ///@{

    /// Unique identifier of the node.
    unsigned int m_id;

    ///@}
};
```

[Here](../../../../../external_data/documentation_guide/example_1/index.html)'s what the generated documentation looks like for the example above. Notice the tags identified with <i>@</i> characters in the annotations; these are doxygen commands that relay formatting directives. They always modify the appearance of the text that follows them, but the number of characters/words they effect varies between commands. For example, the ```@ref``` command creates a link to the documentation of the name that follows it (```@ref Node``` links to ```Node``` in the generated documentation), while ```@param``` describes an argument of the function, taking a parameter name and a possibly multi-line description. There's also ```@name``` that tags every construct in its scope enclosed between ```@{``` and ```@}```. You can browse through the [complete list of doxygen commands](https://www.doxygen.nl/manual/commands.html) to get an idea of what you can work with, but here's a short list often used in <i>Kratos</i>:
- ```@brief```: brief description of the class/function/variable/concept that always gets displayed next to its name.
- ```@details```: detailed description of the class/function/variable/concept that is only displayed on the construct's own documentation page.
- ```@addtogroup```: add all constructs defined in the scope to a doxygen *module* (each *Kratos* application has its own *module*).
- ```@name```: Tag all constructs defined in the scope with a name, creating a separate paragraph for the in the documentation page. This command is used for member constructs within class definitions.
- ```@f(``` and ```@f)```: in-line latex equation.
- ```@f[``` and ```@f]```: latex equation on a new line.
- ```@param```: function argument description.
- ```@tparam```: template argument description
- ```@return```: return value description
- ```@ref```: create a link in the documentation to the name that follows it.
- ```@p```: reference to an argument local to the current annotation block.
- ```@a```: render the following name in italic.
- ```@b```: render the following name in bold.
- ```@c```: render the following name in typewriter font.


Let's take a deep dive and look at a reasonably well-decorated piece of code with lots of doxygen commands. Be sure to compare the source and the [generated documentation](../../../../../external_data/documentation_guide/example_2/index.html)
```cpp
/** @defgroup CompileTimeApplication
 *  @{
 *      @page Kratos Compile-time Application
 *      Utilities for performing calculations and logic
 *      at compile time to help with template programming.
 *  @}
 */

namespace Kratos {

///@addtogroup CompileTimeApplication
///@{

/** @brief Raise any base to an integer power.
 *
 *  @details Compute \f( a^n \f).
 *
 *  @tparam TNumeric: base type; can be any integer or floating point type (example: @c int or @c double).
 *  @param base: base to be raised to a power.
 *  @param exponent: multiply @p base by itself this many times
 *  @note this function can be invoked at compile time,
 *        and is guaranteed no to throw an exception.
 */
template <typename TNumeric>
constexpr TNumeric IntegerPower(TNumeric base, std::size_t exponent) noexcept
{
    if (exponent == 0) {
        return 1;
    }
    else {
        auto number_of_mults = 1;
        TNumeric power = base;
        while (number_of_mults < exponent) {
            if ((exponent - number_of_mults) / number_of_mults) {
                // Enough room to square the current state.
                power *= power;
                number_of_mults *= 2;
            }
            else {
                // Cannot square the current state
                // => multiply once by base
                power *= base;
                ++number_of_mults;
            }
        } // while

        return power;
    } // exponent != 0
} // IntegerPower()


/// @brief Compute the sum of the first @a n terms of a geometric series.
///
/// @details Compute the following expression: \f[ c \cdot \sum_{k=0}^n q^k \f]
///
/// @param coefficient: (@a c) coefficient multiplying each term of the series.
/// @param base: (@a q) base of the geometric series.
/// @param max_terms: (@a n) number of terms to sum up.
/// @note the classic reduced formula is used to compute the sum, unless @c base is @c 1.
///       Power calculation is deferred to @ref IntegerPower<double> which is less efficient
///       than @c std::pow but can be invoked at compile time and does not throw exceptions.
/// @warning be wary of floating point overflows.
constexpr double ComputeGeometricSum(double coefficient, double base, std::size_t max_terms) noexcept
{
    return coefficient * (base == 1 ? base * max_terms : (1 - IntegerPower(base, max_terms)) / (1 - base));
}

///@}

} // namespace Kratos
```

### Advice:
- Don't put annotations **inside function definitions**. Document everything in one block preceding the declaration.
- Define ```@addtogroup``` **inside namespaces** not outside of them, otherwise the constructs will be added to the documentation of the namespace but not the group.
- Always provide a brief description, and tag it with ```@brief```. It's very useful in the documentation.
- When using ```@ref``` to link to constructs in another namespace, make sure to specify the namespace as well. Doxygen won't know which construct to link to if there are more of them with the same name in different namespaces (or different languages - think of C++ classes and their bindings in python).
- If you're writing functions/classes that should not be in the documentation (helpers, implementation details), enclose them in ```namespace``` ```detail``` or ```impl``` if they're in headers, or [unnamed namespaces](https://en.cppreference.com/w/cpp/language/namespace#Unnamed_namespaces) if they're in source files.


