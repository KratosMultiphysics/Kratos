//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

// Forward declaration to avoid the include of tinyexpr in a header file
struct te_expr;

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
 * @class BasicGenericFunctionUtility
 * @ingroup KratosCore
 * @brief This function allows to call a function method of type f(x, y, z, t) implemented in python style.
 * @details The functions can be constructed by providing a python-defined method of the type
 *
 *
 * The object is then insantiated as
 * aux_function = BasicGenericFunctionUtility(aux_object_cpp_callback(self.function_string))
 *
 * Optionally one can specify a rotation matrix and an origin so that the function can be defined in a rotated system of coordinates
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) BasicGenericFunctionUtility
{
public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// Counted pointer of BasicGenericFunctionUtility
    KRATOS_CLASS_POINTER_DEFINITION(BasicGenericFunctionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rFunctionBody The string defining the function
     */
    BasicGenericFunctionUtility(const std::string& rFunctionBody);

    ///Copy constructor
    BasicGenericFunctionUtility(BasicGenericFunctionUtility const& rOther);

    /// Destructor.
    virtual ~BasicGenericFunctionUtility();

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method clones the current function instance
     * @return A clone of the current class
     */
    Pointer Clone()
    {
        return Kratos::make_shared<BasicGenericFunctionUtility>(*this);
    }

    /**
     * @brief This method returns if it depends on space
     * @return True if it uses the local system, false otherwise
     */
    virtual bool UseLocalSystem()
    {
        return false;
    }

    /**
     * @brief This method returns if it depends on space
     * @return True if it depends on space, false otherwise
     */
    bool DependsOnSpace()
    {
        return mDependsOnSpace;
    }

    /**
     * @brief This method returns the function body
     * @return The function body
     */
    std::string FunctionBody();

    /**
     * @brief This method rotates and calls the evaluation function
     * @param x The x coordinate
     * @param y The y coordinate
     * @param z The z coordinate
     * @param t The time variable
     * @param X The initial x coordinate
     * @param Y The initial y coordinate
     * @param Z The initial z coordinate
     */
    virtual double RotateAndCallFunction(
        const double x,
        const double y,
        const double z,
        const double t,
        const double X = 0.0,
        const double Y = 0.0,
        const double Z = 0.0
        );

    /**
     * @brief This calls the evaluation function
     * @param x The x coordinate
     * @param y The y coordinate
     * @param z The z coordinate
     * @param t The time variable
     * @param X The initial x coordinate
     * @param Y The initial y coordinate
     * @param Z The initial z coordinate
     */
    double CallFunction(
        const double x,
        const double y,
        const double z,
        const double t,
        const double X = 0.0,
        const double Y = 0.0,
        const double Z = 0.0
        );

    ///@}

protected:

    ///@name Member Variables
    ///@{

    array_1d<double, 7> mValues = ZeroVector(7);                          /// The variables values considered on the function
    std::vector<te_expr*> mpTinyExpr = std::vector<te_expr*>(1, nullptr); /// The function parser
    std::string mFunctionBody;                                            /// The function body

    bool mDependsOnSpace = true;                                          /// If it depends on space

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method initializes the parser classes
     */
    void InitializeParser();

    ///@}
}; /// GenericFunctionUtility

/**
 * @class GenericFunctionUtility
 * @ingroup KratosCore
 * @brief This function allows to call a function method of type f(x, y, z, t) implemented in python style (with rotations).
 * @details The functions can be constructed by providing a python-defined method of the type
 *
 *
 * The object is then insantiated as
 * aux_function = GenericFunctionUtility(aux_object_cpp_callback(self.function_string))
 *
 * Optionally one can specify a rotation matrix and an origin so that the function can be defined in a rotated system of coordinates
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) GenericFunctionUtility
  : public BasicGenericFunctionUtility
{
public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// Counted pointer of GenericFunctionUtility
    KRATOS_CLASS_POINTER_DEFINITION(GenericFunctionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rFunctionBody The string defining the function
     * @param LocalSystem The parameters defining the local system
     */
    GenericFunctionUtility(
        const std::string& rFunctionBody,
        Parameters LocalSystem = Parameters{}
        );

    ///Copy constructor
    GenericFunctionUtility(GenericFunctionUtility const& rOther);

    /// Destructor.
    ~GenericFunctionUtility() override {};

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method clones the current function instance
     * @return A clone of the current class
     */
    Pointer Clone()
    {
        return Kratos::make_shared<GenericFunctionUtility>(*this);
    }

    /**
     * @brief This method returns if it depends on space
     * @return True if it uses the local system, false otherwise
     */
    bool UseLocalSystem() override
    {
        return mUseLocalSystem;
    }

    /**
     * @brief This method rotates and calls the evaluation function
     * @param x The x coordinate
     * @param y The y coordinate
     * @param z The z coordinate
     * @param t The time variable
     * @param X The initial x coordinate
     * @param Y The initial y coordinate
     * @param Z The initial z coordinate
     */
    double RotateAndCallFunction(
        const double x,
        const double y,
        const double z,
        const double t,
        const double X = 0.0,
        const double Y = 0.0,
        const double Z = 0.0
        ) override;

    ///@}

private:

    ///@name Member Variables
    ///@{

    bool mUseLocalSystem = false;                /// If we use a local system
    BoundedMatrix<double, 3, 3> mRotationMatrix; /// The rotation matrix
    array_1d<double, 3> mCenterCoordinates;      /// The center of coordinates

    ///@}
}; /// GenericFunctionUtility

/**
 * @class FunctionParser
 * @ingroup KratosCore
 * @brief This class parses a string as a std::function
 * @author Vicente Mataix Ferrandiz
 */
class FunctionParser
{
public:
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rFunctionBody The string defining the function
     */
    FunctionParser(const std::string& rFunctionBody)
        : mParser(rFunctionBody)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method parses a string into a std::function
     * @details This is the version with 3 arguments (x,y,z)
     * @param rFunctionBody The string defining the function
     */
    std::function<double(const double, const double, const double)> GetFunctionSpace()
    {
        return [&](const double x, const double y, const double z) -> double {return mParser.CallFunction(x,y,z,0.0);};
    }

    /**
     * @brief This method parses a string into a std::function
     * @details This is the version with 4 arguments (x,y,z,t)
     * @param rFunctionBody The string defining the function
     */
    std::function<double(const double, const double, const double, const double)> GetFunction()
    {
        return [&](const double x, const double y, const double z, const double t) -> double {return mParser.CallFunction(x,y,z,t);};
    }

    /**
     * @brief This method parses a string into a std::function
     * @details This is the version with 7 arguments (x,y,z,t)
     * @param rFunctionBody The string defining the function
     */
    std::function<double(const double, const double, const double, const double, const double, const double, const double)> GetFunctionInitialCoordinates()
    {
        return [&](const double x, const double y, const double z, const double t, const double X, const double Y, const double Z) -> double {return mParser.CallFunction(x,y,z,t,X,Y,Z);};;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    GenericFunctionUtility mParser; /// The parser function utility

    ///@}
}; /// FunctionParser

} /// namespace Kratos
