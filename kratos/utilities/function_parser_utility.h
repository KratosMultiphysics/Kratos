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

#if !defined(KRATOS_FUNCTION_PARSER_UTILITY_H_INCLUDED)
#define  KRATOS_FUNCTION_PARSER_UTILITY_H_INCLUDED

// System includes
#include <cmath>
#include <unordered_map>

// External includes
#include "exprtk/exprtk.hpp"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

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
 * @class GenericFunctionUtility
 * @ingroup KratosCore
 * @brief This function allows to call a function method of type f(x, y, z, t) implemented in python style.
 * @details The functions can be constructed by providing a python-defined method of the type
 *
 *  class aux_object_cpp_callback:
 *    def __init__(self, function_string ):
 *        #TODO: check python version
 *        self.compiled_function = compile(function_string, '', 'eval', optimize=2)
 *
 *    def f(self,x,y,z,t,X,Y,Z):
 *        return  eval(self.compiled_function)
 *
 * The object is then insantiated as
 * aux_function = GenericFunctionUtility(aux_object_cpp_callback(self.function_string))
 *
 * Optionally one can specify a rotation matrix and an origin so that the function can be defined in a rotated system of coordinates
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class GenericFunctionUtility
{
public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// exprtk definitions
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>     expression_t;
    typedef exprtk::parser<double>             parser_t;

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
        ) : mFunctionBody(ReplaceAllSubstrings(rFunctionBody, "**", "^")) // Correcting function from python style to exprtk
    {
        // Defining namespace
        mNameSpace.insert(std::pair<std::string, double>("x", 0.0));
        mNameSpace.insert(std::pair<std::string, double>("y", 0.0));
        mNameSpace.insert(std::pair<std::string, double>("z", 0.0));
        mNameSpace.insert(std::pair<std::string, double>("t", 0.0));
        mNameSpace.insert(std::pair<std::string, double>("X", 0.0));
        mNameSpace.insert(std::pair<std::string, double>("Y", 0.0));
        mNameSpace.insert(std::pair<std::string, double>("Z", 0.0));

        // Defining table
        double& x = mNameSpace["x"];
        double& y = mNameSpace["y"];
        double& z = mNameSpace["z"];
        double& t = mNameSpace["t"];
        double& X = mNameSpace["X"];
        double& Y = mNameSpace["Y"];
        double& Z = mNameSpace["Z"];

        mSymbolTable.add_variable("x",x);
        mSymbolTable.add_variable("y",y);
        mSymbolTable.add_variable("z",z);
        mSymbolTable.add_variable("t",t);
        mSymbolTable.add_variable("X",X);
        mSymbolTable.add_variable("Y",Y);
        mSymbolTable.add_variable("Z",Z);

        // Creating expression
        mExpression.register_symbol_table(mSymbolTable);

        // Compiling expression
        mParser.compile(mFunctionBody, mExpression);

        // Here get the local system if it is provided
        if(LocalSystem.Has("origin")) {
            mUseLocalSystem = true;

            for(IndexType i = 0; i<3; ++i) {
                mCenterCoordinates[i] = LocalSystem["origin"][i].GetDouble();
                for(IndexType j = 0; j<3; ++j) {
                    mRotationMatrix(i, j) = LocalSystem["axes"][i][j].GetDouble();
                }
            }
        } else {
            mUseLocalSystem = false;
        }

        // Check if it depends on space
        if (mFunctionBody.find(std::string("x")) == std::string::npos &&
            mFunctionBody.find(std::string("y")) == std::string::npos &&
            mFunctionBody.find(std::string("z")) == std::string::npos &&
            mFunctionBody.find(std::string("X")) == std::string::npos &&
            mFunctionBody.find(std::string("Y")) == std::string::npos &&
            mFunctionBody.find(std::string("Z")) == std::string::npos) {
            mDependsOnSpace = false;
        }
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns if it depends on space
     * @return True if it uses the local system, false otherwise
     */
    bool UseLocalSystem()
    {
        return mUseLocalSystem;
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
    std::string FunctionBody()
    {
        return ReplaceAllSubstrings(mFunctionBody, "^", "**");
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
        )
    {
        array_1d<double,3> xglobal;
        xglobal[0] = x;
        xglobal[1] = y;
        xglobal[2] = z;
        array_1d<double,3> xlocal = prod(mRotationMatrix, (xglobal - mCenterCoordinates) );
        array_1d<double,3> xglobal_initial;
        xglobal_initial[0] = X;
        xglobal_initial[1] = Y;
        xglobal_initial[2] = Z;
        array_1d<double,3> xlocal_initial = prod(mRotationMatrix, (xglobal_initial - mCenterCoordinates) );
        return CallFunction(xlocal[0],xlocal[1],xlocal[2], t,xlocal_initial[0],xlocal_initial[1],xlocal_initial[2]);
    }

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
        )
    {
        mNameSpace["x"] = x;
        mNameSpace["y"] = y;
        mNameSpace["z"] = z;
        mNameSpace["X"] = X;
        mNameSpace["Y"] = Y;
        mNameSpace["Z"] = Z;
        mNameSpace["t"] = t;

        return mExpression.value();
    }

    ///@}

private:

    ///@name Member Variables
    ///@{

    std::unordered_map<std::string, double>  mNameSpace; /// The variables considered on the function
    symbol_table_t mSymbolTable;                         /// The symbol table of exprtk
    expression_t   mExpression;                          /// The expression of exprtk
    parser_t       mParser;                              /// The parser of exprtk
    std::string mFunctionBody;                           /// The function body

    bool mDependsOnSpace = true;                 /// If it depends on space
    bool mUseLocalSystem = false;                /// If we use a local system
    BoundedMatrix<double, 3, 3> mRotationMatrix; /// The rotation matrix
    array_1d<double, 3> mCenterCoordinates;      /// The center of coordinates

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function replaces from a string all times a certain substring is repeated
     * @param from The original string to be replaced
     * @param to The string which replaces the substring
     */
    std::string ReplaceAllSubstrings(std::string str, const std::string& from, const std::string& to) {
        std::size_t start_pos = 0;
        while((start_pos = str.find(from, start_pos)) != std::string::npos) {
            str.replace(start_pos, from.length(), to);
            start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
        }
        return str;
    }

    ///@}
}; /// GenericFunctionUtility

/**
 * @class ApplyFunctionToNodesUtility
 * @ingroup KratosCore
 * @brief This function applies a givn function to its nodes calling GenericFunctionUtility
 * @details The functions can be constructed by providing a python-defined method of the type
 * @author Riccardo Rossi
 */
class ApplyFunctionToNodesUtility
{
public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// Counted pointer of ApplyFunctionToNodesUtility
    KRATOS_CLASS_POINTER_DEFINITION(ApplyFunctionToNodesUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rNodes The nodes where to set the values
     * @param pFunction The function to set
     */

    ApplyFunctionToNodesUtility(
        ModelPart::NodesContainerType& rNodes,
        GenericFunctionUtility::Pointer pFunction
        ): mrNodes(rNodes),
           mpFunction(pFunction)
    {

    };

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method applies the function for a given variable at a certain time
     * @param rVariable The variable type
     * @param t The time variable
     */
    void ApplyFunction(
        const Variable<double>& rVariable,
        const double t
        )
    {
        // The first node iterator
        const auto it_node_begin = mrNodes.begin();

        if(!mpFunction->UseLocalSystem()) {
            //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->CallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                it_node->FastGetSolutionStepValue(rVariable) = value;
            }
        } else {
            //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->RotateAndCallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                it_node->FastGetSolutionStepValue(rVariable) = value;
            }
        }
    }

    /**
     * @brief This method returns all the evaluated values in a given time
     * @param t The time variable
     */
    std::vector<double> ReturnFunction(const double t)
    {
        // The first node iterator
        const auto it_node_begin = mrNodes.begin();

        // The vector containing the values
        std::vector<double> values(mrNodes.size());

        //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
        if(!mpFunction->UseLocalSystem()) {
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->CallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                values[k] = value;
            }
        } else {
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->RotateAndCallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                values[k] = value;
            }
        }

        return values;
    }

private:

    ///@name Protected member Variables
    ///@{

    ModelPart::NodesContainerType& mrNodes;           /// The nodes where to set the function
    GenericFunctionUtility::Pointer mpFunction; /// The function to set

    ///@}
}; /// ApplyFunctionToNodesUtility

} /// namespace Kratos

#endif // KRATOS_FUNCTION_PARSER_UTILITY_H_INCLUDED
