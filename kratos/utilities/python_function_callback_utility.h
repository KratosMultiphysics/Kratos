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

#if !defined(KRATOS_PYTHON_FUNCTION_CALLBACK_UTILITY_H_INCLUDED)
#define  KRATOS_PYTHON_FUNCTION_CALLBACK_UTILITY_H_INCLUDED

// System includes
#include <cmath>

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/eval.h>

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
 * @class PythonGenericFunctionUtility
 * @ingroup KratosCore
 * @brief This function allows to call a function method of type f(x, y, z, t) implemented in python.
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
 * aux_function = PythonGenericFunctionUtility(aux_object_cpp_callback(self.function_string))
 *
 * Optionally one can specify a rotation matrix and an origin so that the function can be defined in a rotated system of coordinates
 * @note This makes this file to depend on python
 * @author Riccardo Rossi
 */
class PythonGenericFunctionUtility
{
public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// Counted pointer of PythonGenericFunctionUtility
    KRATOS_CLASS_POINTER_DEFINITION(PythonGenericFunctionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rFunctionBody The string defining the function
     * @param LocalSystem The parameters defining the local system
     */
    PythonGenericFunctionUtility(
        const std::string& rFunctionBody,
        Parameters LocalSystem = Parameters{}
        ) : mFunctionBody(rFunctionBody)
    {
        // Compile the function starting from the string function body
        try {
            mMainModule = pybind11::module::import("__main__");
            mMainNameSpace = mMainModule.attr("__dict__");
            pybind11::exec("from math import *", mMainNameSpace);
            // TODO: this is commented. if we are able to compile the function it will be easier anc cheaper to evaluate the function
//             mByteCode = pybind11::object( pybind11::handle<>( (PyObject*) Py_CompileString(rFunctionBody.c_str(), "pyscript", Py_eval_input) ) );
        } catch(pybind11::error_already_set const&) {
            PyErr_Print();
        }

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
        if (rFunctionBody.find(std::string("x")) == std::string::npos &&
            rFunctionBody.find(std::string("y")) == std::string::npos &&
            rFunctionBody.find(std::string("z")) == std::string::npos &&
            rFunctionBody.find(std::string("X")) == std::string::npos &&
            rFunctionBody.find(std::string("Y")) == std::string::npos &&
            rFunctionBody.find(std::string("Z")) == std::string::npos) {
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
        return mFunctionBody;
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
        mMainNameSpace["x"] = x;
        mMainNameSpace["y"] = y;
        mMainNameSpace["z"] = z;
        mMainNameSpace["X"] = X;
        mMainNameSpace["Y"] = Y;
        mMainNameSpace["Z"] = Z;
        mMainNameSpace["t"] = t;

        return pybind11::eval(mFunctionBody, mMainNameSpace).cast<double>();
    }

    ///@}

private:

    ///@name Member Variables
    ///@{

    pybind11::object mMainModule;       /// The main python module
    pybind11::object mMainNameSpace;    /// The main python namespace (the variables considered on the python function)
    const std::string mFunctionBody;    /// The function body
//     pybind11::object mByteCode;         /// Some byte code

    bool mDependsOnSpace = true;                 /// If it depends on space
    bool mUseLocalSystem = false;                /// If we use a local system
    BoundedMatrix<double, 3, 3> mRotationMatrix; /// The rotation matrix
    array_1d<double, 3> mCenterCoordinates;      /// The center of coordinates

    ///@}
}; /// PythonGenericFunctionUtility

/**
 * @class ApplyFunctionToNodesUtility
 * @ingroup KratosCore
 * @brief This function applies a givn function to its nodes calling PythonGenericFunctionUtility
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
        PythonGenericFunctionUtility::Pointer pFunction
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
     * @tparam TVarType The type of variable considered
     */
    template<class TVarType>
    void ApplyFunction(
        const TVarType& rVariable,
        const double t,
        const unsigned int Step = 0)
    {
        // The first node iterator
        const auto it_node_begin = mrNodes.begin();

        if(!mpFunction->UseLocalSystem()) {
            //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->CallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                it_node->FastGetSolutionStepValue(rVariable, Step) = value;
            }
        } else {
            //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
            for (IndexType k = 0; k< mrNodes.size(); k++) {
                auto it_node = it_node_begin + k;
                const double value = mpFunction->RotateAndCallFunction(it_node->X(), it_node->Y(), it_node->Z(), t, it_node->X0(), it_node->Y0(), it_node->Z0());
                it_node->FastGetSolutionStepValue(rVariable, Step) = value;
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
    PythonGenericFunctionUtility::Pointer mpFunction; /// The function to set

    ///@}
}; /// ApplyFunctionToNodesUtility

} /// namespace Kratos

#endif // KRATOS_PYTHON_FUNCTION_CALLBACK_UTILITY_H_INCLUDED
