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

// External includes

// Project includes
#include "utilities/function_parser_utility.h"

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
 * @note This is a legacy file in order keep retrocompatibility
 * @author Riccardo Rossi
 */
class PythonGenericFunctionUtility
    : public GenericFunctionUtility
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
        ) : GenericFunctionUtility(rFunctionBody, LocalSystem)
    {
        KRATOS_WARNING("PythonGenericFunctionUtility") <<  "This is a legacy class. Please use GenericFunctionUtility" << std::endl;
    }

}; /// PythonGenericFunctionUtility

} /// namespace Kratos

#endif // KRATOS_PYTHON_FUNCTION_CALLBACK_UTILITY_H_INCLUDED
