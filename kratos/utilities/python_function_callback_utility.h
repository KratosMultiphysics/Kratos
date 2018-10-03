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
//

#if !defined(KRATOS_PYTHON_FUNCTION_CALLBACK_UTILITY_H_INCLUDED)
#define  KRATOS_PYTHON_FUNCTION_CALLBACK_UTILITY_H_INCLUDED

#include <pybind11/pybind11.h>
#include <pybind11/eval.h>

#include <cmath>
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos
{


/**this function allows to call a function method of type f(x, y, z, t) implemented in python.
* NOTE: this makes this file to depend on python
*
* the functions can be constructed by providing a python-defined method of the type
*
*  class aux_object_cpp_callback:
*    def __init__(self, function_string ):
*        #TODO: check python version
*        self.compiled_function = compile(function_string, '', 'eval', optimize=2)
*
*    def f(self,x,y,z,t):
*        return  eval(self.compiled_function)
*
* the object is then insantiated as
* aux_function = PythonGenericFunctionUtility(aux_object_cpp_callback(self.function_string))
*
* optionally one can specify a rotation matrix and an origin so that the function can be defined in a rotated system of coordinates
*/
class PythonGenericFunctionUtility
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PythonGenericFunctionUtility);

    PythonGenericFunctionUtility(  const std::string& function_body,  Parameters local_system = Parameters{} )
    {
        //compile the function starting from the string function body
        try
        {
            main_module = pybind11::module::import("__main__");
            main_namespace = main_module.attr("__dict__");
            pybind11::exec("from math import *", main_namespace);
            mfunction_body = function_body;
            //mbytecode = pybind11::object( pybind11::handle<>( (PyObject*) Py_CompileString(function_body.c_str(), "pyscript", Py_eval_input) ) );
        }
        catch(pybind11::error_already_set const&)
        {
            PyErr_Print();
        }


        // here get the local system if it is provided
        if(local_system.Has("origin"))
        {
            muse_local_system = true;

            for(unsigned int i = 0; i<3; ++i)
            {
                mxc[i] = local_system["origin"][i].GetDouble();
                for(unsigned int j = 0; j<3; ++j)
                {
                    mR(i, j) = local_system["axes"][i][j].GetDouble();
                }
            }
        }
        else
        {
            muse_local_system = false;
        }

        // check if it depends on space
        if (function_body.find(std::string("x"))==std::string::npos &&
                function_body.find(std::string("y")) ==std::string::npos &&
                function_body.find(std::string("z")) ==std::string::npos)
        {
            mdepends_on_space = false;
        }
    }

//         PythonGenericFunctionUtility(  pybind11::object obj, const Matrix& R, const Vector& xc): mpy_obj(obj), muse_local_system(true), mR(R), mxc(xc)
//         {
//             if(mR.size1() != 3 or mR.size2() != 3)
//                 KRATOS_ERROR << "rotation matrix is expected to have size 3. The matrix given is " << mR;
//             if(mxc.size() != 3)
//                 KRATOS_ERROR << "center position is expected to have size 3. The position given is " << mxc;
//
//             // TODO: check if the matrix R is orthonormal
//         }

    bool UseLocalSystem()
    {
        return muse_local_system;
    }
    bool DependsOnSpace()
    {
        return mdepends_on_space;
    }

    double RotateAndCallFunction(const double x, const double y, const double z, const double t)
    {
        array_1d<double,3> xglobal;
        xglobal[0] = x;
        xglobal[1] = y;
        xglobal[2] = z;
        array_1d<double,3> xlocal = prod(mR, (xglobal - mxc) );
        return CallFunction(xlocal[0],xlocal[1],xlocal[2],t);
    }

    double CallFunction(const double x, const double y, const double z, const double t)
    {
        main_namespace["x"] = x;
        main_namespace["y"] = y;
        main_namespace["z"] = z;
        main_namespace["t"] = t;


//         #if PY_MAJOR_VERSION >= 3
//         PyObject* res = PyEval_EvalCode(mbytecode.ptr(),main_namespace.ptr(),main_namespace.ptr());
//         #else
//         PyObject* res = PyEval_EvalCode((PyCodeObject*)(mbytecode.ptr()),main_namespace.ptr(),main_namespace.ptr());
//         #endif
        return pybind11::eval(mfunction_body, main_namespace).cast<double>();
    }


private:
    pybind11::object main_module;
    pybind11::object main_namespace;
    std::string mfunction_body;
    pybind11::object mbytecode;

    bool mdepends_on_space = true;
    bool muse_local_system = false;
    BoundedMatrix<double, 3, 3> mR;
    array_1d<double, 3> mxc;
};

class ApplyFunctionToNodesUtility
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyFunctionToNodesUtility);

    ApplyFunctionToNodesUtility(  ModelPart::NodesContainerType& rNodes, PythonGenericFunctionUtility::Pointer pfunction): mrNodes(rNodes), mpfunction(pfunction)
    {};

    template < class TVarType >
    void ApplyFunction(const TVarType& rVariable, const double t)
    {
        if(mpfunction->UseLocalSystem() == false)
        {
            //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
            for (int k = 0; k< static_cast<int> (mrNodes.size()); k++)
            {
                ModelPart::NodesContainerType::iterator i = mrNodes.begin() + k;
                const double value = mpfunction->CallFunction(i->X(), i->Y(), i->Z(), t);
                i->FastGetSolutionStepValue(rVariable) = value;
            }
        }
        else
        {
            //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
            for (int k = 0; k< static_cast<int> (mrNodes.size()); k++)
            {
                ModelPart::NodesContainerType::iterator i = mrNodes.begin() + k;
                const double value = mpfunction->RotateAndCallFunction(i->X(), i->Y(), i->Z(), t);
                i->FastGetSolutionStepValue(rVariable) = value;
            }
        }
    }


    std::vector <double> ReturnFunction(const double t)
    {
        std::vector<double> values;
        //WARNING: do NOT put this loop in parallel, the python GIL does not allow you to do it!!
        if(mpfunction->UseLocalSystem() == false)
        {
            for (int k = 0; k< static_cast<int> (mrNodes.size()); k++)
            {
                ModelPart::NodesContainerType::iterator i = mrNodes.begin() + k;
                const double value = mpfunction->CallFunction(i->X(), i->Y(), i->Z(), t);
                values.push_back(value);
            }
        }
        else
        {
            for (int k = 0; k< static_cast<int> (mrNodes.size()); k++)
            {
                ModelPart::NodesContainerType::iterator i = mrNodes.begin() + k;
                const double value = mpfunction->RotateAndCallFunction(i->X(), i->Y(), i->Z(), t);
                values.push_back(value);
            }
        }

        return values;
    }

private:
    ModelPart::NodesContainerType& mrNodes;
    PythonGenericFunctionUtility::Pointer mpfunction;
};

}

#endif // KRATOS_PYTHON_FUNCTION_CALLBACK_UTILITY_H_INCLUDED
