// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#if !defined(KRATOS_SOLUTION_STEP_VARIABLE_INDEXING_PYTHON_H_INCLUDED )
#define  KRATOS_SOLUTION_STEP_VARIABLE_INDEXING_PYTHON_H_INCLUDED



// System includes
//#include <string>
//#include <iostream>
//#include <sstream>


// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;
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

/// Short class definition.
/** Detail class definition.
*/
template<class TContainerType, class TVariableType>
class SolutionStepVariableIndexingPython : public def_visitor<SolutionStepVariableIndexingPython<TContainerType, TVariableType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SolutionStepVariableIndexingPython
    KRATOS_CLASS_POINTER_DEFINITION(SolutionStepVariableIndexingPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SolutionStepVariableIndexingPython() {}

    /// Copy constructor.
    SolutionStepVariableIndexingPython(const SolutionStepVariableIndexingPython& rOther);

    /// Destructor.
    virtual ~SolutionStepVariableIndexingPython() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template <class TClassType>
    void visit(TClassType& ThisClass) const
    {
        ThisClass
        //.def("__contains__", &DataValueContainerHas)
        //.def("__setitem__", &SetBufferValue)
        //.def("__getitem__", &GetBufferValue)
        .def("Has", &DataValueContainerHas)
        .def("SetSolutionStepValue", &SetValue)
        .def("SetSolutionStepValue", &SetBufferValue)
        .def("GetSolutionStepValue", &GetValue)
        .def("GetSolutionStepValue", &GetBufferValue)
        //.def("__delitem__", &DataValueContainerDeleteValue)
        //.def("Erase", &DataValueContainerDeleteValue)
        ;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static void SetValue(TContainerType&  rData, TVariableType const& rV, typename TVariableType::Type const& rValue)
    {
        rData.GetSolutionStepValue(rV) = rValue;
    }

    static void SetBufferValue(TContainerType&  rData, TVariableType const& rV, typename TContainerType::IndexType SolutionStepIndex,typename TVariableType::Type const& rValue)
    {
        rData.GetSolutionStepValue(rV, SolutionStepIndex) = rValue;
    }

    static typename TVariableType::Type GetValue(TContainerType const& rData, TVariableType const& rV)
    {
        return rData.GetSolutionStepValue(rV);
    }

    inline
    static typename TVariableType::Type GetBufferValue(TContainerType const& rData, TVariableType const& rV, typename TContainerType::IndexType SolutionStepIndex)
    {
        return rData.GetSolutionStepValue(rV, SolutionStepIndex);
    }

    //static void DataValueContainerDeleteValue(TContainerType& rData, TVariableType const& rV)
    //{
    //  rData.Erase(rV);
    //}

    static bool DataValueContainerHas(TContainerType const& rData, TVariableType const& rV, typename TContainerType::IndexType SolutionStepIndex)
    {
        return rData.SolutionStepsDataHas(rV);
    }




    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SolutionStepVariableIndexingPython& operator=(const SolutionStepVariableIndexingPython& rOther);



    ///@}

}; // Class SolutionStepVariableIndexingPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_SOLUTION_STEP_VARIABLE_INDEXING_PYTHON_H_INCLUDED  defined 


