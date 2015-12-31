/// Kratos Multi-Physics
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



#if !defined(KRATOS_VARIABLE_INDEXING_PYTHON_H_INCLUDED )
#define  KRATOS_VARIABLE_INDEXING_PYTHON_H_INCLUDED



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
class VariableIndexingPython : public def_visitor<VariableIndexingPython<TContainerType, TVariableType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VariableIndexingPython
    KRATOS_CLASS_POINTER_DEFINITION(VariableIndexingPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VariableIndexingPython() {}

    /// Copy constructor.
    VariableIndexingPython(const VariableIndexingPython& rOther);

    /// Destructor.
    virtual ~VariableIndexingPython() {}


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
        .def("__contains__", &VariableIndexingPython<TContainerType, TVariableType>::DataValueContainerHas)
        .def("__setitem__", &VariableIndexingPython<TContainerType, TVariableType>::DataValueContainerSetValue)
        .def("__getitem__", &VariableIndexingPython<TContainerType, TVariableType>::DataValueContainerGetValue)
        .def("Has", &VariableIndexingPython<TContainerType, TVariableType>::DataValueContainerHas)
        .def("SetValue", &VariableIndexingPython<TContainerType, TVariableType>::DataValueContainerSetValue)
        .def("GetValue", &VariableIndexingPython<TContainerType, TVariableType>::DataValueContainerGetValue)
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

    static void DataValueContainerSetValue(TContainerType&  rData, TVariableType const& rV, typename TVariableType::Type const& rValue)
    {
        rData.SetValue(rV, rValue);
    }

    static typename TVariableType::Type DataValueContainerGetValue(TContainerType const& rData, TVariableType const& rV)
    {
        return rData.GetValue(rV);
    }

    inline
    static typename TVariableType::Type const& DataValueContainerGetReference(TContainerType const& rData, TVariableType const& rV)
    {
        return rData.GetValue(rV);
    }

    //static void DataValueContainerDeleteValue(TContainerType& rData, TVariableType const& rV)
    //{
    //  rData.Erase(rV);
    //}

    static bool DataValueContainerHas(TContainerType const& rData, TVariableType const& rV)
    {
        return rData.Has(rV);
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
    VariableIndexingPython& operator=(const VariableIndexingPython& rOther);


    ///@}

}; // Class VariableIndexingPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_VARIABLE_INDEXING_PYTHON_H_INCLUDED  defined 


