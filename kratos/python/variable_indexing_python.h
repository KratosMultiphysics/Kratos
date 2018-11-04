//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//




#if !defined(KRATOS_VARIABLE_INDEXING_PYTHON_H_INCLUDED )
#define  KRATOS_VARIABLE_INDEXING_PYTHON_H_INCLUDED



// System includes
#include <pybind11/pybind11.h>

// External includes


// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;
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
class VariableIndexingPython //: public def_visitor<VariableIndexingPython<TContainerType, TVariableType> >
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


