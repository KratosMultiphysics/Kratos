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




#if !defined(KRATOS_SOLUTION_STEP_VARIABLE_INDEXING_PYTHON_H_INCLUDED )
#define  KRATOS_SOLUTION_STEP_VARIABLE_INDEXING_PYTHON_H_INCLUDED



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
class SolutionStepVariableIndexingPython //: public def_visitor<SolutionStepVariableIndexingPython<TContainerType, TVariableType> >
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


