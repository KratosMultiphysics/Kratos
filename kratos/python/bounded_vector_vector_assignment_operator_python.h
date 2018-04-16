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


#if !defined(KRATOS_BOUNDED_VECTOR_VECTOR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_BOUNDED_VECTOR_VECTOR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>


// External includes


// Project includes
#include "includes/define_python.h"


namespace Kratos
{

namespace Python
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

/// Short class definition.
/** Detail class definition.
*/
template<class TContainerType, class TOtherContainerType>
class BoundedVectorVectorAssignmentOperatorPython : public def_visitor<BoundedVectorVectorAssignmentOperatorPython<TContainerType, TOtherContainerType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BoundedVectorVectorAssignmentOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(BoundedVectorVectorAssignmentOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BoundedVectorVectorAssignmentOperatorPython() {}

    /// Copy constructor.
    BoundedVectorVectorAssignmentOperatorPython(const BoundedVectorVectorAssignmentOperatorPython& rOther);

    /// Destructor.
    virtual ~BoundedVectorVectorAssignmentOperatorPython() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template <class TClassType>
    void visit(TClassType& ThisClass) const
    {
        implicitly_convertible<TOtherContainerType, TContainerType>();
        ThisClass
        .def("__iadd__", &iadd, return_internal_reference<1>())
        .def("__isub__", &isub, return_internal_reference<1>())
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

    static
    TContainerType&
    iadd(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        if(ThisContainer.size()!=OtherContainer.size())
            KRATOS_THROW_ERROR(std::invalid_argument, "different size of vectors to sum", "");
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            ThisContainer[i] += OtherContainer[i];
        return ThisContainer;
    }

    static
    TContainerType&
    isub(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        if(ThisContainer.size()!=OtherContainer.size())
            KRATOS_THROW_ERROR(std::invalid_argument, "different size of vectors to subtract", "");
        for(typename TContainerType::size_type i = 0 ; i < ThisContainer.size() ; i++)
            ThisContainer[i] -= OtherContainer[i];
        return ThisContainer;
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
    BoundedVectorVectorAssignmentOperatorPython& operator=(const BoundedVectorVectorAssignmentOperatorPython& rOther);


    ///@}

}; // Class BoundedVectorVectorAssignmentOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_BOUNDED_VECTOR_VECTOR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED  defined 


