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




#if !defined(KRATOS_VECTOR_VECTOR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_VECTOR_VECTOR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED



// System includes
#include <string>
#include <iostream>


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
class VectorVectorAssignmentOperatorPython : public def_visitor<VectorVectorAssignmentOperatorPython<TContainerType, TOtherContainerType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VectorVectorAssignmentOperatorPython
    KRATOS_CLASS_POINTER_DEFINITION(VectorVectorAssignmentOperatorPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VectorVectorAssignmentOperatorPython() {}

    /// Copy constructor.
    VectorVectorAssignmentOperatorPython(const VectorVectorAssignmentOperatorPython& rOther);

    /// Destructor.
    virtual ~VectorVectorAssignmentOperatorPython() {}


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
        return ThisContainer += OtherContainer;
    }

    static
    TContainerType&
    isub(TContainerType& ThisContainer, TOtherContainerType const& OtherContainer)
    {
        return ThisContainer -= OtherContainer;
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
    VectorVectorAssignmentOperatorPython& operator=(const VectorVectorAssignmentOperatorPython& rOther);


    ///@}

}; // Class VectorVectorAssignmentOperatorPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_VECTOR_VECTOR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED  defined 


