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
//


#if !defined(KRATOS_WEAK_POINTER_VECTOR_ITERATOR_H_INCLUDED )
#define  KRATOS_WEAK_POINTER_VECTOR_ITERATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/iterator/iterator_adaptor.hpp>


// Project includes
#include "includes/define.h"


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

/// Short class definition.
/** Detail class definition.
*/
template<class TIteratorType, class TDataType>
class WeakPointerVectorIterator
    : public boost::iterator_adaptor<WeakPointerVectorIterator<TIteratorType, TDataType>,
      TIteratorType, TDataType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WeakPointerVectorIterator
    KRATOS_CLASS_POINTER_DEFINITION(WeakPointerVectorIterator);

    typedef boost::iterator_adaptor<WeakPointerVectorIterator,
            TIteratorType, TDataType> BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    WeakPointerVectorIterator() {}

    WeakPointerVectorIterator(TIteratorType NewIterator) :BaseType(NewIterator) {}

    WeakPointerVectorIterator(WeakPointerVectorIterator const & NewIterator) :BaseType(NewIterator.base()) {}

    //template<class TOtherIteratorType>
    //WeakPointerVectorIterator(WeakPointerVectorIterator<TIteratorType> const & NewIterator) :BaseType(NewIterator.base()) {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
    friend class boost::iterator_core_access;
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

    typename BaseType::reference dereference() const
    {
        typename TDataType::Pointer p( (this->base())->lock() );
        return *p;
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


    ///@}

}; // Class WeakPointerVectorIterator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_WEAK_POINTER_VECTOR_ITERATOR_H_INCLUDED  defined 


