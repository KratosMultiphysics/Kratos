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


#if !defined(KRATOS_VECTOR_MAP_ITERATOR_H_INCLUDED )
#define  KRATOS_VECTOR_MAP_ITERATOR_H_INCLUDED



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
class VectorMapIterator
    : public boost::iterator_adaptor<VectorMapIterator<TIteratorType, TDataType>,
      TIteratorType, TDataType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VectorMapIterator
    KRATOS_CLASS_POINTER_DEFINITION(VectorMapIterator);

    typedef boost::iterator_adaptor<VectorMapIterator,
            TIteratorType,  TDataType> BaseType;

    typedef typename TIteratorType::value_type::first_type key_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VectorMapIterator() {}

    VectorMapIterator(TIteratorType NewIterator) :BaseType(NewIterator) {}

    VectorMapIterator(VectorMapIterator const & NewIterator) :BaseType(NewIterator.base()) {}

    //template<class TOtherIteratorType>
    //VectorMapIterator(VectorMapIterator<TIteratorType> const & NewIterator) :BaseType(NewIterator.base()) {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    key_type key()
    {
        return this->base()->first;
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
# if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x551))
        return const_cast<BaseType::reference>(*this->base()->second);
# else
        return this->base()->second;
# endif
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

}; // Class VectorMapIterator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_VECTOR_MAP_ITERATOR_H_INCLUDED  defined 


