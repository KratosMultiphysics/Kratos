//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED )
#define  KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED

// System includes
#include <vector>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/global_pointer.h"
#include "includes/serializer.h"
#include <boost/iterator/indirect_iterator.hpp>


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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
template< class TDataType >
class GlobalPointersVector : public std::vector< GlobalPointer<TDataType> >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GlobalPointersVector
    KRATOS_CLASS_POINTER_DEFINITION(GlobalPointersVector);

    typedef std::vector< GlobalPointer<TDataType> > TContainerType;
    typedef GlobalPointer<TDataType> TPointerType;
    typedef TDataType data_type;
    typedef TPointerType value_type;
    typedef TPointerType pointer;
    typedef const TPointerType const_pointer;
    typedef TDataType& reference;
    typedef const TDataType& const_reference;
    typedef TContainerType ContainerType;
    

    typedef boost::indirect_iterator<typename TContainerType::iterator>                iterator;
    typedef boost::indirect_iterator<typename TContainerType::const_iterator>          const_iterator;
    typedef boost::indirect_iterator<typename TContainerType::reverse_iterator>        reverse_iterator;
    typedef boost::indirect_iterator<typename TContainerType::const_reverse_iterator>  const_reverse_iterator;

    typedef typename TContainerType::size_type size_type;
    typedef typename TContainerType::iterator ptr_iterator;
    typedef typename TContainerType::const_iterator ptr_const_iterator;
    typedef typename TContainerType::reverse_iterator ptr_reverse_iterator;
    typedef typename TContainerType::const_reverse_iterator ptr_const_reverse_iterator;
    typedef typename TContainerType::difference_type difference_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GlobalPointersVector() {}

    GlobalPointersVector(const std::initializer_list<GlobalPointer<TDataType>>& l) 
        : std::vector< GlobalPointer<TDataType> >(l)
    {} 

    /// Destructor.
    virtual ~GlobalPointersVector() {}

    template < class TContainerType >
    void FillFromContainer( TContainerType& container)
    {
        this->reserve(container.size());
        for(auto& item : container)
        {
            this->push_back(GlobalPointer<TDataType>(&item));
        }
    }

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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GlobalPointersVector" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GlobalPointersVector";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    TDataType& operator[](const size_type& i)
    {
        return *(static_cast<TContainerType>(*this))[i];
    }

    TDataType const& operator[](const size_type& i) const
    {
        return *(static_cast<TContainerType>(*this))[i];
    }

    pointer& operator()(const size_type& i)
    {
        return (static_cast<TContainerType>(*this))[i];
    }

    const_pointer& operator()(const size_type& i) const
    {
        return (static_cast<TContainerType>(*this))[i];
    }

    bool operator==( const GlobalPointersVector<TDataType>& r ) const // nothrow
    {
        if( this->size() != r.size() )
            return false;
        else
            return std::equal(this->begin(), this->end(), r.begin(), this->EqualKeyTo());
    }

    ///@}
    ///@name Operations
    ///@{

    iterator                   begin()
    {
        return iterator( (static_cast<TContainerType>(*this)).begin() );
    }
    const_iterator             begin() const
    {
        return const_iterator( (static_cast<TContainerType>(*this)).begin() );
    }
    iterator                   end()
    {
        return iterator( (static_cast<TContainerType>(*this)).end() );
    }
    const_iterator             end() const
    {
        return const_iterator( (static_cast<TContainerType>(*this)).end() );
    }
    reverse_iterator           rbegin()
    {
        return reverse_iterator( (static_cast<TContainerType>(*this)).rbegin() );
    }
    const_reverse_iterator     rbegin() const
    {
        return const_reverse_iterator( (static_cast<TContainerType>(*this)).rbegin() );
    }
    reverse_iterator           rend()
    {
        return reverse_iterator( (static_cast<TContainerType>(*this)).rend() );
    }
    const_reverse_iterator     rend() const
    {
        return const_reverse_iterator( (static_cast<TContainerType>(*this)).rend() );
    }
    ptr_iterator               ptr_begin()
    {
        return (static_cast<TContainerType>(*this)).begin();
    }
    ptr_const_iterator         ptr_begin() const
    {
        return (static_cast<TContainerType>(*this)).cbegin();
    }
    ptr_iterator               ptr_end()
    {
        return (static_cast<TContainerType>(*this)).end();
    }
    ptr_const_iterator         ptr_end() const
    {
        return (static_cast<TContainerType>(*this)).cend();
    }
    ptr_reverse_iterator       ptr_rbegin()
    {
        return (static_cast<TContainerType>(*this)).rbegin();
    }
    ptr_const_reverse_iterator ptr_rbegin() const
    {
        return (static_cast<TContainerType>(*this)).rbegin();
    }
    ptr_reverse_iterator       ptr_rend()
    {
        return (static_cast<TContainerType>(*this)).rend();
    }
    ptr_const_reverse_iterator ptr_rend() const
    {
        return (static_cast<TContainerType>(*this)).rend();
    }

    reference        front()       /* nothrow */
    {
        assert( !this->empty() );
        return *((static_cast<TContainerType>(*this)).front());
    }
    const_reference  front() const /* nothrow */
    {
        assert( !this->empty() );
        return *((static_cast<TContainerType>(*this)).front());
    }
    reference        back()        /* nothrow */
    {
        assert( !this->empty() );
        return *((static_cast<TContainerType>(*this)).back());
    }
    const_reference  back() const  /* nothrow */
    {
        assert( !this->empty() );
        return *((static_cast<TContainerType>(*this)).back());
    }


    iterator insert(iterator Position, const TPointerType pData)
    {
        return iterator((static_cast<TContainerType>(*this)).insert(Position, pData));
    }

    template <class InputIterator>
    void insert(InputIterator First, InputIterator Last)
    {
        for(; First != Last; ++First)
            insert(*First);
    }


    iterator erase(iterator pos)
    {
        return iterator((static_cast<TContainerType>(*this)).erase(pos.base()));
    }

    iterator erase( iterator first, iterator last )
    {
        return iterator( (static_cast<TContainerType>(*this)).erase( first.base(), last.base() ) );
    }



    ///@}
    ///@name Access
    ///@{

    /** Gives a reference to underly normal container. */
    TContainerType& GetContainer()
    {
        return static_cast<TContainerType>(*this);
    }

    /** Gives a constant reference to underly normal container. */
    const TContainerType& GetContainer() const
    {
        return static_cast<TContainerType>(*this);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
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
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        if(rSerializer.Is(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION))
        {
            std::size_t pointer_size = sizeof(GlobalPointer<TDataType> );

            std::string data;
            data.resize( this->size() * pointer_size);
            for(std::size_t i=0; i<this->size(); ++i)
            {
                (*this)(i).save(&data[0]+i*pointer_size);
            }

            rSerializer.save("Size", this->size());
            rSerializer.save("Data", data);
        }
        else //SERIALIZING THE POINTER CONTENT TOO
        {
            rSerializer.save("Size", this->size());
            for(auto it = this->ptr_begin(); it!=this->ptr_end(); ++it)
            {
                TPointerType p = *it;
                rSerializer.save("Gp", p);
            }
        }
    }



    void load(Serializer& rSerializer)
    {
        if(rSerializer.Is(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION))
        {
            std::size_t pointer_size = sizeof(GlobalPointer<TDataType> );

            std::size_t size;
            rSerializer.load("Size", size);
            this->reserve(size);

            std::string tmp;
            rSerializer.load("Data", tmp);

            for(std::size_t i = 0; i<size; ++i)
            {
                GlobalPointer<TDataType> p(nullptr);
                p.load(&tmp[0]+i*pointer_size);
                this->push_back(p);
            }
        }
        else //SERIALIZING THE POINTER CONTENT TOO
        {
            std::size_t size;
            rSerializer.load("Size", size);
            this->reserve(size);
            for(std::size_t i = 0; i<size; ++i)
            {
                GlobalPointer<TDataType> p(nullptr);
                rSerializer.load("Gp", p);
                this->push_back(p);
            }
        }

    }


    ///@}
    ///@name Private Operations
    ///@{


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

}; // Class GlobalPointersVector

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDataType >
inline std::istream& operator >> (std::istream& rIStream,
                                  GlobalPointersVector<TDataType>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDataType >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GlobalPointersVector<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED  defined 


