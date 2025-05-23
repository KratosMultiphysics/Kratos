//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <array>
#include <algorithm>
#include <initializer_list>

// External	includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/detail/vector_assign.hpp>

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type	Definitions
///@{

///@}
///@name	Enum's
///@{

///@}
///@name	Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short	class definition.
/** Detail class definition.
*/
template<class T,	std::size_t	N>
class	array_1d	: public boost::numeric::ublas::vector_expression< array_1d<T, N> >
{
public:
//#ifndef	BOOST_UBLAS_NO_PROXY_SHORTCUTS
//		BOOST_UBLAS_USING vector_expression<array_1d<T, N> >::operator ();
//#endif

    ///@name Type	Definitions
    ///@{

    /// Pointer definition of	array_1d
    KRATOS_CLASS_POINTER_DEFINITION(array_1d);

    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using const_reference = typename boost::numeric::ublas::type_traits<T>::const_reference;
    using reference = T&;
    using array_type = std::array<T,N>;
    using pointer = T*;
    using self_type = array_1d<T, N>;
    using const_closure_type = const boost::numeric::ublas::vector_reference<const self_type>;
    using closure_type = boost::numeric::ublas::vector_reference<self_type>;
    using vector_temporary_type = self_type;
    using storage_category = boost::numeric::ublas::dense_tag;
//		using simd_category = concrete_tag; //removed for the new ublas

    ///@}
    ///@name Life	Cycle
    ///@{

    /// Default constructor.
    BOOST_UBLAS_INLINE
    array_1d ():
        vector_expression<self_type> ()
    {
        // intentionally does not initialize the contents for performance reasons
    }

    explicit BOOST_UBLAS_INLINE
    array_1d (size_type array_size):
        vector_expression<self_type> ()
    {
        // intentionally does not initialize the contents for performance reasons
    }

    explicit BOOST_UBLAS_INLINE
    array_1d (size_type array_size, value_type v):
        vector_expression<self_type> ()
    {
        KRATOS_DEBUG_ERROR_IF(array_size>N) << "Given size is greater than the size of the array!" << std::endl;

        std::fill(data().begin(), data().begin() + array_size, v);
        // intentionally does not initialize the remaining entries for performance reasons
    }

    explicit BOOST_UBLAS_INLINE
    array_1d (const std::initializer_list<value_type>& rInitList):
        vector_expression<self_type> ()
    {
        KRATOS_DEBUG_ERROR_IF(rInitList.size()>N) << "Size of list greater than the size of the array!" << std::endl;

        std::copy(rInitList.begin(), rInitList.end(), data().begin()); // copy content of initializer list
        // intentionally does not initialize the remaining entries for performance reasons
    }

    BOOST_UBLAS_INLINE
    array_1d (size_type array_size,	const array_type & rdata):
        vector_expression<self_type> (),
        data_ (rdata) {}

    BOOST_UBLAS_INLINE
    array_1d (const array_1d &v):
        vector_expression<self_type> (),
        data_ (v.data_)	{}

//		template<class AE>
//		BOOST_UBLAS_INLINE
//		array_1d (const vector_expression<AE>	&ae):
//			vector_expression<self_type> ()	{
//			vector_assign (scalar_assign<reference,	typename AE::value_type> (), *this,	ae);
//		template<class AE> //boost 1.33.1
//		BOOST_UBLAS_INLINE
//		array_1d (const vector_expression<AE> &ae):
//            vector_expression<self_type> () {
//            vector_assign<scalar_assign> (*this, ae);

    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d (const boost::numeric::ublas::vector_expression<AE> &ae)
    {
        boost::numeric::ublas::vector_assign<boost::numeric::ublas::scalar_assign> (*this, ae);
    }

    ///@}
    ///@name Operators
    ///@{

    // Element access
    BOOST_UBLAS_INLINE
    const_reference	operator ()	(size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i>=N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return data_[i];
    }
    BOOST_UBLAS_INLINE
    reference operator () (size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i>=N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return data_[i];
    }

    BOOST_UBLAS_INLINE
    const_reference	operator []	(size_type i) const
    {
        KRATOS_DEBUG_ERROR_IF(i>=N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return data_[i];
    }
    BOOST_UBLAS_INLINE
    reference operator [] (size_type i)
    {
        KRATOS_DEBUG_ERROR_IF(i>=N) << "Index greater than the size of the array - index is i = " << i << std::endl;
        return data_[i];
    }

    // Assignment
    BOOST_UBLAS_INLINE
    array_1d &operator = (const array_1d &v)
    {
        data_ 	= v.data_;
        return *this;
    }

    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &operator = (const boost::numeric::ublas::vector_expression<AE>	&ae)
    {
        return assign (self_type	(ae));
        //self_type temporary	(ae);
        //return assign_temporary	(temporary);
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &operator +=	(const boost::numeric::ublas::vector_expression<AE> &ae)
    {
        return assign (self_type	(*this + ae));
        //self_type temporary	(*this + ae);
        //return assign_temporary	(temporary);
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &operator -=	(const boost::numeric::ublas::vector_expression<AE> &ae)
    {
        return assign (self_type	(*this - ae));
        //self_type temporary	(*this - ae);
        //return assign_temporary	(temporary);
    }
    template<class AT>
    BOOST_UBLAS_INLINE
    array_1d &operator /=	(const AT &at)
    {
        vector_assign_scalar<scalar_divides_assign> (*this, at); //included for ublas 1.33.1
        return *this;
    }

    /**
     * @brief Compares whether this array_1d is equal to the given array_1d.
     * @param v the array_1d to compare to
     * @return true if the two arrays are equal, false otherwise
     */
    BOOST_UBLAS_INLINE
    bool operator == (const array_1d &v) const
    {
        return std::equal (data_.begin(), data_.end(), v.data_.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    // Resizing
    BOOST_UBLAS_INLINE
    void resize	(size_type array_size, bool preserve = true)
    {
        if (!preserve)
            std::fill (data_.begin(), data_.end(), value_type	());
    }

    BOOST_UBLAS_INLINE
    array_1d &assign_temporary (array_1d &v)
    {
        swap (v);
        return *this;
    }


    template<class AT>
    BOOST_UBLAS_INLINE
    array_1d &operator *=	(const AT &at)
    {
        vector_assign_scalar<scalar_multiplies_assign> (*this, at); //included for ublas 1.33.1
        return *this;
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &plus_assign	(const boost::numeric::ublas::vector_expression<AE> &ae)
    {
        vector_assign<scalar_plus_assign> (*this, ae); //included for ublas 1.33.1
        //vector_assign (scalar_plus_assign<reference, typename AE::value_type> (), *this, ae);
        return *this;
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &assign (const boost::numeric::ublas::vector_expression<AE>	&ae)
    {
        vector_assign<scalar_assign> (*this, ae); //included for ublas 1.33.1
        //vector_assign (scalar_assign<reference,	typename AE::value_type> (), *this,	ae);
        return *this;
    }
    // Swapping
    BOOST_UBLAS_INLINE
    void swap (array_1d &v)
    {
        if (this !=	&v)
        {
            data ().swap (v.data ());
        }
    }
#ifndef	BOOST_UBLAS_NO_MEMBER_FRIENDS
    BOOST_UBLAS_INLINE
    friend void	swap (array_1d &v1, array_1d &v2)
    {
        v1.swap	(v2);
    }
#endif

    // Element insertion and erasure
    // These functions should work with	std::vector.
    // Thanks to Kresimir Fresl	for	spotting this.
//		BOOST_UBLAS_INLINE
//		void insert	(size_type i, const_reference t) {
    // FIXME: only works for EqualityComparable	value types.
    // BOOST_UBLAS_CHECK (data () [i] == value_type	(0), bad_index ());
    // Previously: data	().insert (data	().begin ()	+ i, t);
//			data ()	[i]	= t;
//		}
//		BOOST_UBLAS_INLINE
//		void erase (size_type i) {
//			// Previously: data	().erase (data ().begin	() + i);
//			data ()	[i]	= value_type (0);
//		}
    // Element assignment
    BOOST_UBLAS_INLINE
    reference insert_element (size_type i, const_reference t)
    {
        BOOST_UBLAS_CHECK (i < N, bad_index ());
        return (data_ [i] = t);
    }
    BOOST_UBLAS_INLINE
    void erase_element (size_type i)
    {
        BOOST_UBLAS_CHECK (i < N, bad_index ());
        data_ [i] = value_type/*zero*/();
    }
    BOOST_UBLAS_INLINE
    void clear ()
    {
        // Previously: data	().clear ();
        std::fill (data	().begin (), data ().end (), value_type	(0));
    }

    ///@}
    ///@name Access
    ///@{

    // Iterator	types
private:
    // Use the storage array1 iterator
    using const_iterator_type = typename array_type::const_iterator;
    using iterator_type = typename array_type::iterator;

public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
    using iterator = indexed_iterator<self_type,	dense_random_access_iterator_tag>;
    using const_iterator = indexed_const_iterator<self_type, dense_random_access_iterator_tag>;
#else
    class const_iterator;
    class iterator;
#endif

    // Element lookup
    BOOST_UBLAS_INLINE
    const_iterator find	(size_type i) const
    {
#ifndef	BOOST_UBLAS_USE_INDEXED_ITERATOR
        return const_iterator (*this, data ().begin	() + i);
#else
        return const_iterator (*this, i);
#endif
    }
    BOOST_UBLAS_INLINE
    iterator find (size_type i)
    {
#ifndef	BOOST_UBLAS_USE_INDEXED_ITERATOR
        return iterator	(*this,	data ().begin () + i);
#else
        return iterator	(*this,	i);
#endif
    }
    BOOST_UBLAS_INLINE
    size_type size () const
    {
        return N;
    }
    template<class AE>
    BOOST_UBLAS_INLINE
    array_1d &minus_assign (const	boost::numeric::ublas::vector_expression<AE> &ae)
    {
        vector_assign<scalar_minus_assign>(*this,ae);
        //vector_assign (scalar_minus_assign<reference, typename AE::value_type> (), *this, ae);
        return *this;
    }
    BOOST_UBLAS_INLINE
    const array_type &data () const
    {
        return data_;
    }
    BOOST_UBLAS_INLINE
    array_type &data ()
    {
        return data_;
    }

#ifndef	BOOST_UBLAS_USE_INDEXED_ITERATOR
    class const_iterator:
        public container_const_reference<array_1d>,
        public random_access_iterator_base<dense_random_access_iterator_tag,
        const_iterator, value_type, difference_type>
    {
    public:
        using iterator_category = dense_random_access_iterator_tag;
#ifdef BOOST_MSVC_STD_ITERATOR
        using reference = const_reference;
#else
        using difference_type = typename array_1d::difference_type;
        using value_type = typename array_1d::value_type;
        using reference = typename array_1d::const_reference;
        using pointer = const typename array_1d::pointer;
#endif

        // Construction	and	destruction
        BOOST_UBLAS_INLINE
        const_iterator ():
            container_const_reference<self_type> (), it_ ()	{}
        BOOST_UBLAS_INLINE
        const_iterator (const self_type	&v,	const const_iterator_type &it):
            container_const_reference<self_type> (v), it_ (it) {}
        BOOST_UBLAS_INLINE
#ifndef	BOOST_UBLAS_QUALIFIED_TYPENAME
        const_iterator (const iterator &it):
#else
        const_iterator (const typename self_type::iterator &it):
#endif
            container_const_reference<self_type> (it ()), it_ (it.it_) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        const_iterator &operator ++	()
        {
            ++ it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator &operator --	()
        {
            -- it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator &operator +=	(difference_type n)
        {
            it_	+= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        const_iterator &operator -=	(difference_type n)
        {
            it_	-= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type	operator - (const const_iterator &it) const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ - it.it_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        const_reference	operator * () const
        {
            BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index	());
            return *it_;
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index	() const
        {
            BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_, bad_index	());
            return it_ - (*this) ().begin ().it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        const_iterator &operator = (const const_iterator &it)
        {
            container_const_reference<self_type>::assign (&it ());
            it_	= it.it_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const	const_iterator &it)	const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ == it.it_;
        }
        BOOST_UBLAS_INLINE
        bool operator <	(const const_iterator &it) const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ < it.it_;
        }

    private:
        const_iterator_type	it_;

        friend class iterator;
    };
#endif

#ifndef	BOOST_UBLAS_USE_INDEXED_ITERATOR
    class iterator:
        public container_reference<array_1d>,
        public random_access_iterator_base<dense_random_access_iterator_tag,
        iterator, value_type, difference_type>
    {
    public:
        using iterator_category = dense_random_access_iterator_tag;
#ifndef	BOOST_MSVC_STD_ITERATOR
        using difference_type = typename array_1d::difference_type;
        using value_type = typename array_1d::value_type;
        using reference = typename array_1d::reference;
        using pointer = typename array_1d::pointer;
#endif

        // Construction	and	destruction
        BOOST_UBLAS_INLINE
        iterator ():
            container_reference<self_type> (), it_ () {}
        BOOST_UBLAS_INLINE
        iterator (self_type	&v,	const iterator_type	&it):
            container_reference<self_type> (v),	it_	(it) {}

        // Arithmetic
        BOOST_UBLAS_INLINE
        iterator &operator ++ ()
        {
            ++ it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        iterator &operator -- ()
        {
            -- it_;
            return *this;
        }
        BOOST_UBLAS_INLINE
        iterator &operator += (difference_type n)
        {
            it_	+= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        iterator &operator -= (difference_type n)
        {
            it_	-= n;
            return *this;
        }
        BOOST_UBLAS_INLINE
        difference_type	operator - (const iterator &it)	const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ - it.it_;
        }

        // Dereference
        BOOST_UBLAS_INLINE
        reference operator * ()	const
        {
            BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_ , bad_index ());
            return *it_;
        }

        // Index
        BOOST_UBLAS_INLINE
        size_type index	() const
        {
            BOOST_UBLAS_CHECK (it_ >= (*this) ().begin ().it_ && it_ < (*this) ().end ().it_ , bad_index ());
            return it_ - (*this) ().begin ().it_;
        }

        // Assignment
        BOOST_UBLAS_INLINE
        iterator &operator = (const	iterator &it)
        {
            container_reference<self_type>::assign (&it	());
            it_	= it.it_;
            return *this;
        }

        // Comparison
        BOOST_UBLAS_INLINE
        bool operator == (const	iterator &it) const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ == it.it_;
        }
        BOOST_UBLAS_INLINE
        bool operator <	(const iterator	&it) const
        {
            BOOST_UBLAS_CHECK (&(*this)	() == &it (), external_logic ());
            return it_ < it.it_;
        }

    private:
        iterator_type it_;

        friend class const_iterator;
    };
#endif


    BOOST_UBLAS_INLINE
    const_iterator begin ()	const
    {
        return find	(0);
    }
    BOOST_UBLAS_INLINE
    const_iterator end () const
    {
        return find	(data_.size	());
    }

    BOOST_UBLAS_INLINE
    iterator begin ()
    {
        return find	(0);
    }
    BOOST_UBLAS_INLINE
    iterator end ()
    {
        return find	(data_.size	());
    }

    // Reverse iterator

#ifdef BOOST_MSVC_STD_ITERATOR
    using const_reverse_iterator = reverse_iterator_base<const_iterator, value_type, const_reference>;
#else
    using const_reverse_iterator = boost::numeric::ublas::reverse_iterator_base<const_iterator>;
#endif

    BOOST_UBLAS_INLINE
    const_reverse_iterator rbegin () const
    {
        return const_reverse_iterator (end ());
    }
    BOOST_UBLAS_INLINE
    const_reverse_iterator rend	() const
    {
        return const_reverse_iterator (begin ());
    }

#ifdef BOOST_MSVC_STD_ITERATOR
    using reverse_iterator = reverse_iterator_base<iterator,	value_type,	reference>;
#else
    using reverse_iterator = boost::numeric::ublas::reverse_iterator_base<iterator>;
#endif

    BOOST_UBLAS_INLINE
    reverse_iterator rbegin	()
    {
        return reverse_iterator	(end ());
    }
    BOOST_UBLAS_INLINE
    reverse_iterator rend ()
    {
        return reverse_iterator	(begin ());
    }
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

protected:
    ///@name Protected static	Member Variables
    ///@{


    ///@}
    ///@name Protected member	Variables
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

    array_type data_;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private	Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class	array_1d

///@}
///@name Type	Definitions
///@{

///@}
///@name Input and output
///@{

///@}


}  // namespace	Kratos.

namespace AuxiliaryHashCombine
{
    /**
     * @brief This method creates an "unique" hash for the input value
     * @details It comes from boost, taken from here: https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
     * @tparam TClassType The type of class to be hashed
     * @param rSeed This is the seed used to create the hash
     * @param rValue This is the value to be hashed
     * @todo Once the hashers and comparors are moved from key_hash.h, include key_hash and remove this. Right now there is a cross inclussion
     */
    template <class TClassType>
    inline void HashCombine(
        std::size_t& rSeed,
        const TClassType& rValue
        )
    {
        std::hash<TClassType> hasher;
        rSeed ^= hasher(rValue) + 0x9e3779b9 + (rSeed<<6) + (rSeed>>2);
    }
} /// namespace

namespace std
{
template<class T, std::size_t N>
struct hash<Kratos::array_1d<T,N>>
{
    std::size_t operator()(const Kratos::array_1d<T,N>& rArray) {
            std::size_t seed = 0;
            for (auto component : rArray) {AuxiliaryHashCombine::HashCombine(seed, component);}
            return seed;
        }
};
} // namespace std.
