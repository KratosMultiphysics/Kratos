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


#if	!defined(KRATOS_ARRAY_1D_H_INCLUDED	)
#define	 KRATOS_ARRAY_1D_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <array>


// External	includes

// Project includes
#include "includes/define.h"
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
#include "includes/amatrix_interface.h"
#else
#include "includes/ublas_interface.h"

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/detail/vector_assign.hpp>
#endif // ifdef KRATOS_USE_AMATRIX

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

#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
template <typename TDataType, std::size_t TSize> using array_1d = Internals::Matrix<TDataType,TSize, 1>;
#else
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

    typedef std::size_t size_type;
    typedef	std::ptrdiff_t difference_type;
    typedef	T value_type;
    typedef	typename boost::numeric::ublas::type_traits<T>::const_reference const_reference;
    typedef	T &reference;
    typedef	std::array<T,N> array_type;
    typedef	T *pointer;
    typedef	array_1d<T, N> self_type;
    typedef	const boost::numeric::ublas::vector_reference<const self_type>	const_closure_type;
    typedef	boost::numeric::ublas::vector_reference<self_type>	closure_type;
    typedef	self_type vector_temporary_type;
    typedef	boost::numeric::ublas::dense_tag storage_category;
//		typedef	concrete_tag simd_category; //removed for the new ublas

    ///@}
    ///@name Life	Cycle
    ///@{

    /// Default constructor.
    BOOST_UBLAS_INLINE
    array_1d ():
        vector_expression<self_type> ()
    {
        //sin esto no funciona en windows!!
        //std::fill (data().begin(), data().end(), value_type	());
    }
    explicit BOOST_UBLAS_INLINE
    array_1d (size_type array_size):
        vector_expression<self_type> ()
    {
        //sin esto no funciona en windows!!
        //std::fill (data().begin(), data().end(), value_type	());

        /* 				std::fill (data().begin(), data().end(), value_type	()); */
    }

    explicit BOOST_UBLAS_INLINE
    array_1d (size_type array_size, value_type v):
        vector_expression<self_type> ()
    {
        std::fill (data().begin(), data().begin() + array_size, v);
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
#ifndef NDEBUG
        if(i>=N) KRATOS_THROW_ERROR(std::argument_error,"index greater than the size of the array - index is i = ", i);
#endif
        return data_[i];
    }
    BOOST_UBLAS_INLINE
    reference operator () (size_type i)
    {
#ifndef NDEBUG
        if(i>=N) KRATOS_THROW_ERROR(std::argument_error,"index greater than the size of the array - index is i = ", i);
#endif
        return data_[i];
    }

    BOOST_UBLAS_INLINE
    const_reference	operator []	(size_type i) const
    {
#ifndef NDEBUG
        if(i>=N) KRATOS_THROW_ERROR(std::argument_error,"index greater than the size of the array - index is i = ", i);
#endif
        return data_[i];
    }
    BOOST_UBLAS_INLINE
    reference operator [] (size_type i)
    {
#ifndef NDEBUG
        if(i>=N) KRATOS_THROW_ERROR(std::argument_error,"index greater than the size of the array - index is i = ", i);
#endif
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
    typedef	typename array_type::const_iterator const_iterator_type;
    typedef	typename array_type::iterator iterator_type;

public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
    typedef	indexed_iterator<self_type,	dense_random_access_iterator_tag> iterator;
    typedef	indexed_const_iterator<self_type, dense_random_access_iterator_tag>	const_iterator;
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
        typedef	dense_random_access_iterator_tag iterator_category;
#ifdef BOOST_MSVC_STD_ITERATOR
        typedef	const_reference	reference;
#else
        typedef	typename array_1d::difference_type difference_type;
        typedef	typename array_1d::value_type	value_type;
        typedef	typename array_1d::const_reference reference;
        typedef	const typename array_1d::pointer pointer;
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
        typedef	dense_random_access_iterator_tag iterator_category;
#ifndef	BOOST_MSVC_STD_ITERATOR
        typedef	typename array_1d::difference_type difference_type;
        typedef	typename array_1d::value_type	value_type;
        typedef	typename array_1d::reference reference;
        typedef	typename array_1d::pointer pointer;
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
    typedef	reverse_iterator_base<const_iterator, value_type, const_reference> const_reverse_iterator;
#else
    typedef	boost::numeric::ublas::reverse_iterator_base<const_iterator> const_reverse_iterator;
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
    typedef	reverse_iterator_base<iterator,	value_type,	reference> reverse_iterator;
#else
    typedef	boost::numeric::ublas::reverse_iterator_base<iterator>	reverse_iterator;
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
#endif // ifndef KRATOS_USE_AMATRIX

///@}

///@name Type	Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace	Kratos.

#endif // KRATOS_ARRAY_1D_H_INCLUDED  defined
