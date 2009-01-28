/*
 * 
 * Copyright (c) Kresimir Fresl 2002 
 *
 * Permission to copy, modify, use and distribute this software 
 * for any non-commercial or commercial purpose is granted provided 
 * that this license appear on all copies of the software source code.
 *
 * Author assumes no responsibility whatsoever for its use and makes 
 * no guarantees about its quality, correctness or reliability.
 *
 * Author acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_SYMM_HERM_TRAITS_HPP
#define BOOST_NUMERIC_BINDINGS_TRAITS_SYMM_HERM_TRAITS_HPP

#include <boost/numeric/bindings/traits/type.hpp> 
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/static_assert.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS

namespace boost { namespace numeric { namespace bindings { namespace traits {

  namespace detail {

    // select symmetric or hermitian matrix structure 

    template <typename T> 
    struct symm_herm_t {}; 
    template<>
    struct symm_herm_t<float> {
      typedef symmetric_t type;
    }; 
    template<>
    struct symm_herm_t<double> {
      typedef symmetric_t type;
    }; 
    template<>
    struct symm_herm_t<complex_f> {
      typedef hermitian_t type;
    }; 
    template<>
    struct symm_herm_t<complex_d> {
      typedef hermitian_t type;
    }; 

    template <typename T> 
    struct symm_herm_pack_t {}; 
    template<>
    struct symm_herm_pack_t<float> {
      typedef symmetric_packed_t type;
    }; 
    template<>
    struct symm_herm_pack_t<double> {
      typedef symmetric_packed_t type;
    }; 
    template<>
    struct symm_herm_pack_t<complex_f> {
      typedef hermitian_packed_t type;
    }; 
    template<>
    struct symm_herm_pack_t<complex_d> {
      typedef hermitian_packed_t type;
    }; 


    template <class T, class S>
    struct symm_herm_compatible {
      BOOST_STATIC_CONSTANT( bool, value=false ) ;
    };

    template <class T>
    struct symm_herm_compatible< T, hermitian_t > {
      BOOST_STATIC_CONSTANT( bool, value=true ) ;
    };

    template <class T>
    struct symm_herm_compatible< T, symmetric_t > {
      BOOST_STATIC_CONSTANT( bool, value=true ) ;
    };

    template <class T>
    struct symm_herm_compatible< std::complex<T>, symmetric_t > {
      BOOST_STATIC_CONSTANT( bool, value=false ) ;
    };

  }

}}}}

#endif // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_SYMM_HERM_TRAITS_HPP 
