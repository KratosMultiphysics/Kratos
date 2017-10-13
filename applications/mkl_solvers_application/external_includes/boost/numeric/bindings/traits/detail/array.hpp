/*
 * 
 * Copyright (c) Kresimir Fresl 2003
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

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_ARRAY_HPP
#define BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_ARRAY_HPP

#include <boost/numeric/bindings/traits/vector_traits.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#include <boost/numeric/bindings/traits/detail/array_impl.hpp>

namespace boost { namespace numeric { namespace bindings { namespace traits {

  template <typename T>
  struct vector_traits<detail::array<T> > {
    typedef T value_type;
    typedef T* pointer; 

    static pointer storage (detail::array<T>& a) { return a.storage(); }
    static int size (detail::array<T>& a) { return a.size(); } 
  }; 
  

}}}}

#endif // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS 

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_ARRAY_HPP
