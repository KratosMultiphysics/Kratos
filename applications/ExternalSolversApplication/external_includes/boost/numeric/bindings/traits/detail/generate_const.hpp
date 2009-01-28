/*
 * 
 * Copyright (c) 2002, 2003 Kresimir Fresl, Toon Knapen and Karl Meerbergen
 *
 * Permission to copy, modify, use and distribute this software 
 * for any non-commercial or commercial purpose is granted provided 
 * that this license appear on all copies of the software source code.
 *
 * Authors assume no responsibility whatsoever for its use and makes 
 * no guarantees about its quality, correctness or reliability.
 *
 * First author acknowledges the support of the Faculty of Civil 
 * Engineering, University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_GENERATE_CONST_HPP
#define BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_GENERATE_CONST_HPP

#include <boost/numeric/bindings/traits/config.hpp> 

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS

namespace boost { namespace numeric { namespace bindings { namespace traits { namespace detail {
 
  /// Copy const from V to X if present

  template <typename V, typename X>
  struct generate_const {
     typedef X type; 
  };

  template <typename V, typename X>
  struct generate_const< const V, X > {
     typedef X const type; 
  };

  template <typename V, typename X>
  struct generate_const< V const, X const > {
     typedef X const type; 
  };

}}}}}

#endif // BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_GENERATE_CONST_HPP
