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

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_UBLAS_UPLO_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_UBLAS_UPLO_H

#include <boost/numeric/ublas/fwd.hpp> 

namespace boost { namespace numeric { namespace bindings { namespace traits {

  namespace detail {

    template <typename UpLoTag>
    struct ublas_uplo {};
    
    template<> 
    struct ublas_uplo<boost::numeric::ublas::lower> {
      typedef lower_t type; 
    };
    template<> 
    struct ublas_uplo<boost::numeric::ublas::upper> {
      typedef upper_t type; 
    };

  }

}}}}

#endif 
