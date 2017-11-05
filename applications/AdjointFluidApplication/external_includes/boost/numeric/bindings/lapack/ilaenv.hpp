/*
 *
 * Copyright (c) Toon Knapen, Kresimir Fresl and Matthias Troyer 2003
 *
 * Permission to copy, modify, use and distribute this software
 * for any non-commercial or commercial purpose is granted provided
 * that this license appear on all copies of the software source code.
 *
 * Authors assume no responsibility whatsoever for its use and makes
 * no guarantees about its quality, correctness or reliability.
 *
 * KF acknowledges the support of the Faculty of Civil Engineering,
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_ILAENV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_ILAENV_HPP

#include <cstring> 
#include <boost/numeric/bindings/lapack/lapack.h>

namespace boost { namespace numeric { namespace bindings { namespace lapack {

  /*
   * ilaenv() is called from the LAPACK routines to choose  
   * problem-dependent parameters such as the block sizes 
   * for the local environment.
   */
  
  inline
  int ilaenv (int const ispec, const char* name, const char* opts,
              int const n1 = -1, int const n2 = -1, 
              int const n3 = -1, int const n4 = -1)
  {
    return ::LAPACK_ILAENV (&ispec, name, opts, &n1, &n2, &n3, &n4,
                          std::strlen (name), std::strlen (opts));
  }


}}}}

#endif

