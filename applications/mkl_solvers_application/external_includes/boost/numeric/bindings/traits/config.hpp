//
//  Copyright (c) 2002-2003
//  Toon Knapen, Kresimir Fresl, Joerg Walter
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//
//

#ifndef BOOST_NUMERIC_BINDINGS_CONFIG_HPP
#define BOOST_NUMERIC_BINDINGS_CONFIG_HPP

#include <boost/config.hpp> 

// Microsoft Visual C++
#if defined (BOOST_MSVC)
// .. version 6.0 & 7.0
#  if BOOST_MSVC <= 1300
#    define BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
#  endif

#elif defined(__ICC)

#elif defined(__IBMCPP__)

#elif defined(__GNUC__)

#elif defined(__COMO__)

#elif defined(__KCC)

#elif defined(__sgi)

#else
#error bindings do not recognise compiler
#endif
 

#if defined (BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS)

// structure checks require proper traits
#  define BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 

// type checks require proper traits 
#  define BOOST_NUMERIC_BINDINGS_NO_TYPE_CHECK 

#endif 

#endif

