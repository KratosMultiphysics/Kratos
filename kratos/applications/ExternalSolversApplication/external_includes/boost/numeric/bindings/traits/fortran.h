//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//  Copyright (C) 2002, 2003 Si-Lab b.v.b.a., Toon Knapen and Kresimir Fresl 

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_FORTRAN_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_FORTRAN_H

// First we need to know what the conventions for linking
// C with Fortran is on this platform/toolset
#if defined(__GNUC__) || defined(__ICC) || defined(__sgi) || defined(__COMO__) || defined(__KCC)
#define BIND_FORTRAN_LOWERCASE_UNDERSCORE
#elif defined(__IBMCPP__) || defined(_MSC_VER)
#define BIND_FORTRAN_LOWERCASE
#else
#error do not know how to link with fortran for the given platform
#endif

// Next we define macro's to convert our symbols to 
// the current convention
#if defined(BIND_FORTRAN_LOWERCASE_UNDERSCORE)
#define FORTRAN_ID( id ) id##_
#elif defined(BIND_FORTRAN_LOWERCASE)
#define FORTRAN_ID( id ) id
#else
#error do not know how to bind to fortran calling convention
#endif

#endif // BOOST_NUMERIC_BINDINGS_TRAITS_FORTRAN_H
