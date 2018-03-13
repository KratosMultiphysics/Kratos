#ifndef boost_numeric_bindings_type_hpp
#define boost_numeric_bindings_type_hpp

// This header provides typedefs to float complex and double complex.
// This makes it possible to redefine the complex class being used.

#include <complex>

namespace boost { namespace numeric { namespace bindings { namespace traits {

  typedef std::complex< float >  complex_f ;
  typedef std::complex< double > complex_d ; 

  template <typename T> 
  T real (std::complex<T> const& c) { return std::real (c); }
  template <typename T> 
  T imag (std::complex<T> const& c) { return std::imag (c); }

}}}}

#endif // boost_numeric_bindings_type_hpp
