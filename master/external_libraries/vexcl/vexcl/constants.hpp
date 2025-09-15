#ifndef VEXCL_CONSTANTS_HPP
#define VEXCL_CONSTANTS_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   vexcl/constants.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Constants for use in vector expressions.
 */

#include <sstream>
#include <iomanip>
#include <string>
#include <type_traits>

#include <vexcl/backend.hpp>

#include <vexcl/types.hpp>
#include <vexcl/operations.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/version.hpp>

namespace vex {

//---------------------------------------------------------------------------
// std::integral_constant
//---------------------------------------------------------------------------
template <class T, T v>
struct is_cl_native< std::integral_constant<T, v> > : std::true_type {};

namespace traits {

template <class T, T v>
struct kernel_param_declaration< std::integral_constant<T, v> >
{
    static void get(backend::source_generator&,
            const std::integral_constant<T, v>&,
            const backend::command_queue&, const std::string &/*prm_name*/,
            detail::kernel_generator_state_ptr)
    { }
};

template <class T, T v>
struct partial_vector_expr< std::integral_constant<T, v> >
{
    static void get(backend::source_generator &src,
            const std::integral_constant<T, v>&,
            const backend::command_queue&, const std::string &/*prm_name*/,
            detail::kernel_generator_state_ptr)
    {
        src << v;
    }
};

template <class T, T v>
struct kernel_arg_setter< std::integral_constant<T, v> >
{
    static void set(const std::integral_constant<T, v>&,
            backend::kernel&, unsigned/*part*/, size_t/*index_offset*/,
            detail::kernel_generator_state_ptr)
    {
    }
};

} // namespace traits

//---------------------------------------------------------------------------
// boost::math::constants wrappers
//---------------------------------------------------------------------------
template <class Impl>
struct user_constant {
    typedef typename Impl::value_type value_type;
};

template <class Impl>
struct is_cl_native< user_constant<Impl> > : std::true_type {};

namespace traits {

template <class Impl>
struct kernel_param_declaration< user_constant<Impl> >
{
    static void get(backend::source_generator&,
            const user_constant<Impl>&,
            const backend::command_queue&, const std::string &/*prm_name*/,
            detail::kernel_generator_state_ptr)
    { }
};

template <class Impl>
struct partial_vector_expr< user_constant<Impl> >
{
    static void get(backend::source_generator &src,
            const user_constant<Impl>&,
            const backend::command_queue&, const std::string &/*prm_name*/,
            detail::kernel_generator_state_ptr)
    {
        src << Impl::get();
    }
};

template <class Impl>
struct kernel_arg_setter< user_constant<Impl> >
{
    static void set(const user_constant<Impl>&,
            backend::kernel&, unsigned/*part*/, size_t/*index_offset*/,
            detail::kernel_generator_state_ptr)
    {
    }
};

} // namespace traits

/**
 * Creates user-defined constan functor for use in VexCL expressions.
 * ``value`` will be copied verbatim into kernel source.
 */
#define VEX_CONSTANT(name, value)                                              \
  struct constant_##name {                                                     \
    typedef decltype(value) value_type;                                        \
    static std::string get() {                                                 \
      static const value_type v = value;                                       \
      std::ostringstream s;                                                    \
      s << "( " << std::scientific << std::setprecision(16) << v << " )";      \
      return s.str();                                                          \
    }                                                                          \
    decltype(boost::proto::as_expr<vex::vector_domain>(                        \
          vex::user_constant<constant_##name>()))                              \
    operator()() const {                                                       \
      return boost::proto::as_expr<vex::vector_domain>(                        \
          vex::user_constant<constant_##name>());                              \
    }                                                                          \
    operator value_type() const {                                              \
      static const value_type v = value;                                       \
      return v;                                                                \
    }                                                                          \
  };                                                                           \
  const constant_##name name = {}

/// Mathematical constants.
namespace constants {

VEX_CONSTANT( pi, boost::math::constants::pi<double>() );
VEX_CONSTANT( root_pi, boost::math::constants::root_pi<double>() );
VEX_CONSTANT( root_half_pi, boost::math::constants::root_half_pi<double>() );
VEX_CONSTANT( root_two_pi, boost::math::constants::root_two_pi<double>() );
VEX_CONSTANT( root_ln_four, boost::math::constants::root_ln_four<double>() );
VEX_CONSTANT( e, boost::math::constants::e<double>() );
VEX_CONSTANT( half, boost::math::constants::half<double>() );
VEX_CONSTANT( euler, boost::math::constants::euler<double>() );
VEX_CONSTANT( root_two, boost::math::constants::root_two<double>() );
VEX_CONSTANT( ln_two, boost::math::constants::ln_two<double>() );
VEX_CONSTANT( ln_ln_two, boost::math::constants::ln_ln_two<double>() );
VEX_CONSTANT( third, boost::math::constants::third<double>() );
VEX_CONSTANT( twothirds, boost::math::constants::twothirds<double>() );
VEX_CONSTANT( pi_minus_three, boost::math::constants::pi_minus_three<double>() );
VEX_CONSTANT( four_minus_pi, boost::math::constants::four_minus_pi<double>() );
VEX_CONSTANT( two_pi, boost::math::constants::two_pi<double>() );
VEX_CONSTANT( half_root_two, boost::math::constants::half_root_two<double>() );
VEX_CONSTANT( exp_minus_half, boost::math::constants::exp_minus_half<double>() );
#if (BOOST_VERSION >= 105000) || defined(DOXYGEN)
VEX_CONSTANT( one_div_two_pi, boost::math::constants::one_div_two_pi<double>() );
VEX_CONSTANT( catalan, boost::math::constants::catalan<double>() );
VEX_CONSTANT( cbrt_pi, boost::math::constants::cbrt_pi<double>() );
VEX_CONSTANT( cosh_one, boost::math::constants::cosh_one<double>() );
VEX_CONSTANT( cos_one, boost::math::constants::cos_one<double>() );
VEX_CONSTANT( degree, boost::math::constants::degree<double>() );
VEX_CONSTANT( e_pow_pi, boost::math::constants::e_pow_pi<double>() );
VEX_CONSTANT( euler_sqr, boost::math::constants::euler_sqr<double>() );
VEX_CONSTANT( extreme_value_skewness, boost::math::constants::extreme_value_skewness<double>() );
VEX_CONSTANT( four_thirds_pi, boost::math::constants::four_thirds_pi<double>() );
VEX_CONSTANT( glaisher, boost::math::constants::glaisher<double>() );
VEX_CONSTANT( half_pi, boost::math::constants::half_pi<double>() );
VEX_CONSTANT( khinchin, boost::math::constants::khinchin<double>() );
VEX_CONSTANT( ln_phi, boost::math::constants::ln_phi<double>() );
VEX_CONSTANT( ln_ten, boost::math::constants::ln_ten<double>() );
VEX_CONSTANT( log10_e, boost::math::constants::log10_e<double>() );
VEX_CONSTANT( one_div_cbrt_pi, boost::math::constants::one_div_cbrt_pi<double>() );
VEX_CONSTANT( one_div_euler, boost::math::constants::one_div_euler<double>() );
VEX_CONSTANT( one_div_ln_phi, boost::math::constants::one_div_ln_phi<double>() );
VEX_CONSTANT( one_div_log10_e, boost::math::constants::one_div_log10_e<double>() );
VEX_CONSTANT( one_div_root_pi, boost::math::constants::one_div_root_pi<double>() );
VEX_CONSTANT( one_div_root_two, boost::math::constants::one_div_root_two<double>() );
VEX_CONSTANT( one_div_root_two_pi, boost::math::constants::one_div_root_two_pi<double>() );
VEX_CONSTANT( phi, boost::math::constants::phi<double>() );
VEX_CONSTANT( pi_cubed, boost::math::constants::pi_cubed<double>() );
VEX_CONSTANT( pi_pow_e, boost::math::constants::pi_pow_e<double>() );
VEX_CONSTANT( pi_sqr, boost::math::constants::pi_sqr<double>() );
VEX_CONSTANT( pi_sqr_div_six, boost::math::constants::pi_sqr_div_six<double>() );
VEX_CONSTANT( radian, boost::math::constants::radian<double>() );
VEX_CONSTANT( rayleigh_kurtosis, boost::math::constants::rayleigh_kurtosis<double>() );
VEX_CONSTANT( rayleigh_kurtosis_excess, boost::math::constants::rayleigh_kurtosis_excess<double>() );
VEX_CONSTANT( rayleigh_skewness, boost::math::constants::rayleigh_skewness<double>() );
VEX_CONSTANT( root_e, boost::math::constants::root_e<double>() );
VEX_CONSTANT( root_one_div_pi, boost::math::constants::root_one_div_pi<double>() );
VEX_CONSTANT( root_three, boost::math::constants::root_three<double>() );
VEX_CONSTANT( root_two_div_pi, boost::math::constants::root_two_div_pi<double>() );
VEX_CONSTANT( sinh_one, boost::math::constants::sinh_one<double>() );
VEX_CONSTANT( sin_one, boost::math::constants::sin_one<double>() );
VEX_CONSTANT( sixth_pi, boost::math::constants::sixth_pi<double>() );
VEX_CONSTANT( third_pi, boost::math::constants::third_pi<double>() );
VEX_CONSTANT( three_quarters, boost::math::constants::three_quarters<double>() );
VEX_CONSTANT( three_quarters_pi, boost::math::constants::three_quarters_pi<double>() );
VEX_CONSTANT( two_div_pi, boost::math::constants::two_div_pi<double>() );
VEX_CONSTANT( two_thirds, boost::math::constants::two_thirds<double>() );
VEX_CONSTANT( two_thirds_pi, boost::math::constants::two_thirds_pi<double>() );
VEX_CONSTANT( zeta_three, boost::math::constants::zeta_three<double>() );
VEX_CONSTANT( zeta_two, boost::math::constants::zeta_two<double>() );
#endif


} //namespace constants

} // namespace vex

#endif
