//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef DEFINE_CHECKS_HPP
#define DEFINE_CHECKS_HPP

//// Project includes
#include "queso/includes/define.hpp"

namespace queso {

#define QuESo_CHECK(IsTrue) if(!(IsTrue)) QuESo_ERROR << "Check failed because " << #IsTrue << " is not true\n";
#define QuESo_CHECK_IS_FALSE(IsFalse) if(IsFalse) QuESo_ERROR  << "Check failed because " << #IsFalse << " is not false\n";

#define QuESo_CHECK_EQUAL(a,b) if(!((a) == (b))) QuESo_ERROR << "Check failed because " << #a << " = " << a << " is not equal to " << #b << " = " << b << "\n";
#define QuESo_CHECK_NOT_EQUAL(a,b) if((a) == (b)) QuESo_ERROR << "Check failed because " << #a << " = " << a << " is equal to " << #b << " = " << b << "\n";

#define QuESo_CHECK_LT(a,b) if(!((a) < (b))) QuESo_ERROR << "Check failed because " << #a << " = " << a << " is not less than " << #b << " = " << b << "\n";
#define QuESo_CHECK_GT(a,b) if(!((a) > (b))) QuESo_ERROR << "Check failed because " << #a << " = " << a << " is not greater than " << #b << " = " << b << "\n";

#define QuESo_CHECK_NEAR(a, b, tolerance) if(!(std::abs(a - b) <= tolerance)) QuESo_ERROR << "Check failed because " << #a << " = " << a << \
" is not near to " << #b << " = " << b << " within the tolerance " << tolerance << "\n";

#define QuESo_CHECK_RELATIVE_NEAR(a, b, tolerance) if(!(std::abs(b) <= std::numeric_limits<double>::epsilon())) { QuESo_ERROR_IF(!(std::abs((a - b)/b) <= tolerance)) << "Check failed because " << #a << " = " << a << \
" is not near to " << #b << " = " << b << " within the relative tolerance " << tolerance << std::endl; } else {QuESo_CHECK_NEAR(a,b,tolerance);}

#define QuESo_CHECK_POINT_NEAR(a, b, tolerance)                          \
    QuESo_CHECK_NEAR(a[0], b[0], tolerance );                            \
    QuESo_CHECK_NEAR(a[1], b[1], tolerance );                            \
    QuESo_CHECK_NEAR(a[2], b[2], tolerance );

#define QuESo_CHECK_POINT_RELATIVE_NEAR(a, b, tolerance)                 \
    QuESo_CHECK_RELATIVE_NEAR(a[0], b[0], tolerance );                   \
    QuESo_CHECK_RELATIVE_NEAR(a[1], b[1], tolerance );                   \
    QuESo_CHECK_RELATIVE_NEAR(a[2], b[2], tolerance );

#define QuESo_CHECK_Vector3i_EQUAL(a, b)                                 \
    QuESo_CHECK_EQUAL(a[0], b[0] );                                      \
    QuESo_CHECK_EQUAL(a[1], b[1] );                                      \
    QuESo_CHECK_EQUAL(a[2], b[2] );

} // End namespace queso

#endif // DEFINE_CHECKS_HPP