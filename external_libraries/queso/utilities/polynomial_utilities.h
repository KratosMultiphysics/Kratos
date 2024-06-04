//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef POLYNOMIAL_UTILITIES_INCLUDE_H
#define POLYNOMIAL_UTILITIES_INCLUDE_H

// STL includes
#include <vector>
#include <variant>
#include <cmath>

namespace queso {

///@name QuESo Classes
///@{

///
/**
 * @class  Polynomial
 * @author Manuel Messmer
 * @brief  Provides evaluations of Legendre Polynomials.
 * @todo Check if this is slower or faster than a standard switch.
*/
class Polynomial {

public:
    ///@name Operations
    ///@{

    /// Legendre Polynomial defined on (a,b)
    static double inline f_x( double x, int order, double a, double b ){
        return f_x_visit( mLegendre[order], x, a, b);
    }

    /// Returns integral of Legendre polynomial defined on (a,b)
    static double inline f_x_int( double x, int order, double a, double b ){
        return f_x_int_visit( mLegendre[order], x, a, b);
    }

    ///@}
private:

    ///@name Type Definitions
    ///@{
    typedef struct {} p0;
    typedef struct {} p1;
    typedef struct {} p2;
    typedef struct {} p3;
    typedef struct {} p4;
    typedef struct {} p5;
    typedef struct {} p6;
    typedef struct {} p7;
    typedef struct {} p8;

    using Lp = std::variant<p0,p1,p2,p3,p4,p5,p6,p7,p8>;
    using Legendre = std::vector<Lp>;

    /// @brief F_x provides implementation of Legendre polynomials
    typedef struct F_x_ {
        // Shifted legendre polynomial (p=0)
        double inline operator()(const p0& p) const {
            return 1;
        }
        // Shifted legendre polynomial (p=1)
        double inline operator()(const p1& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return tmp_x;
        }
        // Shifted legendre polynomial (p=2)
        double inline operator()(const p2& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/2.0*(3.0*power(tmp_x,2)-1.0);
        }
        // Shifted legendre polynomial (p=3)
        double inline operator()(const p3& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/2.0*(5.0*power(tmp_x,3) - 3.0*tmp_x);
        }
        // Shifted legendre polynomial (p=4)
        double inline operator()(const p4& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/8.0*(35.0*power(tmp_x,4)-30.0*power(tmp_x,2) +3.0);
        }
        // Shifted legendre polynomial (p=5)
        double inline operator()(const p5& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/8.0*(63.0*power(tmp_x,5)-70.0*power(tmp_x,3)+15.0*tmp_x);
        }
        // Shifted legendre polynomial (p=6)
        double inline operator()(const p6& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/16.0*(231.0*power(tmp_x,6)-315.0*power(tmp_x,4)+105.0*power(tmp_x,2)-5.0);
        }
        // Shifted legendre polynomial (p=7)
        double inline operator()(const p7& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/16.0*(429.0*power(tmp_x,7)-693.0*power(tmp_x,5)+315.0*power(tmp_x,3)-35.0*tmp_x);
        }
        // Shifted legendre polynomial (p=8)
        double inline operator()(const p8& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/128.0*(6435.0*power(tmp_x,8) - 12012.0*power(tmp_x,6)+6930.0*power(tmp_x,4)-1260.0*power(tmp_x,2)+35.0);
        }
        double x{};
        double a{};
        double b{};
    } F_x;

    /// @brief F_x_int provides implementation of integrals of Legendre polynomials
    /// @todo: power() is still relatively slow.
    typedef struct F_x_int_ {
        // Integral of shifted legendre polynomial (p=0)
        double operator()(const p0& p) const {
            return x;
        }
        // Integral of shifted legendre polynomial (p=1)
        double operator()(const p1& p) const {
            return -power((a + b - 2.0*x),2)/(4.0*(a - b));;
        }
        // Integral of shifted legendre polynomial (p=2)
        double operator()(const p2& p) const {
            return -x/2.0 - power((a + b - 2.0*x),3)/(4.0*power((a - b),2));
        }
        // Integral of shifted legendre polynomial (p=3)
        double operator()(const p3& p) const {
            return (3.0*power( (a + b - 2.0*x),2) )/(8.0*(a - b)) - (5*power((a + b - 2.0*x),4))/(16*power((a - b),3));
        }
        // Integral of shifted legendre polynomial (p=4)
        double operator()(const p4& p) const {
            return (3.0*x)/8.0 + (5.0*power((a + b - 2.0*x),3))/(8*power((a - b),2)) - (7.0*power((a + b - 2*x),5))/(16*power((a - b),4));
        }
        // Integral of shifted legendre polynomial (p=5)
        double operator()(const p5& p) const {
            return (35*power((a + b - 2*x),4))/(32*power((a - b),3)) - (15*power((a + b - 2*x),2))/(32*(a - b)) - (21*power((a + b - 2*x),6))/(32*power((a - b),5));
        }
        // Integral of shifted legendre polynomial (p=6)
        double operator()(const p6& p) const {
            return (63*power((a + b - 2*x),5))/(32*power((a - b),4)) - (35*power((a + b - 2*x),3))/(32*power((a - b),2)) - (5*x)/16 - (33*power((a + b - 2*x),7))/(32*power((a - b),6));
        }
        // Integral of shifted legendre polynomial (p=7)
        double operator()(const p7& p) const {
            return (35.0*power( (a + b - 2*x),2))/(64*(a - b)) - (315*power((a + b - 2*x),4))/(128.0*power((a - b),3)) + (231.0*power((a + b - 2*x),6))/(64.0*power((a - b),5)) - (429.0*power((a + b - 2*x),8))/(256.0*power((a - b),7));
        }
        // Integral of shifted legendre polynomial (p=8)
        double operator()(const p8& p) const {
            return (35.0*x)/128.0 + (105.0*power( (a + b - 2*x),3) )/(64.0*power( (a - b),2) ) - (693.0*power( (a + b - 2*x),5) )/(128.0*power((a - b),4) ) + (429.0*power((a + b - 2*x),7))/(64.0*power((a - b),6)) - (715.0*power( (a + b - 2*x),9) )/(256.0*power((a - b),8));
        }
        double x{};
        double a{};
        double b{};
    } F_x_int;

    ///@}
    ///@name Private operations
    ///@{

    /// Caller for F_x
    static double inline f_x_visit( const Lp& s, double x, double a, double b ) {
        return std::visit( F_x{x, a, b}, s );
    }

    /// Caller for F_x_int
    static double inline f_x_int_visit( const Lp& s, double x, double a, double b ){
        return std::visit( F_x_int{x, a, b}, s );
    }

    ///@brief Simple Power functions
    ///@details For gcc (without --ffast-math compiler flag) this is faster than std::pow().
    static double inline power( double x, std::size_t p){
        double result = 1.0;
        while( p > 0UL ) {
            result = result * x;
            p -= 1;
        }
        return result;
    }
    ///@}
    ///@name Private operations
    ///@{

    inline static Legendre mLegendre = {Polynomial::p0{}, Polynomial::p1{}, Polynomial::p2{}, Polynomial::p3{}, Polynomial::p4{}, Polynomial::p5{}, Polynomial::p6{}, Polynomial::p7{}, Polynomial::p8{}};
    ///@}
}; // End Class

///@}

} // End namespace queso

#endif