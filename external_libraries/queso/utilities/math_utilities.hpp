// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MATH_UTILITIES_HPP
#define MATH_UTILITIES_HPP

//// Project includes
#include "includes/define.hpp"

namespace queso {

namespace Math {
    ///@brief Simple Power function
    ///@param x Value
    ///@param p Order
    ///@details For gcc (without --ffast-math compiler flag) this is faster than std::pow().
    static double inline Pow( double x, std::size_t p){
        double result = 1.0;
        while( p > 0UL ) {
            result = result * x;
            p -= 1;
        }
        return result;
    }

    ///@brief Dot product of two vectors.
    ///@param Vector3d rLHs
    ///@param Vector3d rRHs
    ///@return double.
    static double inline Dot( const Vector3d& rLhs, const Vector3d& rRhs) {
        return (rLhs[0]*rRhs[0] + rLhs[1]*rRhs[1] + rLhs[2]*rRhs[2]);
    }

    ///@brief Cross product of two vectors.
    ///@param Vector3d rLHs
    ///@param Vector3d rRHs
    ///@return double.
    static Vector3d inline Cross( const Vector3d& rLhs, const Vector3d& rRhs) {
        return Vector3d( rLhs[1]*rRhs[2] - rLhs[2]*rRhs[1],
                         rLhs[2]*rRhs[0] - rLhs[0]*rRhs[2],
                         rLhs[0]*rRhs[1] - rLhs[1]*rRhs[0] );
    }

    ///@brief Norm of vector.
    ///@param Vector3d rLHs
    ///@return double.
    static double inline Norm( const Vector3d& rLhs ) {
        return std::sqrt( rLhs[0]*rLhs[0] + rLhs[1]*rLhs[1] + rLhs[2]*rLhs[2] );
    }

    /// @brief Returns max value of vector
    /// @tparam T
    /// @param rVector
    /// @return T
    template<typename T>
    static T inline Max( const Vector3<T>& rVector ){
        return std::max<T>(std::max<T>(rVector[0], rVector[1]), rVector[2]);
    }

    /// @brief Returns min value of vector.
    /// @tparam T
    /// @param rVector
    /// @return T
    template<typename T>
    static T inline Min( const Vector3<T>& rVector ){
        return std::min<T>(std::min<T>(rVector[0], rVector[1]), rVector[2]);
    }

} // End namespace math

} // End namespace queso

#endif