// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MATH_UTILITIES_HPP
#define MATH_UTILITIES_HPP

/// STL includes
#include <cmath>
//// Project includes
#include "queso/includes/define.hpp"

namespace queso {

namespace Math {

    ///@brief Simple Power function
    ///@param x Value
    ///@param p Order
    ///@details For gcc (without --ffast-math compiler flag) this is faster than std::pow().
    static inline double Pow( double x, std::size_t p){
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

    ///@brief Returns: rLHs + rRhs
    ///@param Vector3d rLHs
    ///@param Vector3d rRHs
    ///@return Vector3d.
    static inline Vector3d Add( const Vector3d& rLhs, const Vector3d& rRhs) {
        return {rLhs[0]+rRhs[0], rLhs[1]+rRhs[1], rLhs[2]+rRhs[2]};
    }

    /// @brief Performs: rLHs += rRhs
    /// @param rLhs
    /// @param rRhs
    static inline void AddSelf( Vector3d& rLhs, const Vector3d& rRhs) {
        rLhs[0] += rRhs[0];
        rLhs[1] += rRhs[1];
        rLhs[2] += rRhs[2];
    }

    /// @brief Returns rValue*(rLHs + rRhs)
    /// @param rValue
    /// @param rLhs
    /// @param rRhs
    /// @return Vector3d
    static inline Vector3d AddAndMult( double rValue, const Vector3d& rLhs, const Vector3d& rRhs) {
        return {rValue*(rLhs[0]+rRhs[0]), rValue*(rLhs[1]+rRhs[1]), rValue*(rLhs[2]+rRhs[2])};
    }

    /// @brief Returns: rLHs - rRhs
    /// @param rLhs
    /// @param rRhs
    /// @return Vector3d
    static inline Vector3d Subtract( const Vector3d& rLhs, const Vector3d& rRhs) {
        return {rLhs[0]-rRhs[0], rLhs[1]-rRhs[1], rLhs[2]-rRhs[2]};
    }

    /// @brief Performs: rLHs += rRhs
    /// @param rLhs
    /// @param rRhs
    static inline void SubstractSelf( Vector3d& rLhs, const Vector3d& rRhs) {
        rLhs[0] -= rRhs[0];
        rLhs[1] -= rRhs[1];
        rLhs[2] -= rRhs[2];
    }

    /// @brief Returns: rValue*(rLhs - rRhs)
    /// @param rValue
    /// @param rLhs
    /// @param rRhs
    /// @return
    static inline Vector3d SubstractAndMult( double rValue, const Vector3d& rLhs, const Vector3d& rRhs) {
        return {rValue*(rLhs[0]-rRhs[0]), rValue*(rLhs[1]-rRhs[1]), rValue*(rLhs[2]-rRhs[2])};
    }

    /// @brief Returns rValue*rRhs
    /// @param rValue
    /// @param rRhs
    /// @return Vector3d
    static inline Vector3d Mult( double rValue, const Vector3d& rRhs) {
        return {rValue*rRhs[0], rValue*rRhs[1], rValue*rRhs[2]};
    }

    /// @brief Performs: rLhs += rValue
    /// @param rLhs
    /// @param rValue
    static inline void MultSelf( Vector3d& rLhs, double rValue) {
        rLhs[0] *= rValue;
        rLhs[1] *= rValue;
        rLhs[2] *= rValue;
    }

    /// @brief Returns rLhs*rRhs (Elementwise)
    /// @param rLhs
    /// @param rRhs
    /// @return Vector3d
    static inline Vector3d MultElementWise( const Vector3d& rLhs, const Vector3d& rRhs) {
        return {rLhs[0]*rRhs[0], rLhs[1]*rRhs[1], rLhs[2]*rRhs[2]};
    }

    /// @brief Returns rLhs / rValue
    /// @param rLhs
    /// @param rValue
    /// @return Vector3d
    static inline Vector3d Divide( const Vector3d& rLhs, double rValue ) {
        return {rLhs[0]/rValue, rLhs[1]/rValue, rLhs[2]/rValue};
    }

    /// @brief Performs: rLhs /= rValue
    /// @param rLhs
    /// @param rValue
    static inline void DivideSelf( Vector3d& rLhs, double rValue) {
        rLhs[0] /= rValue;
        rLhs[1] /= rValue;
        rLhs[2] /= rValue;
    }

    ///@brief Cross product of two vectors.
    ///@param Vector3d rLHs
    ///@param Vector3d rRHs
    ///@return double.
    static inline Vector3d Cross( const Vector3d& rLhs, const Vector3d& rRhs) {
        return { rLhs[1]*rRhs[2] - rLhs[2]*rRhs[1],
                 rLhs[2]*rRhs[0] - rLhs[0]*rRhs[2],
                 rLhs[0]*rRhs[1] - rLhs[1]*rRhs[0] };
    }

    ///@brief Norm of vector.
    ///@param Vector3d rLHs
    ///@return double.
    static inline double Norm( const Vector3d& rLhs ) {
        return std::sqrt( rLhs[0]*rLhs[0] + rLhs[1]*rLhs[1] + rLhs[2]*rLhs[2] );
    }

    /// @brief Returns max value of vector
    /// @tparam T
    /// @param rVector
    /// @return T
    template<typename T>
    static inline T Max( const std::array<T, 3>& rVector ){
        return std::max<T>(std::max<T>(rVector[0], rVector[1]), rVector[2]);
    }

    /// @brief Returns min value of vector.
    /// @tparam T
    /// @param rVector
    /// @return T
    template<typename T>
    static T inline Min( const std::array<T, 3>& rVector ){
        return std::min<T>(std::min<T>(rVector[0], rVector[1]), rVector[2]);
    }

} // End namespace math

} // End namespace queso

#endif