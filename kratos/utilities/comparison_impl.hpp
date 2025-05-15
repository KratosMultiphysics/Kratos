//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// System includes
#include <type_traits> // std::is_same_v
#include <algorithm> // std::min
#include <limits> // std::numeric_limits
#include <cmath> // std::abs


namespace Kratos::Impl {


template <class T>
struct IntegerComparison
{
    static_assert(
           std::is_same_v<T,short>
        || std::is_same_v<T,int>
        || std::is_same_v<T,long>
        || std::is_same_v<T,long long>
        || std::is_same_v<T,unsigned short>
        || std::is_same_v<T,unsigned>
        || std::is_same_v<T,unsigned long>
        || std::is_same_v<T,unsigned long long>,
        "attempting to instantiate IntegerComparison for unsupported type"
    );

    struct Equal {
        constexpr bool operator()(T Left, T Right) const noexcept {
            return Left == Right;
        }
    }; // struct Equal

    struct Less {
        constexpr bool operator()(T Left, T Right) const noexcept {
            return Left < Right;
        }
    }; // struct Less
}; // struct IntegerComparison


/// @brief Relaxed floating point equality comparison and ordering.
/// @details This class implements a relaxed equality comparison for floating point
///          numbers that switches between relative and absolute tolerances depending
///          on the magnitude of numbers to be compared. This provides greater
///          robustness on a wider range.
///          The ordering builds on top of the equality operator.
/// @see https://stackoverflow.com/a/32334103/12350793
template <class T>
struct FloatComparison {
    static_assert(
           std::is_same_v<T,float>
        || std::is_same_v<T,double>
        || std::is_same_v<T,long double>,
        "attempting to instantiate FloatComparison for unsupported type"
    );

    class Less;

    class Equal {
    public:
        constexpr Equal() noexcept
            : Equal(static_cast<T>(0), static_cast<T>(0))
        {}

        constexpr Equal(T AbsoluteTolerance, T RelativeTolerance) noexcept
            : mAbsoluteTolerance(AbsoluteTolerance),
              mRelativeTolerance(RelativeTolerance)
        {}

        constexpr bool operator()(T Left, T Right) const noexcept {
            const T norm = std::min(std::abs(Left) + std::abs(Right),
                                    std::numeric_limits<T>::max());
            return std::abs(Left - Right) < std::max(mAbsoluteTolerance, mRelativeTolerance * norm);
        }

    private:
        T mAbsoluteTolerance, mRelativeTolerance;
    }; // struct Equal

    class Less {
    public:
        constexpr Less() noexcept
            : Less(static_cast<T>(0), static_cast<T>(0))
        {}

        constexpr Less(T AbsoluteTolerance, T RelativeTolerance) noexcept
            : mEqualityComparison(AbsoluteTolerance, RelativeTolerance)
        {}

        constexpr bool operator()(T Left, T Right) const noexcept {
            return Left < Right && !mEqualityComparison(Left, Right);
        }

    private:
        FloatComparison<T>::Equal mEqualityComparison;
    }; // class Less
}; // struct FloatComparison


} // namespace Kratos::Impl
