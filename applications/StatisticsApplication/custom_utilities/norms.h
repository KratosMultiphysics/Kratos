//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#pragma once

// System includes
#include <cmath>
#include <string>
#include <type_traits>
#include <variant>

// Application incldues
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
namespace Norms
{

using IndexType = std::size_t;

class L2 {
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

    KRATOS_CLASS_POINTER_DEFINITION(L2);

    ///@}
    ///@name Public operations
    ///@{

    template <class TDataType>
    inline double Evaluate(const TDataType& rValue) const
    {
        if constexpr(std::is_arithmetic_v<TDataType>) {
            return std::abs(rValue);
        } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>> ||
                            std::is_same_v<TDataType, array_1d<double, 4>> ||
                            std::is_same_v<TDataType, array_1d<double, 6>> ||
                            std::is_same_v<TDataType, array_1d<double, 9>> ||
                            std::is_same_v<TDataType, Vector>) {
            return norm_2(rValue);
        } else if constexpr(std::is_same_v<TDataType, Matrix>) {
            return norm_frobenius(rValue);
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>,
                          "Unsupported data type.");
            return 0;
        }
    }

    std::string Info() const { return "L2"; }

    ///@}
    ///@name Public static operations
    ///@{

    static std::string TypeInfo() { return "L2"; }

    ///@}
};

class Infinity {
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

    KRATOS_CLASS_POINTER_DEFINITION(Infinity);

    ///@}
    ///@name Public operations
    ///@{

    template <class TDataType>
    inline double Evaluate(const TDataType& rValue) const
    {
        if constexpr(std::is_arithmetic_v<TDataType>) {
            return std::abs(rValue);
        } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>> ||
                            std::is_same_v<TDataType, array_1d<double, 4>> ||
                            std::is_same_v<TDataType, array_1d<double, 6>> ||
                            std::is_same_v<TDataType, array_1d<double, 9>> ||
                            std::is_same_v<TDataType, Vector> ||
                            std::is_same_v<TDataType, Matrix>) {
            return norm_inf(rValue);
        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>,
                          "Unsupported data type.");
            return 0;
        }
    }

    std::string Info() const { return "Infinity"; }

    ///@}
    ///@name Public static operations
    ///@{

    static std::string TypeInfo() { return "Infinity"; }

    ///@}
};

class P
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

    KRATOS_CLASS_POINTER_DEFINITION(P);

    ///@}
    ///@name Life cycle
    ///@{

    P(const double P): mP(P)
    {
        KRATOS_ERROR_IF(mP < 1.0)
            << "p-norm only supports p >= 1 values. [ " << mP << " !>= 1 ].\n";
    }

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    inline double Evaluate(const TDataType& rValue) const
    {
        if constexpr(std::is_same_v<TDataType, array_1d<double, 3>> ||
                     std::is_same_v<TDataType, array_1d<double, 4>> ||
                     std::is_same_v<TDataType, array_1d<double, 6>> ||
                     std::is_same_v<TDataType, array_1d<double, 9>> ||
                     std::is_same_v<TDataType, Vector> ||
                     std::is_same_v<TDataType, Matrix>) {

            double value = 0.0;

            if constexpr(!std::is_same_v<TDataType, Matrix>) {
                for (IndexType i = 0; i < rValue.size(); ++i) {
                    value += std::pow(std::abs(rValue[i]), mP);
                }
            } else {
                const auto& r_data = rValue.data();
                for (IndexType i = 0; i < r_data.size(); ++i) {
                    value += std::pow(std::abs(r_data[i]), mP);
                }
            }

            return std::pow(value, 1 / mP);

        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>,
                          "Unsupported data type.");
            return 0;
        }
    }

    std::string Info() const { return "P_(" + std::to_string(mP) + ")"; }

    ///@}
    ///@name Public static operations
    ///@{

    static std::string TypeInfo() { return "P"; }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mP;

    ///@}
};

class Trace
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

    KRATOS_CLASS_POINTER_DEFINITION(Trace);

    ///@}
    ///@name Public operations
    ///@{

    inline double Evaluate(const Matrix& rValue) const
    {
        const IndexType n1 = rValue.size1();
        const IndexType n2 = rValue.size2();

        KRATOS_ERROR_IF(n1 != n2)
            << "Trace is only supported for square matrices.\n";

        double result = 0.0;

        for (IndexType i = 0; i < n1; ++i)
        {
            result += rValue(i, i);
        }
        return result;
    }

    std::string Info() const { return "Trace"; }

    ///@}
    ///@name Public static operations
    ///@{

    static std::string TypeInfo() { return "Trace"; }

    ///@}
};

class LPQ
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

    KRATOS_CLASS_POINTER_DEFINITION(LPQ);

    ///@}
    ///@name Life cycle
    ///@{

    LPQ(const double P,
        const double Q)
        : mP(P),
          mQ(Q)
    {
        KRATOS_ERROR_IF(mP < 1.0)
            << "lpqnorm only supports p >= 1 values. [ " << mP << " !>= 1 ].\n";
        KRATOS_ERROR_IF(mQ < 1.0)
            << "lpqnorm only supports q >= 1 values. [ " << mQ << " !>= 1 ].\n";
    }

    ///@}
    ///@name Public operations
    ///@{

    inline double Evaluate(const Matrix& rValue) const
    {
        const IndexType n1 = rValue.size1();
        const IndexType n2 = rValue.size2();
        const double coeff = mQ / mP;

        double result = 0.0;
        for (IndexType j = 0; j < n2; ++j)
        {
            double col_sum = 0.0;
            for (IndexType i = 0; i < n1; ++i)
            {
                col_sum += std::pow(std::abs(rValue(i, j)), mP);
            }
            result += std::pow(col_sum, coeff);
        }

        return std::pow(result, 1.0 / mQ);
    }

    std::string Info() const { return "LPQ_(" + std::to_string(mP) + "," + std::to_string(mQ) + ")"; }

    ///@}
    ///@name Public static operations
    ///@{

    static std::string TypeInfo() { return "LPQ"; }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mP;

    const double mQ;

    ///@}
};

using AllNormTypes = std::variant<
                                Norms::L2,
                                Norms::Infinity,
                                Norms::P,
                                Norms::Trace,
                                Norms::LPQ
                            >;

template<class TDataType>
struct NormType {};

template<>
struct NormType<int>
{
    using type = std::variant<Norms::L2, Norms::Infinity>;
};

template<>
struct NormType<double>
{
    using type = std::variant<Norms::L2, Norms::Infinity>;
};

template<class TDataType, std::size_t TSize>
struct NormType<array_1d<TDataType, TSize>>
{
    using type = std::variant<Norms::L2, Norms::Infinity, Norms::P>;
};

template<>
struct NormType<Vector>
{
    using type = std::variant<Norms::L2, Norms::Infinity, Norms::P>;
};

template<>
struct NormType<Matrix>
{
    using type = std::variant<Norms::L2, Norms::Infinity, Norms::P, Norms::Trace, Norms::LPQ>;
};

} // namespace Norms
} // namespace Kratos