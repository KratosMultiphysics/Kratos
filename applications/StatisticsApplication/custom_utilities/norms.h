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
#include "utilities/string_utilities.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
namespace Norms
{

using IndexType = std::size_t;

class Value {
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = TDataType;

    ///@}
    ///@name Public operations
    ///@{

    template <class TDataType>
    inline TDataType Evaluate(const TDataType& rValue) const
    {
        return rValue;
    }

    ///@}
};

class L2 {
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

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

    ///@}
};


class Infinity {
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

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

    ///@}
};

class VectorComponent {
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

    ///@}
    ///@name Life cycle
    ///@{

    VectorComponent(const IndexType ComponentIndex) : mComponentIndex(ComponentIndex) {}

    ///@}
    ///@name Public operations
    ///@{

    template <class TDataType>
    inline double Evaluate(const TDataType& rValue) const
    {
        if constexpr(std::is_same_v<TDataType, array_1d<double, 3>> ||
                     std::is_same_v<TDataType, array_1d<double, 4>> ||
                     std::is_same_v<TDataType, array_1d<double, 6>> ||
                     std::is_same_v<TDataType, array_1d<double, 9>> ||
                     std::is_same_v<TDataType, Vector>) {

            KRATOS_ERROR_IF(mComponentIndex >= rValue.size())
                << "Component index is larger than the size of the vector [ Component index = "
                << mComponentIndex << ", vector size  = " << rValue.size() << " ].\n";

            return rValue[mComponentIndex];

        } else {
            static_assert(!std::is_same_v<TDataType, TDataType>,
                          "Unsupported data type.");
            return 0;
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const IndexType mComponentIndex;

    ///@}
};

class MatrixComponent {
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

    ///@}
    ///@name Life cycle
    ///@{

    MatrixComponent(
        const IndexType RowIndex,
        const IndexType ColIndex)
        : mRowIndex(RowIndex),
          mColIndex(ColIndex)
    {
    }

    ///@}
    ///@name Public operations
    ///@{

    inline double Evaluate(const Matrix& rValue) const
    {
        KRATOS_ERROR_IF(mRowIndex >= rValue.size1())
            << "Row index is larger than the size1 of the matrix [ Row index = "
            << mRowIndex << ", matrix size1  = " << rValue.size1() << " ].\n";

        KRATOS_ERROR_IF(mColIndex >= rValue.size2())
            << "Col index is larger than the size1 of the matrix [ Col index = "
            << mColIndex << ", matrix size2  = " << rValue.size2() << " ].\n";

        return rValue(mRowIndex, mColIndex);
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const IndexType mRowIndex;

    const IndexType mColIndex;

    ///@}
};

class P
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

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

    ///@}
};

class LPQ
{
public:
    ///@name Type definitions
    ///@{

    template<class TDataType>
    using ResultantValueType = double;

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

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mP;

    const double mQ;

    ///@}
};

template<class TDataType,
        std::enable_if_t<std::disjunction_v<
                                std::is_same<TDataType, int>,
                                std::is_same<TDataType, double>>, bool> = true>
static std::variant<Norms::L2, Norms::Infinity> GetNorm(const std::string& rNormType)
{
    if (rNormType == "l2") {
        return L2();
    } else if (rNormType == "magnitude") {
        return L2();
    } else if (rNormType == "infinity") {
        return Infinity();
    } else {
        KRATOS_ERROR << "Unsupported norm requested for arithmetic data type [ norm_type = \""
                     << rNormType << "\"]. Followings are the supported norm types:"
                     << "\n\t" << "value"
                     << "\n\t" << "l2"
                     << "\n\t" << "magnitude"
                     << "\n\t" << "infinity";
    }
}

template<class TDataType,
        std::enable_if_t<std::disjunction_v<
                                std::is_same<TDataType, array_1d<double, 3>>,
                                std::is_same<TDataType, array_1d<double, 4>>,
                                std::is_same<TDataType, array_1d<double, 6>>,
                                std::is_same<TDataType, array_1d<double, 9>>,
                                std::is_same<TDataType, Vector>>, bool> = true>
static std::variant<Norms::L2, Norms::Infinity, Norms::P, Norms::VectorComponent> GetNorm(const std::string& rNormType)
{
    if (rNormType == "l2") {
        return L2();
    } else if (rNormType == "magnitude") {
        return L2();
    } else if (rNormType == "euclidean") {
        return L2();
    } else if (rNormType == "infinity") {
        return Infinity();
    } else if (rNormType == "component_x") {
        return VectorComponent(0);
    } else if (rNormType == "component_y") {
        return VectorComponent(1);
    } else if (rNormType == "component_z") {
        return VectorComponent(2);
    } else if (rNormType.size() > 6 && rNormType.substr(0, 6) == "pnorm_") {
        const std::string p_str = rNormType.substr(6, rNormType.size() - 6);
        return P(std::stod(p_str));
    } else if (rNormType.size() > 6 && rNormType.substr(0, 6) == "index_") {
        const std::string index_str = rNormType.substr(6, rNormType.size() - 6);
        return VectorComponent(std::stoul(index_str));
    } else {
        KRATOS_ERROR << "Unsupported norm requested for vector data type [ norm_type = \""
                     << rNormType << "\"]. Followings are the supported norm types:"
                     << "\n\t" << "value"
                     << "\n\t" << "l2"
                     << "\n\t" << "magnitude"
                     << "\n\t" << "euclidean"
                     << "\n\t" << "infinity"
                     << "\n\t" << "pnorm_p"
                     << "\n\t" << "component_x"
                     << "\n\t" << "component_y"
                     << "\n\t" << "component_z"
                     << "\n\t" << "index_i";
        return Infinity();
    }
}

template<class TDataType,
        std::enable_if_t<std::is_same_v<TDataType, Matrix>, bool> = true>
static std::variant<Norms::L2, Norms::Infinity, Norms::P, Norms::MatrixComponent, Norms::Trace, Norms::LPQ> GetNorm(const std::string& rNormType)
{
    if (rNormType == "magnitude") {
        return L2();
    } else if (rNormType == "frobenius") {
        return L2();
    } else if (rNormType == "infinity") {
        return Infinity();
    } else if (rNormType.size() > 6 && rNormType.substr(0, 6) == "pnorm_") {
        const std::string p_str = rNormType.substr(6, rNormType.size() - 6);
        return P(std::stod(p_str));
    } else if (rNormType.size() > 7 && rNormType.substr(0, 7) == "index_(") {
        const std::string& r_values_str = rNormType.substr(
            7, rNormType.size() - std::min(8, static_cast<int>(rNormType.size() - 1)));
        const auto& str_indices = StringUtilities::SplitStringByDelimiter(r_values_str, ',');
        return MatrixComponent(std::stoul(str_indices[0]), std::stoul(str_indices[1]));
    } else if (rNormType == "trace") {
        return Trace();
    } else if (rNormType.size() > 9 && rNormType.substr(0, 9) == "lpqnorm_(") {
        const std::string& r_values_str = rNormType.substr(
            9, rNormType.size() - std::min(10, static_cast<int>(rNormType.size() - 1)));
        const auto& str_indices = StringUtilities::SplitStringByDelimiter(r_values_str, ',');
        return LPQ(std::stod(str_indices[0]), std::stod(str_indices[1]));
    } else {
        KRATOS_ERROR << "Unsupported norm requested for matrix data type [ norm_type = \""
                     << rNormType << "\"]. Followings are the supported norm types:"
                     << "\n\t" << "value"
                     << "\n\t" << "frobenius"
                     << "\n\t" << "magnitude"
                     << "\n\t" << "infinity"
                     << "\n\t" << "pnorm_p"
                     << "\n\t" << "index_(i,j)"
                     << "\n\t" << "trace"
                     << "\n\t" << "lpqnorm_(p,q)";
        return Infinity();
    }
}

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
    using type = std::variant<Norms::L2, Norms::Infinity, Norms::P, Norms::VectorComponent>;
};

template<>
struct NormType<Vector>
{
    using type = std::variant<Norms::L2, Norms::Infinity, Norms::P, Norms::VectorComponent>;
};

template<>
struct NormType<Matrix>
{
    using type = std::variant<Norms::L2, Norms::Infinity, Norms::P, Norms::MatrixComponent, Norms::Trace, Norms::LPQ>;
};

} // namespace Norms
} // namespace Kratos