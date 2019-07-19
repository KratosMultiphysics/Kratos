//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

#if !defined(KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED )
#define  KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED

#include "nurbs_curve_shape_functions.h"
#include "nurbs_utility.h"

#include <vector>

namespace Kratos {

class NurbsSurfaceShapeFunction
{
private:    // variables
    int m_order;
    NurbsCurveShapeFunction m_shape_u;
    NurbsCurveShapeFunction m_shape_v;
    std::vector<double> m_weighted_sums;
    std::vector<double> m_values;
    int m_first_nonzero_pole_u;
    int m_first_nonzero_pole_v;

public:     // static methods
    static constexpr inline int GetNbShapeFunctions(const int Order) noexcept
    {
        return (1 + Order) * (2 + Order) / 2;
    }

    static constexpr inline int GetShapeIndex(const int DerivativeU,
        const int DerivativeV) noexcept
    {
        return DerivativeV + (DerivativeU + DerivativeV) * (1 + DerivativeU +
            DerivativeV) / 2;
    }

private:    // methods
    double& GetWeightedSum(const int Index)
    {
        return m_weighted_sums[Index];
    }

    double& GetWeightedSum(const int DerivativeU, const int DerivativeV)
    {
        const int index = GetShapeIndex(DerivativeU, DerivativeV);

        return GetWeightedSum(index);
    }

    inline int GetIndex(const int Derivative, const int PoleU, const int PoleV)
        const
    {
        const int pole = NurbsUtility::GetSingleIndex(GetNbNonzeroPolesU(),
            GetNbNonzeroPolesV(), PoleU, PoleV);
        const int index = NurbsUtility::GetSingleIndex(GetNbShapeFunctions(),
            GetNbNonzeroPoles(), Derivative, pole);

        return index;
    }

    double& GetValue(const int Derivative, const int pole)
    {
        const int index = NurbsUtility::GetSingleIndex(GetNbShapeFunctions(),
            GetNbNonzeroPoles(), Derivative, pole);

        return m_values[index];
    }

    double& GetValue(const int Derivative, const int PoleU, const int PoleV)
    {
        const int index = this->GetIndex(Derivative, PoleU, PoleV);

        return m_values[index];
    }

public:     // constructors
    NurbsSurfaceShapeFunction()
    {
    }

    NurbsSurfaceShapeFunction(const int DegreeU, const int DegreeV,
        const int Order)
    {
        Resize(DegreeU, DegreeV, Order);
    }

public:     // methods
    void Resize(const int DegreeU, const int DegreeV, const int Order)
    {
        const int nb_shape_functions = this->GetNbShapeFunctions(Order);
        const int nb_nonzero_poles = (DegreeU + 1) * (DegreeV + 1);

        m_shape_u.Resize(DegreeU, Order);
        m_shape_v.Resize(DegreeV, Order);
        m_values.resize(nb_shape_functions * nb_nonzero_poles);
        m_weighted_sums.resize(nb_shape_functions);

        m_order = Order;
    }

    int GetDegreeU() const
    {
        return m_shape_u.GetDegree();
    }

    int GetDegreeV() const
    {
        return m_shape_v.GetDegree();
    }

    int GetOrder() const
    {
        return m_order;
    }

    int GetNbShapeFunctions() const
    {
        return GetNbShapeFunctions(GetOrder());
    }

    int GetNbNonzeroPolesU() const
    {
        return m_shape_u.GetNbNonzeroPoles();
    }

    int GetNbNonzeroPolesV() const
    {
        return m_shape_v.GetNbNonzeroPoles();
    }

    int GetNbNonzeroPoles() const
    {
        return GetNbNonzeroPolesU() * GetNbNonzeroPolesV();
    }

    std::vector<std::pair<int, int>> GetNonzeroPoleIndices() const
    {
        std::vector<std::pair<int, int>> indices(GetNbNonzeroPoles());

        for (int i = 0; i < GetNbNonzeroPolesU(); i++) {
            for (int j = 0; j < GetNbNonzeroPolesV(); j++) {
                int poleIndex = NurbsUtility::GetSingleIndex(GetNbNonzeroPolesU(),
                    GetNbNonzeroPolesV(), i, j);

                int poleU = GetFirstNonzeroPoleU() + i;
                int poleV = GetFirstNonzeroPoleV() + j;

                indices[poleIndex] = {poleU, poleV};
            }
        }

        return indices;
    }

    const double GetValue(const int Derivative, const int PoleU, const int PoleV)
        const
    {
        const int index = this->GetIndex(Derivative, PoleU, PoleV);

        return m_values[index];
    }

    const double GetValue(const int Derivative, const int pole) const
    {
        const int index = NurbsUtility::GetSingleIndex(GetNbShapeFunctions(),
            GetNbNonzeroPoles(), Derivative, pole);

        return m_values[index];
    }

    double operator()(const int Derivative, const int pole) const
    {
        return GetValue(Derivative, pole);
    }

    double operator()(const int Derivative, const int PoleU,
        const int PoleV) const
    {
        return GetValue(Derivative, PoleU, PoleV);
    }

    int GetFirstNonzeroPoleU() const
    {
        return m_first_nonzero_pole_u;
    }

    int GetLastNonzeroPoleU() const
    {
        return GetFirstNonzeroPoleU() + GetDegreeU();
    }

    int GetFirstNonzeroPoleV() const
    {
        return m_first_nonzero_pole_v;
    }

    int GetLastNonzeroPoleV() const
    {
        return GetFirstNonzeroPoleV() + GetDegreeV();
    }

    void ComputeAtSpan(const std::vector<double>& rKnotsU,
        const std::vector<double>& rKnotsV, const int SpanU, const int SpanV,
        const double ParameterU, const double ParameterV)
    {
        const int nbvalues = GetNbShapeFunctions() * GetNbNonzeroPoles();

        std::fill(m_values.begin(), m_values.begin() + nbvalues, 0);

        m_first_nonzero_pole_u = SpanU - GetDegreeU() + 1;
        m_first_nonzero_pole_v = SpanV - GetDegreeV() + 1;

        // compute 1D shape functions

        m_shape_u.ComputeAtSpan(rKnotsU, SpanU, ParameterU);
        m_shape_v.ComputeAtSpan(rKnotsV, SpanV, ParameterV);

        // compute 2D shape functions

        for (int i = 0; i <= GetOrder(); i++) {
            for (int j = 0; j <= GetOrder() - i; j++) {
                for (int a = 0; a < GetNbNonzeroPolesU(); a++) {
                    for (int b = 0; b < GetNbNonzeroPolesV(); b++) {
                        const int index = GetShapeIndex(i, j);

                        GetValue(index, a, b) = m_shape_u(i, a) * m_shape_v(j, b);
                    }
                }
            }
        }
    }

    void Compute(const std::vector<double>& rKnotsU,
        const std::vector<double>& rKnotsV, const double ParameterU,
        const double ParameterV)
    {
        const int SpanU = NurbsUtility::GetLowerSpan(GetDegreeU(), rKnotsU, ParameterU);
        const int SpanV = NurbsUtility::GetLowerSpan(GetDegreeV(), rKnotsV, ParameterV);

        ComputeAtSpan(rKnotsU, rKnotsV, SpanU, SpanV, ParameterU, ParameterV);
    }

    template <typename TWeights>
    void ComputeAtSpan(const std::vector<double>& rKnotsU,
        const std::vector<double>& rKnotsV, const int SpanU, const int SpanV,
        const TWeights& Weights, const double ParameterU,
        const double ParameterV)
    {
        // compute B-Spline shape

        ComputeAtSpan(rKnotsU, rKnotsV, SpanU, SpanV, ParameterU, ParameterV);

        // apply weights

        for (int shape = 0; shape < GetNbShapeFunctions(); shape++) {
            GetWeightedSum(shape) = double(0);

            for (int i = 0; i < GetNbNonzeroPolesU(); i++) {
                for (int j = 0; j < GetNbNonzeroPolesV(); j++) {
                    const int poleU = GetFirstNonzeroPoleU() + i;
                    const int poleV = GetFirstNonzeroPoleV() + j;

                    const double weight = Weights(poleU, poleV);
                    GetValue(shape, i, j) *= weight;
                    GetWeightedSum(shape) += GetValue(shape, i, j);
                }
            }
        }

        for (int k = 0; k <= GetOrder(); k++) {
            for (int l = 0; l <= GetOrder() - k; l++) {
                const int shape = GetShapeIndex(k, l);

                for (int j = 1; j <= l; j++) {
                    const int index = GetShapeIndex(k, l - j);

                    double a = NurbsUtility::Binom(l, j) * GetWeightedSum(0, j);

                    for (int p = 0; p < GetNbNonzeroPoles(); p++) {
                        GetValue(shape, p) -= a * GetValue(index, p);
                    }
                }

                for (int i = 1; i <= k; i++) {
                    const int index = GetShapeIndex(k - i, l);

                    double a = NurbsUtility::Binom(k, i) * GetWeightedSum(i, 0);

                    for (int p = 0; p < GetNbNonzeroPoles(); p++) {
                        GetValue(shape, p) -= a * GetValue(index, p);
                    }
                }

                for (int i = 1; i <= k; i++) {
                    const double a = NurbsUtility::Binom(k, i);

                    for (int j = 1; j <= l; j++) {
                        const int index = GetShapeIndex(k - i, l - j);

                        const double b = a * NurbsUtility::Binom(l, j) *
                            GetWeightedSum(i, j);

                        for (int p = 0; p < GetNbNonzeroPoles(); p++) {
                            GetValue(shape, p) -= b * GetValue(index, p);
                        }
                    }
                }

                for (int p = 0; p < GetNbNonzeroPoles(); p++) {
                    GetValue(shape, p) /= GetWeightedSum(0);
                }
            }
        }
    }

    template <typename TWeights>
    void Compute(const std::vector<double>& rKnotsU,
        const std::vector<double>& rKnotsV, const TWeights& Weights,
        const double ParameterU, const double ParameterV)
    {
        const int SpanU = NurbsUtility::GetLowerSpan(GetDegreeU(), rKnotsU, ParameterU);
        const int SpanV = NurbsUtility::GetLowerSpan(GetDegreeV(), rKnotsV, ParameterV);

        ComputeAtSpan(rKnotsU, rKnotsV, SpanU, SpanV, Weights, ParameterU, ParameterV);
    }
}; // class NurbsSurfaceShapeFunction

} // namespace Kratos

#endif // KRATOS_NURBS_SURFACE_SHAPE_FUNCTIONS_H_INCLUDED defined