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

#if !defined(KRATOS_NURBS_SURFACE_GEOMETRY_H_INCLUDED )
#define  KRATOS_NURBS_SURFACE_GEOMETRY_H_INCLUDED

#include "nurbs_interval.h"
#include "nurbs_surface_shape_functions.h"
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

#include <stdexcept>
#include <vector>

namespace Kratos {

template <int TDimension, typename TPoleContainer =
    std::vector<array_1d<double, TDimension>>>
class NurbsSurfaceGeometry
{
public:     // types
    using Type = NurbsSurfaceGeometry<TDimension>;
    using Vector = array_1d<double, TDimension>;

private:    // variables
    int m_degree_u;
    int m_degree_v;
    int m_nb_poles_u;
    int m_nb_poles_v;
    std::vector<double> m_knots_u;
    std::vector<double> m_knots_v;
    TPoleContainer m_poles;
    std::vector<double> m_weights;

public:     // constructor
    NurbsSurfaceGeometry(
        const int DegreeU,
        const int DegreeV,
        const int NbPolesU,
        const int NbPolesV,
        const bool IsRational
    ) : m_degree_u(DegreeU),
        m_degree_v(DegreeV),
        m_nb_poles_u(NbPolesU),
        m_nb_poles_v(NbPolesV),
        m_knots_u(NbPolesU + DegreeU - 1),
        m_knots_v(NbPolesV + DegreeV - 1),
        m_poles(NbPolesU * NbPolesV),
        m_weights(IsRational ? NbPolesU * NbPolesV : 0)
    {
        static_assert(TDimension > 0, "Invalid dimension");
    }

public:     // static methods
    static constexpr int GetDimension()
    {
        return TDimension;
    }

public:     // methods
    int ToSingleIndex(const int IndexU, const int IndexV) const
    {
        return IndexU * GetNbPolesV() + IndexV;
    }

    std::pair<int, int> ToDoubleIndex(const int Index) const
    {
        const int index_u = Index / GetNbPolesV();
        const int index_v = Index % GetNbPolesV();
        return {index_u, index_v};
    }

    Vector GetPole(const int IndexU, const int IndexV) const
    {
        const int index = ToSingleIndex(IndexU, IndexV);
        return GetPole(index);
    }

    Vector& GetPole(const int IndexU, const int IndexV)
    {
        const int index = ToSingleIndex(IndexU, IndexV);
        return GetPole(index);
    }

    void SetPole(const int IndexU, const int IndexV, const Vector& Value)
    {
        const int index = ToSingleIndex(IndexU, IndexV);
        m_poles[index] = Value;
    }

    double GetWeight(const int IndexU, const int IndexV) const
    {
        if (IsRational()) {
            const int index = ToSingleIndex(IndexU, IndexV);
            return m_weights[index];
        } else {
            return 1;
        }
    }

    double& GetWeight(const int IndexU, const int IndexV)
    {
        if (IsRational()) {
            const int index = ToSingleIndex(IndexU, IndexV);
            return m_weights[index];
        } else {
            throw std::runtime_error("");
        }
    }

    void SetWeight(const int IndexU, const int IndexV, const double Value)
    {
        if (IsRational()) {
            const int index = ToSingleIndex(IndexU, IndexV);
            m_weights[index] = Value;
        } else {
            throw std::invalid_argument("Geometry is not rational");
        }
    }

    bool IsRational() const
    {
        return m_weights.size() != 0;
    }

    int GetDegreeU() const
    {
        return m_degree_u;
    }

    int GetDegreeV() const
    {
        return m_degree_v;
    }

    Interval GetDomainU() const
    {
        double u0 = GetKnotU(GetDegreeU() - 1);
        double u1 = GetKnotU(GetNbKnotsU() - GetDegreeU());

        return Interval(u0, u1);
    }

    Interval GetDomainV() const
    {
        double v0 = GetKnotV(GetDegreeV() - 1);
        double v1 = GetKnotV(GetNbKnotsV() - GetDegreeV());

        return Interval(v0, v1);
    }

    int GetNbKnotsU() const
    {
        return static_cast<int>(m_knots_u.size());
    }

    double GetKnotU(const int Index) const
    {
        return m_knots_u[Index];
    }

    void SetKnotU(const int Index, const double Value)
    {
        m_knots_u[Index] = Value;
    }

    const std::vector<double>& GetKnotsU() const
    {
        return m_knots_u;
    }

    int GetNbKnotsV() const
    {
        return static_cast<int>(m_knots_v.size());
    }

    double GetKnotV(const int Index) const
    {
        return m_knots_v[Index];
    }

    void SetKnotV(const int Index, const double Value)
    {
        m_knots_v[Index] = Value;
    }

    const std::vector<double>& GetKnotsV() const
    {
        return m_knots_v;
    }

    int GetNbPolesU() const
    {
        return GetNbKnotsU() - GetDegreeU() + 1;
    }

    int GetNbPolesV() const
    {
        return GetNbKnotsV() - GetDegreeV() + 1;
    }

    int GetNbPoles() const
    {
        return GetNbPolesU() * GetNbPolesV();
    }

    Vector GetPole(const int Index) const
    {
        return m_poles[Index];
    }

    Vector& GetPole(const int Index)
    {
        return m_poles[Index];
    }

    void SetPole(const int Index, const Vector& Value)
    {
        m_poles[Index] = Value;
    }

    const std::vector<Vector>& GetPoles() const
    {
        return m_poles;
    }

    double GetWeight(const int Index) const
    {
        return m_weights[Index];
    }

    double& GetWeight(const int Index)
    {
        return m_weights[Index];
    }

    void SetWeight(const int Index, const double Value)
    {
        const int index_u = Index / GetNbPolesV();
        const int index_v = Index % GetNbPolesV();

        SetWeight(index_u, index_v, Value);
    }

    const std::vector<double>& GetWeights() const
    {
        return m_weights;
    }

    Vector GetPointAt(const double ParameterU, const double ParameterV) const
    {
        // compute shape functions

        NurbsSurfaceShapeFunction shape(GetDegreeU(), GetDegreeV(), 0);

        if (IsRational()) {
            shape.Compute(GetKnotsU(), GetKnotsV(), [&](int i, int j) {
            return GetWeight(i, j); }, ParameterU, ParameterV);
        } else {
            shape.Compute(GetKnotsU(), GetKnotsV(), ParameterU, ParameterV);
        }

        // compute value

        Vector result;

        for (int i = 0; i <= GetDegreeU(); i++) {
            for (int j = 0; j <= GetDegreeV(); j++) {
                int pole_u = shape.GetFirstNonzeroPoleU() + i;
                int pole_v = shape.GetFirstNonzeroPoleV() + j;

                Vector value = GetPole(pole_u, pole_v) * shape(0, i, j);

                if (i == 0 && j == 0) {
                    result = value;
                } else {
                    result += value;
                }
            }
        }

        return result;
    }

    std::vector<Vector> GetDerivativesAt(const double ParameterU,
        const double ParameterV, const int Order) const
    {
        // compute shape functions

        NurbsSurfaceShapeFunction shape_functions(GetDegreeU(), GetDegreeV(),
            Order);

        if (IsRational()) {
            shape_functions.Compute(GetKnotsU(), GetKnotsV(),
                [&](int i, int j) { return GetWeight(i, j); },
                ParameterU, ParameterV);
        } else {
            shape_functions.Compute(GetKnotsU(), GetKnotsV(), ParameterU,
                ParameterV);
        }

        // compute derivatives

        const int nb_shape_functions = shape_functions.GetNbShapeFunctions();

        std::vector<Vector> result(nb_shape_functions);

        for (int k = 0; k < nb_shape_functions; k++) {
            for (int i = 0; i <= GetDegreeU(); i++) {
                for (int j = 0; j <= GetDegreeV(); j++) {
                    const int pole_u =
                        shape_functions.GetFirstNonzeroPoleU() + i;
                    const int pole_v =
                        shape_functions.GetFirstNonzeroPoleV() + j;

                    const Vector value =
                        GetPole(pole_u, pole_v) * shape_functions(k, i, j);

                    if (i == 0 && j == 0) {
                        result[k] = value;
                    } else {
                        result[k] += value;
                    }
                }
            }
        }

        return result;
    }

    std::pair<std::vector<int>, Matrix> GetShapeFunctionsAt(
        const double ParameterU, const double ParameterV, const int Order) const
    {
        NurbsSurfaceShapeFunction shape(GetDegreeU(), GetDegreeV(), Order);

        if (IsRational()) {
            shape.Compute(GetKnotsU(), GetKnotsV(), [&](int i, int j) -> double {
                return GetWeight(i, j);
            }, ParameterU, ParameterV);
        } else {
            shape.Compute(GetKnotsU(), GetKnotsV(), ParameterU, ParameterV);
        }

        Matrix shapeFunctions(shape.GetNbShapeFunctions(),
            shape.nb_nonzero_poles());

        for (int i = 0; i < shape.GetNbShapeFunctions(); i++) {
            for (int j = 0; j < shape.nb_nonzero_poles(); j++) {
                shapeFunctions(i, j) = shape(i, j);
            }
        }

        std::vector<int> indices(shape.nb_nonzero_poles());
        auto it = indices.begin();

        for (int i = 0; i < shape.nb_nonzero_poles_u(); i++) {
            for (int j = 0; j < shape.nb_nonzero_poles_v(); j++) {
                const int poleIndex = Math::single_index(GetNbPolesU(),
                    GetNbPolesV(), shape.GetFirstNonzeroPoleU() + i,
                    shape.GetFirstNonzeroPoleV() + j);

                *(it++) = poleIndex;
            }
        }

        return {indices, shapeFunctions};
    }

    std::vector<Interval> GetSpansU() const
    {
        const int first_span = GetDegreeU() - 1;
        const int last_span = GetNbKnotsU() - GetDegreeU() - 1;

        const int nb_spans = last_span - first_span + 1;

        std::vector<Interval> result(nb_spans);

        for (int i = 0; i < nb_spans; i++) {
            const double t0 = GetKnotU(first_span + i);
            const double t1 = GetKnotU(first_span + i + 1);

            result[i] = Interval(t0, t1);
        }

        return result;
    }

    std::vector<Interval> GetSpansV() const
    {
        const int first_span = GetDegreeV() - 1;
        const int last_span = GetNbKnotsV() - GetDegreeV() - 1;

        const int nb_spans = last_span - first_span + 1;

        std::vector<Interval> result(nb_spans);

        for (int i = 0; i < nb_spans; i++) {
            const double t0 = GetKnotV(first_span + i);
            const double t1 = GetKnotV(first_span + i + 1);

            result[i] = Interval(t0, t1);
        }

        return result;
    }
}; // class NurbsSurfaceGeometry

} // namespace Kratos

#endif // KRATOS_NURBS_SURFACE_GEOMETRY_H_INCLUDED defined