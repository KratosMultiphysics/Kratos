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
        const int degree_u,
        const int degree_v,
        const int nb_poles_u,
        const int nb_poles_v,
        const bool is_rational
    ) : m_degree_u(degree_u),
        m_degree_v(degree_v),
        m_nb_poles_u(nb_poles_u),
        m_nb_poles_v(nb_poles_v),
        m_knots_u(nb_poles_u + degree_u - 1),
        m_knots_v(nb_poles_v + degree_v - 1),
        m_poles(nb_poles_u * nb_poles_v),
        m_weights(is_rational ? nb_poles_u * nb_poles_v : 0)
    {
        static_assert(TDimension > 0, "Invalid dimension");
    }

public:     // static methods
    static constexpr int dimension()
    {
        return TDimension;
    }

public:     // methods
    int to_single_index(const int index_u, const int index_v) const
    {
        return index_u * nb_poles_v() + index_v;
    }

    std::pair<int, int> to_double_index(const int index) const
    {
        const int index_u = index / nb_poles_v();
        const int index_v = index % nb_poles_v();
        return {index_u, index_v};
    }

    Vector pole(const int index_u, const int index_v) const
    {
        const int index = to_single_index(index_u, index_v);
        return pole(index);
    }

    Vector& pole(const int index_u, const int index_v)
    {
        const int index = to_single_index(index_u, index_v);
        return pole(index);
    }

    void set_pole(const int index_u, const int index_v, const Vector& value)
    {
        const int index = to_single_index(index_u, index_v);
        m_poles[index] = value;
    }

    double weight(const int index_u, const int index_v) const
    {
        if (is_rational()) {
            const int index = to_single_index(index_u, index_v);
            return m_weights[index];
        } else {
            return 1;
        }
    }

    double& weight(const int index_u, const int index_v)
    {
        if (is_rational()) {
            const int index = to_single_index(index_u, index_v);
            return m_weights[index];
        } else {
            throw std::runtime_error("");
        }
    }

    void set_weight(const int index_u, const int index_v, const double value)
    {
        if (is_rational()) {
            const int index = to_single_index(index_u, index_v);
            m_weights[index] = value;
        } else {
            throw std::invalid_argument("Geometry is not rational");
        }
    }

    bool is_rational() const
    {
        return m_weights.size() != 0;
    }

    int degree_u() const
    {
        return m_degree_u;
    }

    int degree_v() const
    {
        return m_degree_v;
    }

    Interval domain_u() const
    {
        double u0 = knot_u(degree_u() - 1);
        double u1 = knot_u(nb_knots_u() - degree_u());

        return Interval(u0, u1);
    }

    Interval domain_v() const
    {
        double v0 = knot_v(degree_v() - 1);
        double v1 = knot_v(nb_knots_v() - degree_v());

        return Interval(v0, v1);
    }

    int nb_knots_u() const
    {
        return static_cast<int>(m_knots_u.size());
    }

    double knot_u(const int index) const
    {
        return m_knots_u[index];
    }

    void set_knot_u(const int index, const double value)
    {
        m_knots_u[index] = value;
    }

    const std::vector<double>& knots_u() const
    {
        return m_knots_u;
    }

    int nb_knots_v() const
    {
        return static_cast<int>(m_knots_v.size());
    }

    double knot_v(const int index) const
    {
        return m_knots_v[index];
    }

    void set_knot_v(const int index, const double value)
    {
        m_knots_v[index] = value;
    }

    const std::vector<double>& knots_v() const
    {
        return m_knots_v;
    }

    int nb_poles_u() const
    {
        return nb_knots_u() - degree_u() + 1;
    }

    int nb_poles_v() const
    {
        return nb_knots_v() - degree_v() + 1;
    }

    int nb_poles() const
    {
        return nb_poles_u() * nb_poles_v();
    }

    Vector pole(const int index) const
    {
        return m_poles[index];
    }

    Vector& pole(const int index)
    {
        return m_poles[index];
    }

    void set_pole(const int index, const Vector& value)
    {
        m_poles[index] = value;
    }

    const std::vector<Vector>& poles() const
    {
        return m_poles;
    }

    double weight(const int index) const
    {
        return m_weights[index];
    }

    double& weight(const int index)
    {
        return m_weights[index];
    }

    void set_weight(const int index, const double value)
    {
        const int index_u = index / nb_poles_v();
        const int index_v = index % nb_poles_v();

        set_weight(index_u, index_v, value);
    }

    const std::vector<double>& weights() const
    {
        return m_weights;
    }

    Vector point_at(const double u, const double v) const
    {
        // compute shape functions

        NurbsSurfaceShapeFunction shape(degree_u(), degree_v(), 0);

        if (is_rational()) {
            shape.compute(knots_u(), knots_v(), [&](int i, int j) {
            return weight(i, j); }, u, v);
        } else {
            shape.compute(knots_u(), knots_v(), u, v);
        }

        // compute value

        Vector result;

        for (int i = 0; i <= degree_u(); i++) {
            for (int j = 0; j <= degree_v(); j++) {
                int pole_u = shape.first_nonzero_pole_u() + i;
                int pole_v = shape.first_nonzero_pole_v() + j;

                Vector value = pole(pole_u, pole_v) * shape(0, i, j);

                if (i == 0 && j == 0) {
                    result = value;
                } else {
                    result += value;
                }
            }
        }

        return result;
    }

    std::vector<Vector> derivatives_at(const double u, const double v,
        const int order) const
    {
        // compute shape functions

        NurbsSurfaceShapeFunction shape_functions(degree_u(), degree_v(),
            order);

        if (is_rational()) {
            shape_functions.compute(knots_u(), knots_v(), [&](int i, int j) {
                return weight(i, j); }, u, v);
        } else {
            shape_functions.compute(knots_u(), knots_v(), u, v);
        }

        // compute derivatives

        const int nb_shape_functions = shape_functions.nb_shape_functions();

        std::vector<Vector> result(nb_shape_functions);

        for (int k = 0; k < nb_shape_functions; k++) {
            for (int i = 0; i <= degree_u(); i++) {
                for (int j = 0; j <= degree_v(); j++) {
                    const int pole_u =
                        shape_functions.first_nonzero_pole_u() + i;
                    const int pole_v =
                        shape_functions.first_nonzero_pole_v() + j;

                    const Vector value =
                        pole(pole_u, pole_v) * shape_functions(k, i, j);

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

    std::pair<std::vector<int>, Matrix> shape_functions_at(const double u,
        const double v, const int order) const
    {
        NurbsSurfaceShapeFunction shape(degree_u(), degree_v(), order);

        if (is_rational()) {
            shape.compute(knots_u(), knots_v(), [&](int i, int j) -> double {
                return weight(i, j);
            }, u, v);
        } else {
            shape.compute(knots_u(), knots_v(), u, v);
        }

        Matrix shapeFunctions(shape.nb_shape_functions(),
            shape.nb_nonzero_poles());

        for (int i = 0; i < shape.nb_shape_functions(); i++) {
            for (int j = 0; j < shape.nb_nonzero_poles(); j++) {
                shapeFunctions(i, j) = shape(i, j);
            }
        }

        std::vector<int> indices(shape.nb_nonzero_poles());
        auto it = indices.begin();

        for (int i = 0; i < shape.nb_nonzero_poles_u(); i++) {
            for (int j = 0; j < shape.nb_nonzero_poles_v(); j++) {
                const int poleIndex = Math::single_index(nb_poles_u(),
                    nb_poles_v(), shape.first_nonzero_pole_u() + i,
                    shape.first_nonzero_pole_v() + j);

                *(it++) = poleIndex;
            }
        }

        return {indices, shapeFunctions};
    }

    std::vector<Interval> spans_u() const
    {
        const int first_span = degree_u() - 1;
        const int last_span = nb_knots_u() - degree_u() - 1;

        const int nb_spans = last_span - first_span + 1;

        std::vector<Interval> result(nb_spans);

        for (int i = 0; i < nb_spans; i++) {
            const double t0 = knot_u(first_span + i);
            const double t1 = knot_u(first_span + i + 1);

            result[i] = Interval(t0, t1);
        }

        return result;
    }

    std::vector<Interval> spans_v() const
    {
        const int first_span = degree_v() - 1;
        const int last_span = nb_knots_v() - degree_v() - 1;

        const int nb_spans = last_span - first_span + 1;

        std::vector<Interval> result(nb_spans);

        for (int i = 0; i < nb_spans; i++) {
            const double t0 = knot_v(first_span + i);
            const double t1 = knot_v(first_span + i + 1);

            result[i] = Interval(t0, t1);
        }

        return result;
    }

    void reparametrize(Interval domain_u, Interval domain_v)
    {
        const auto old_domain_u = this->domain_u();

        for (auto& knot : m_knots_u) {
            knot = domain_u.parameter_at_normalized(
                old_domain_u.normalized_at(knot));
        }

        const auto old_domain_v = this->domain_v();

        for (auto& knot : m_knots_v) {
            knot = domain_u.parameter_at_normalized(
                old_domain_v.normalized_at(knot));
        }
    }
}; // class NurbsSurfaceGeometry

} // namespace Kratos

#endif // KRATOS_NURBS_SURFACE_GEOMETRY_H_INCLUDED defined