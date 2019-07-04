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

#if !defined(KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED )
#define  KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED

#include "nurbs_curve_shape_functions.h"
#include "nurbs_interval.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"

#include <stdexcept>
#include <vector>

namespace Kratos {

template <int TDimension, typename TPoleContainer =
    std::vector<array_1d<double, TDimension>>>
struct NurbsCurveGeometry
{
public:     // types
    using Type = NurbsCurveGeometry<TDimension>;
    using Vector = array_1d<double, TDimension>;

private:    // variables
    const int m_degree;
    std::vector<double> m_knots;
    TPoleContainer m_poles;
    std::vector<double> m_weights;

public:     // constructors
    NurbsCurveGeometry(
        const int degree,
        int nb_poles,
        bool is_rational
    ) : m_degree(degree),
        m_poles(nb_poles),
        m_weights(is_rational ? nb_poles : 0),
        m_knots(NurbsUtility::nb_knots(degree, nb_poles))
    {
        static_assert(TDimension > 0, "Invalid dimension");
    }

    NurbsCurveGeometry(
        const int degree,
        const std::vector<double>& knots,
        const std::vector<Vector>& poles
    ) : m_degree(degree),
        m_knots(knots),
        m_poles(poles),
        m_weights()
    {
        static_assert(TDimension > 0, "Invalid dimension");

        if (knots.size() != NurbsUtility::nb_knots(degree, poles.size())) {
            throw std::runtime_error("Number of knots and poles do not match");
        }
    }

    NurbsCurveGeometry(
        const int degree,
        const std::vector<double>& knots,
        const std::vector<Vector>& poles,
        const std::vector<double>& weights
    ) : m_degree(degree),
        m_knots(knots),
        m_poles(poles),
        m_weights(weights)
    {
        static_assert(TDimension > 0);

        if (knots.size() != NurbsUtility::nb_knots(degree, poles.size())) {
            throw std::runtime_error("Number of knots and poles do not match");
        }

        if (weights.size() != poles.size()) {
            throw std::runtime_error(
                "Number of poles and weights do not match");
        }
    }

public:     // static methods
    static constexpr int dimension()
    {
        return TDimension;
    }

public:     // methods
    int degree() const
    {
        return m_degree;
    }

    std::vector<Vector> derivatives_at(const double t, const int order) const
    {
        NurbsCurveShapeFunction shape_function;

        shape_function.resize(m_degree, order);

        if (m_weights.size() > 0) {
            shape_function.compute(m_knots, [&](int i) { return weight(i); }, t);
        } else {
            shape_function.compute(m_knots, t);
        }

        std::vector<Vector> derivatives(shape_function.nb_shapes());

        for (int order = 0; order < shape_function.nb_shapes(); order++) {
            for (int i = 0; i < shape_function.nb_nonzero_poles(); i++) {
                int index = shape_function.first_nonzero_pole() + i;

                if (i == 0) {
                    derivatives[order] = pole(index) * shape_function(order, i);
                } else {
                    derivatives[order] += pole(index) * shape_function(order, i);
                }
            }
        }

        return derivatives;
    }

    Interval domain() const
    {
        return Interval(m_knots[m_degree - 1], m_knots[nb_knots() - m_degree]);
    }

    bool is_rational() const
    {
        return m_weights.size() != 0;
    }

    double knot(int i) const
    {
        return m_knots[i];
    }

    const std::vector<double>& knots() const
    {
        return m_knots;
    }

    int nb_knots() const
    {
        return static_cast<int>(m_knots.size());
    }

    int nb_poles() const
    {
        return static_cast<int>(m_poles.size());
    }

    Vector point_at(const double t) const
    {
        NurbsCurveShapeFunction shape_function;

        shape_function.resize(m_degree, 0);

        if (m_weights.size() > 0) {
            shape_function.compute(m_knots,
                [&](int i) { return weight(i); }, t);
        } else {
            shape_function.compute(m_knots, t);
        }

        Vector point;

        for (int i = 0; i < shape_function.nb_nonzero_poles(); i++) {
            const int index = shape_function.first_nonzero_pole() + i;

            if (i == 0) {
                point = pole(index) * shape_function(0, i);
            } else {
                point += pole(index) * shape_function(0, i);
            }
        }

        return point;
    }

    Vector pole(int i) const
    {
        return m_poles.at(i);
    }

    const TPoleContainer& poles() const
    {
        return m_poles;
    }
    
    void set_knot(int i, double value)
    {
        m_knots[i] = value;
    }
    
    void set_pole(int i, Vector value)
    {
        m_poles.at(i) = value;
    }
    
    void set_weight(int i, double value)
    {
        m_weights.at(i) = value;
    }

    std::pair<std::vector<int>, Matrix> shape_functions_at(const double t,
        const int order) const
    {
        NurbsCurveShapeFunction shape_function(degree(), order);

        if (is_rational()) {
            shape_function.compute(knots(), [&](int i) {
                return weight(i); }, t);
        } else {
            shape_function.compute(knots(), t);
        }

        Matrix values(shape_function.nb_shapes(),
            shape_function.nb_nonzero_poles());

        for (int i = 0; i < shape_function.nb_shapes(); i++) {
            for (int j = 0; j < shape_function.nb_nonzero_poles(); j++) {
                values(i, j) = shape_function(i, j);
            }
        }

        return {shape_function.nonzero_pole_indices(), values};
    }

    std::vector<Interval> spans() const
    {
        const int first_span = degree() - 1;
        const int last_span = nb_knots() - degree() - 1;

        const int nb_spans = last_span - first_span + 1;

        std::vector<Interval> result(nb_spans);

        for (int i = 0; i < nb_spans; i++) {
            const double t0 = knot(first_span + i);
            const double t1 = knot(first_span + i + 1);

            result[i] = Interval(t0, t1);
        }

        return result;
    }
    
    double weight(int i) const
    {
        return m_weights.at(i);
    }
    
    const std::vector<double>& weights() const
    {
        return m_weights;
    }
}; // class NurbsCurveGeometry

} // namespace Kratos

#endif // KRATOS_NURBS_CURVE_GEOMETRY_H_INCLUDED defined