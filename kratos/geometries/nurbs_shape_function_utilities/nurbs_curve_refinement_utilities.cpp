//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//

#include "nurbs_curve_refinement_utilities.h"

namespace Kratos {

void NurbsCurveRefinementUtilities::KnotRefinement(
    NurbsCurveGeometryType& rGeometry,
    std::vector<double>& rKnotsToInsert,
    PointerVector<NodeType>& rPointsRefined,
    Vector& rKnotsRefined,
    Vector& rWeightsRefined)
{
    const SizeType degree = rGeometry.PolynomialDegree(0);

    const Vector knots_old = CreateFullKnotVector(rGeometry.Knots(), degree);

    SortAndFilter(rKnotsToInsert, knots_old);

    if (rKnotsToInsert.empty()) {
        KRATOS_WARNING("::NurbsCurveRefinementUtilities::KnotRefinement")
            << "No refinement applied as knot vector is empty." << std::endl;
        return;
    }

    const SizeType nb_cp_old = rGeometry.size();
    const SizeType nb_knots_old = knots_old.size();
    const SizeType nb_knots_to_insert = rKnotsToInsert.size();

    KRATOS_ERROR_IF_NOT(nb_knots_old == nb_cp_old + degree + 1)
        << "Invalid full knot vector. knots.size() = " << nb_knots_old
        << ", points.size() = " << nb_cp_old
        << ", degree = " << degree
        << ", expected = " << nb_cp_old + degree + 1
        << std::endl;

    Vector weights_old = rGeometry.Weights();
    if (weights_old.size() != nb_cp_old) {
        weights_old.resize(nb_cp_old);
        std::fill(weights_old.begin(), weights_old.end(), 1.0);
    }

    const SizeType n = nb_cp_old - 1;
    const SizeType m = nb_knots_old - 1;
    const SizeType r = nb_knots_to_insert - 1;

    const SizeType a = NurbsUtilities::GetUpperSpan(
        degree,
        knots_old,
        rKnotsToInsert.front()
    );

    SizeType b = NurbsUtilities::GetUpperSpan(
        degree,
        knots_old,
        rKnotsToInsert.back()
    );

    b += 1;

    const SizeType nb_cp_refined = nb_cp_old + nb_knots_to_insert;
    const SizeType nb_knots_refined = nb_knots_old + nb_knots_to_insert;

    rPointsRefined.resize(nb_cp_refined);
    rWeightsRefined.resize(nb_cp_refined, false);
    rKnotsRefined.resize(nb_knots_refined, false);

    rWeightsRefined = ZeroVector(nb_cp_refined);
    rKnotsRefined = ZeroVector(nb_knots_refined);

    for (IndexType i = 0; i <= a - degree; ++i) {
        const array_1d<double, 3> cp_coordinates =
            rGeometry[i] * weights_old[i];

        rPointsRefined(i) =
            Kratos::make_intrusive<NodeType>(0, cp_coordinates);

        rWeightsRefined[i] = weights_old[i];
    }

    for (IndexType i = b - 1; i <= n; ++i) {
        const IndexType refined_index = i + nb_knots_to_insert;

        KRATOS_ERROR_IF(refined_index >= nb_cp_refined)
            << "Out-of-bounds refined_index = " << refined_index
            << ", nb_cp_refined = " << nb_cp_refined << std::endl;

        const array_1d<double, 3> cp_coordinates =
            rGeometry[i] * weights_old[i];

        rPointsRefined(refined_index) =
            Kratos::make_intrusive<NodeType>(0, cp_coordinates);

        rWeightsRefined[refined_index] = weights_old[i];
    }

    for (IndexType i = 0; i <= a; ++i) {
        rKnotsRefined[i] = knots_old[i];
    }

    for (IndexType i = b + degree; i <= m; ++i) {
        const IndexType refined_index = i + nb_knots_to_insert;

        KRATOS_ERROR_IF(refined_index >= nb_knots_refined)
            << "Out-of-bounds refined knot index = " << refined_index
            << ", nb_knots_refined = " << nb_knots_refined << std::endl;

        rKnotsRefined[refined_index] = knots_old[i];
    }

    IndexType i = b + degree - 1;
    IndexType k = b + degree + r;

    for (int j = static_cast<int>(r); j >= 0; --j) {
        while (i > a && rKnotsToInsert[j] <= knots_old[i]) {
            const IndexType old_index = i - degree - 1;
            const IndexType refined_index = k - degree - 1;

            const array_1d<double, 3> cp_coordinates =
                rGeometry[old_index] * weights_old[old_index];

            rPointsRefined(refined_index) =
                Kratos::make_intrusive<NodeType>(0, cp_coordinates);

            rWeightsRefined[refined_index] = weights_old[old_index];
            rKnotsRefined[k] = knots_old[i];

            --k;
            --i;
        }

        rPointsRefined(k - degree - 1) =
            rPointsRefined(k - degree);

        rWeightsRefined[k - degree - 1] =
            rWeightsRefined[k - degree];

        for (IndexType l = 1; l <= degree; ++l) {
            const IndexType index = k - degree + l;
            const IndexType old_knot_index = i - degree + l;

            double alpha =
                rKnotsRefined[k + l] - rKnotsToInsert[j];

            if (std::abs(alpha) < 1.0e-7) {
                rPointsRefined(index - 1) =
                    rPointsRefined(index);

                rWeightsRefined[index - 1] =
                    rWeightsRefined[index];
            } else {
                const double denominator =
                    rKnotsRefined[k + l] - knots_old[old_knot_index];

                KRATOS_ERROR_IF(std::abs(denominator) < 1.0e-14)
                    << "Zero denominator in alpha computation." << std::endl;

                alpha /= denominator;

                const array_1d<double, 3> cp_coordinates =
                    alpha * rPointsRefined[index - 1]
                    + (1.0 - alpha) * rPointsRefined[index];

                rPointsRefined(index - 1) =
                    Kratos::make_intrusive<NodeType>(0, cp_coordinates);

                rWeightsRefined[index - 1] =
                    alpha * rWeightsRefined[index - 1]
                    + (1.0 - alpha) * rWeightsRefined[index];
            }
        }

        rKnotsRefined[k] = rKnotsToInsert[j];

        --k;
    }

    for (IndexType i_cp = 0; i_cp < nb_cp_refined; ++i_cp) {
        rPointsRefined(i_cp)->Coordinates() =
            rPointsRefined(i_cp)->Coordinates()
            / rWeightsRefined[i_cp];
    }

    rKnotsRefined = CreateReducedKnotVector(rKnotsRefined, degree);
}

void NurbsCurveRefinementUtilities::DegreeElevation(
    NurbsCurveGeometryType& rGeometry,
    SizeType rDegreeToElevate,
    PointerVector<NodeType>& rPointsRefined,
    Vector& rKnotsRefined,
    Vector& rWeightsRefined)
{
    const Vector& knots_old = rGeometry.Knots();

    if (rDegreeToElevate == 0) {
        KRATOS_WARNING("::NurbsCurveRefinementUtilities::DegreeElevation")
            << "No elevation applied as degree is zero." << std::endl;
        return;
    }

    const SizeType degree_old = rGeometry.PolynomialDegree(0);
    const SizeType degree = degree_old + rDegreeToElevate;

    const SizeType nb_cp_old = rGeometry.size();
    const SizeType nb_knots_old = knots_old.size();

    Vector weights_old = rGeometry.Weights();
    if (weights_old.size() != nb_cp_old) {
        weights_old.resize(nb_cp_old);
        std::fill(weights_old.begin(), weights_old.end(), 1.0);
    }

    SizeType number_of_non_zero_spans = 0;
    for (IndexType i = 0; i < nb_knots_old - 1; ++i) {
        if (knots_old[i] != knots_old[i + 1]) {
            ++number_of_non_zero_spans;
        }
    }

    const SizeType nb_knots_refined =
        nb_knots_old + rDegreeToElevate * (number_of_non_zero_spans + 1);

    const SizeType nb_cp_refined =
        nb_knots_refined - degree + 1;

    rKnotsRefined.resize(nb_knots_refined, false);
    rWeightsRefined.resize(nb_cp_refined, false);
    rPointsRefined.resize(nb_cp_refined);

    rKnotsRefined = ZeroVector(nb_knots_refined);
    rWeightsRefined = ZeroVector(nb_cp_refined);

    Matrix bezier_alphas = ZeroMatrix(degree + 1, degree_old + 1);

    PointerVector<NodeType> bezier_points;
    bezier_points.resize(degree_old + 1);

    Vector bezier_points_weights = ZeroVector(degree_old + 1);

    PointerVector<NodeType> elevate_bezier_points;
    elevate_bezier_points.resize(degree + 1);

    Vector elevate_bezier_points_weights = ZeroVector(degree + 1);

    PointerVector<NodeType> next_bezier_points;
    if (degree_old > 1) {
        next_bezier_points.resize(degree_old - 1);
    }

    Vector next_bezier_points_weights =
        degree_old > 1 ? ZeroVector(degree_old - 1) : ZeroVector(0);

    bezier_alphas(0, 0) = 1.0;
    bezier_alphas(degree, degree_old) = 1.0;

    for (IndexType i = 1; i < degree / 2 + 1; ++i) {
        const double inv =
            1.0 / NurbsUtilities::GetBinomCoefficient(degree, i);

        const SizeType min_degree = std::min(degree_old, i);
        const int index = std::max(
            0,
            static_cast<int>(i - rDegreeToElevate)
        );

        for (IndexType j = index; j < min_degree + 1; ++j) {
            bezier_alphas(i, j) =
                inv
                * NurbsUtilities::GetBinomCoefficient(degree_old, j)
                * NurbsUtilities::GetBinomCoefficient(rDegreeToElevate, i - j);

            bezier_alphas(degree - i, degree_old - j) =
                bezier_alphas(i, j);
        }
    }

    SizeType a = degree_old - 1;
    SizeType b = degree_old;

    int r = -1;

    IndexType knot_index = degree + 1;
    IndexType control_point_index = 1;

    double knot_old_a = knots_old[a];

    for (IndexType i = 0; i < degree; ++i) {
        rKnotsRefined[i] = knot_old_a;
    }

    {
        const array_1d<double, 3> cp_coordinates =
            rGeometry[0] * weights_old[0];

        rPointsRefined(0) =
            Kratos::make_intrusive<NodeType>(0, cp_coordinates);

        rWeightsRefined[0] = weights_old[0];
    }

    for (IndexType i = 0; i < degree_old + 1; ++i) {
        const array_1d<double, 3> cp_coordinates =
            rGeometry[i] * weights_old[i];

        bezier_points(i) =
            Kratos::make_intrusive<NodeType>(0, cp_coordinates);

        bezier_points_weights[i] = weights_old[i];
    }

    while (b < nb_knots_old) {
        IndexType i = b;

        while (b < nb_knots_old - 1 && knots_old[b] == knots_old[b + 1]) {
            ++b;
        }

        if (b + 1 == nb_knots_old) {
            ++b;
        }

        const IndexType mult = b - i + 1;
        const double knot_old_b = knots_old[i];

        const int r_old = r;
        r = static_cast<int>(degree_old - mult);

        IndexType left_bezier;
        IndexType right_bezier;

        if (r_old > 0) {
            left_bezier = (r_old + 2) / 2;
        } else {
            left_bezier = 1;
        }

        if (r > 0) {
            right_bezier = degree - ((r + 1) / 2);
        } else {
            right_bezier = degree;
        }

        if (r > 0) {
            Vector alphas = ZeroVector(degree_old - 1);

            for (IndexType k = degree_old; k > mult; --k) {
                alphas[k - mult - 1] =
                    (knot_old_b - knot_old_a)
                    / (knots_old[a + k] - knot_old_a);
            }

            for (IndexType j = 1; j < static_cast<IndexType>(r + 1); ++j) {
                const int save = r - static_cast<int>(j);
                const IndexType s = mult + j;

                for (IndexType k = degree_old; k >= s; --k) {
                    const array_1d<double, 3> cp_coordinates =
                        alphas[k - s] * bezier_points[k]
                        + (1.0 - alphas[k - s]) * bezier_points[k - 1];

                    bezier_points(k) =
                        Kratos::make_intrusive<NodeType>(0, cp_coordinates);

                    bezier_points_weights[k] =
                        alphas[k - s] * bezier_points_weights[k]
                        + (1.0 - alphas[k - s]) * bezier_points_weights[k - 1];

                    if (k == s) {
                        break;
                    }
                }

                if (degree_old > 1) {
                    next_bezier_points(save) = bezier_points(degree_old);
                    next_bezier_points_weights[save] =
                        bezier_points_weights[degree_old];
                }
            }
        }

        for (IndexType i_elev = left_bezier; i_elev < degree + 1; ++i_elev) {
            const array_1d<double, 3> zero_coordinates = ZeroVector(3);

            elevate_bezier_points(i_elev) =
                Kratos::make_intrusive<NodeType>(0, zero_coordinates);

            elevate_bezier_points_weights[i_elev] = 0.0;

            const SizeType min_degree = std::min(degree_old, i_elev);
            const int index = std::max(
                0,
                static_cast<int>(i_elev - rDegreeToElevate)
            );

            for (IndexType j = index; j < min_degree + 1; ++j) {
                const array_1d<double, 3> cp_coordinates =
                    elevate_bezier_points[i_elev]
                    + bezier_alphas(i_elev, j) * bezier_points[j];

                elevate_bezier_points(i_elev) =
                    Kratos::make_intrusive<NodeType>(0, cp_coordinates);

                elevate_bezier_points_weights[i_elev] +=
                    bezier_alphas(i_elev, j) * bezier_points_weights[j];
            }
        }

        if (r_old > 1) {
            IndexType first = knot_index - 2;
            IndexType last = knot_index;

            const double beta =
                (knot_old_b - rKnotsRefined[knot_index - 2])
                / (knot_old_b - knot_old_a);

            for (IndexType t = 1; t < static_cast<IndexType>(r_old); ++t) {
                IndexType i_remove = first;
                IndexType j_remove = last;
                IndexType k_remove = j_remove - knot_index + 1;

                while (j_remove - i_remove > t) {
                    if (i_remove < control_point_index) {
                        const double alpha =
                            (knot_old_b - rKnotsRefined[i_remove - 1])
                            / (knot_old_a - rKnotsRefined[i_remove - 1]);

                        const array_1d<double, 3> cp_coordinates =
                            alpha * rPointsRefined[i_remove]
                            + (1.0 - alpha) * rPointsRefined[i_remove - 1];

                        rPointsRefined(i_remove) =
                            Kratos::make_intrusive<NodeType>(0, cp_coordinates);

                        rWeightsRefined[i_remove] =
                            alpha * rWeightsRefined[i_remove]
                            + (1.0 - alpha) * rWeightsRefined[i_remove - 1];
                    }

                    if (j_remove >= left_bezier) {
                        if ((j_remove - t) <= (knot_index - degree + r_old)) {
                            const double gamma =
                                (knot_old_b - rKnotsRefined[j_remove - t - 1])
                                / (knot_old_b - knot_old_a);

                            const array_1d<double, 3> cp_coordinates =
                                gamma * elevate_bezier_points[k_remove]
                                + (1.0 - gamma) * elevate_bezier_points[k_remove + 1];

                            elevate_bezier_points(k_remove) =
                                Kratos::make_intrusive<NodeType>(0, cp_coordinates);

                            elevate_bezier_points_weights[k_remove] =
                                gamma * elevate_bezier_points_weights[k_remove]
                                + (1.0 - gamma)
                                * elevate_bezier_points_weights[k_remove + 1];
                        } else {
                            const array_1d<double, 3> cp_coordinates =
                                beta * elevate_bezier_points[k_remove]
                                + (1.0 - beta) * elevate_bezier_points[k_remove + 1];

                            elevate_bezier_points(k_remove) =
                                Kratos::make_intrusive<NodeType>(0, cp_coordinates);

                            elevate_bezier_points_weights[k_remove] =
                                beta * elevate_bezier_points_weights[k_remove]
                                + (1.0 - beta)
                                * elevate_bezier_points_weights[k_remove + 1];
                        }
                    }

                    ++i_remove;
                    --j_remove;
                    --k_remove;
                }

                --first;
                ++last;
            }
        }

        if ((a + 1) != degree_old) {
            for (IndexType i_knot = 0; i_knot < degree - r_old; ++i_knot) {
                rKnotsRefined[knot_index - 1] = knot_old_a;
                ++knot_index;
            }
        }

        for (IndexType j = left_bezier; j < right_bezier + 1; ++j) {
            rPointsRefined(control_point_index) =
                elevate_bezier_points(j);

            rWeightsRefined[control_point_index] =
                elevate_bezier_points_weights[j];

            ++control_point_index;
        }

        if (b < nb_knots_old) {
            for (IndexType j = 0; j < static_cast<IndexType>(r); ++j) {
                bezier_points(j) = next_bezier_points(j);
                bezier_points_weights[j] = next_bezier_points_weights[j];
            }

            for (IndexType j = r; j < degree_old + 1; ++j) {
                const IndexType old_index =
                    b - degree_old + j + 1;

                const array_1d<double, 3> cp_coordinates =
                    rGeometry[old_index] * weights_old[old_index];

                bezier_points(j) =
                    Kratos::make_intrusive<NodeType>(0, cp_coordinates);

                bezier_points_weights[j] =
                    weights_old[old_index];
            }

            a = b;
            ++b;
            knot_old_a = knot_old_b;
        } else {
            for (IndexType i_knot = 0; i_knot < degree; ++i_knot) {
                rKnotsRefined[knot_index - 1 + i_knot] = knot_old_b;
            }
        }
    }

    for (IndexType i = 0; i < nb_cp_refined; ++i) {
        rPointsRefined(i)->Coordinates() =
            rPointsRefined(i)->Coordinates() / rWeightsRefined[i];
    }
}

void NurbsCurveRefinementUtilities::SortAndFilter(
    std::vector<double>& rKnotsToInsert,
    const Vector& rKnotsOld)
{
    std::sort(rKnotsToInsert.begin(), rKnotsToInsert.end());

    const double lower_bound =
        std::min(rKnotsOld[0], rKnotsOld[rKnotsOld.size() - 1]);

    const double upper_bound =
        std::max(rKnotsOld[0], rKnotsOld[rKnotsOld.size() - 1]);

    auto start =
        std::lower_bound(
            rKnotsToInsert.begin(),
            rKnotsToInsert.end(),
            lower_bound);

    auto end =
        std::upper_bound(
            rKnotsToInsert.begin(),
            rKnotsToInsert.end(),
            upper_bound);

    rKnotsToInsert =
        std::vector<double>(start, end);
}

Vector NurbsCurveRefinementUtilities::CreateFullKnotVector(
    const Vector& rReducedKnots,
    const SizeType Degree)
{
    const double first = rReducedKnots[0];
    const double last = rReducedKnots[rReducedKnots.size() - 1];

    SizeType first_multiplicity = 0;
    for (IndexType i = 0; i < rReducedKnots.size(); ++i) {
        if (std::abs(rReducedKnots[i] - first) < 1.0e-12) {
            ++first_multiplicity;
        } else {
            break;
        }
    }

    SizeType last_multiplicity = 0;
    for (int i = static_cast<int>(rReducedKnots.size()) - 1; i >= 0; --i) {
        if (std::abs(rReducedKnots[i] - last) < 1.0e-12) {
            ++last_multiplicity;
        } else {
            break;
        }
    }

    const SizeType required_multiplicity = Degree + 1;

    const SizeType add_first =
        required_multiplicity > first_multiplicity
        ? required_multiplicity - first_multiplicity
        : 0;

    const SizeType add_last =
        required_multiplicity > last_multiplicity
        ? required_multiplicity - last_multiplicity
        : 0;

    Vector full_knots(rReducedKnots.size() + add_first + add_last);

    IndexType index = 0;

    for (IndexType i = 0; i < add_first; ++i) {
        full_knots[index++] = first;
    }

    for (IndexType i = 0; i < rReducedKnots.size(); ++i) {
        full_knots[index++] = rReducedKnots[i];
    }

    for (IndexType i = 0; i < add_last; ++i) {
        full_knots[index++] = last;
    }

    return full_knots;
}

Vector NurbsCurveRefinementUtilities::CreateReducedKnotVector(
    const Vector& rFullKnots,
    const SizeType Degree)
{
    const double first = rFullKnots[0];
    const double last = rFullKnots[rFullKnots.size() - 1];

    const SizeType remove_each_side = 1;

    Vector reduced_knots(rFullKnots.size() - 2 * remove_each_side);

    IndexType index = 0;

    for (IndexType i = remove_each_side; i < rFullKnots.size() - remove_each_side; ++i) {
        reduced_knots[index++] = rFullKnots[i];
    }

    return reduced_knots;
}

} // namespace Kratos