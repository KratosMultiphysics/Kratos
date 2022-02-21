//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "nurbs_surface_refinement_utilities.h"

namespace Kratos {

    ///@name Operations
    ///@{

    void NurbsSurfaceRefinementUtilities::KnotRefinementU(
        NurbsSurfaceGeometryType& rGeometry,
        std::vector<double>& rKnotsUToInsert,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsURefined,
        Vector& rWeightsRefined)
    {
        const Vector& knots_u_old = rGeometry.KnotsU();

        SortAndFilter(rKnotsUToInsert, knots_u_old);
        if (rKnotsUToInsert.size() > 0) {
            const SizeType nb_knots_u_to_insert = rKnotsUToInsert.size();

            const SizeType degree_u = rGeometry.PolynomialDegree(0);

            const SizeType nb_cp_u_old = rGeometry.PointsNumberInDirection(0);
            const SizeType nb_cp_v = rGeometry.PointsNumberInDirection(1);

            Vector weights_old = rGeometry.Weights();
            if (weights_old.size() != rGeometry.size()) {
                weights_old.resize(rGeometry.size());
                std::fill(weights_old.begin(), weights_old.end(), 1.0);
            }

            const SizeType nb_knots_u_old = knots_u_old.size();

            const SizeType a = NurbsUtilities::GetUpperSpan(degree_u, knots_u_old, knots_u_old[0]);
            const SizeType b = NurbsUtilities::GetUpperSpan(degree_u, knots_u_old, knots_u_old[knots_u_old.size() - 1]);

            const SizeType nb_cp_u_refined = nb_cp_u_old + nb_knots_u_to_insert;
            const SizeType nb_knots_u_refined = nb_knots_u_old + nb_knots_u_to_insert;

            // Resize reference variables
            if (rKnotsURefined.size() != nb_knots_u_refined) {
                rKnotsURefined.resize(nb_knots_u_refined);
            }
            rKnotsURefined = ZeroVector(nb_knots_u_refined);
            if (rWeightsRefined.size() != nb_cp_u_refined * nb_cp_v) {
                rWeightsRefined.resize(nb_cp_u_refined * nb_cp_v);
            }
            rWeightsRefined = ZeroVector(nb_cp_u_refined * nb_cp_v);
            if (rPointsRefined.size() != nb_cp_u_refined * nb_cp_v) {
                rPointsRefined.resize(nb_cp_u_refined * nb_cp_v);
            }

            // Create new node vector
            for (IndexType i = 0; i < a - degree_u + 2; ++i) {
                for (IndexType m = 0; m < nb_cp_v; ++m) {
                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_refined, nb_cp_v, i, m);
                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_old, nb_cp_v, i, m);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                }
            }

            for (IndexType i = b; i < nb_cp_u_old; ++i) {
                for (IndexType m = 0; m < nb_cp_v; ++m) {
                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_refined, nb_cp_v, i + nb_knots_u_to_insert, m);
                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_old, nb_cp_v, i, m);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                }
            }

            // Create new knot span vector
            for (IndexType i = 0; i < a + 1; ++i) {
                rKnotsURefined[i] = knots_u_old[i];
            }
            for (IndexType i = b + degree_u - 1; i < nb_knots_u_old; ++i) {
                rKnotsURefined[i + nb_knots_u_to_insert] = knots_u_old[i];
            }

            const IndexType r = nb_knots_u_to_insert - 1;

            IndexType i = b + 2 + degree_u - 1;
            IndexType k = b + 2 + degree_u + r;

            for (int j = r; j >= 0; --j) {
                while (rKnotsUToInsert[j] <= knots_u_old[i - 1] && i > a + 1) {
                    for (IndexType m = 0; m < nb_cp_v; ++m) {
                        IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u_refined, nb_cp_v, k - degree_u - 1, m);
                        IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u_old, nb_cp_v, i - degree_u - 1, m);

                        const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                        rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                        rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                    }

                    rKnotsURefined[k - 1] = knots_u_old[i - 1];

                    k -= 1;
                    i -= 1;
                }

                for (IndexType m = 0; m < nb_cp_v; m++) {
                    IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_refined, nb_cp_v, k - degree_u, m);
                    IndexType cp_index_refined_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_refined, nb_cp_v, k - degree_u - 1, m);

                    rPointsRefined(cp_index_refined_after) = rPointsRefined(cp_index_refined_before);
                    rWeightsRefined[cp_index_refined_after] = rWeightsRefined[cp_index_refined_before];
                }

                for (IndexType l = 1; l < degree_u + 1; ++l) {
                    const IndexType index = k - degree_u + l;
                    auto alpha = rKnotsURefined[k + l - 1] - rKnotsUToInsert[j];

                    if (std::abs(alpha) < 1e-7) {
                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u_refined, nb_cp_v, index, m);
                            IndexType cp_index_refined_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u_refined, nb_cp_v, index - 1, m);

                            rPointsRefined(cp_index_refined_after) = rPointsRefined(cp_index_refined_before);
                            rWeightsRefined[cp_index_refined_after] = rWeightsRefined[cp_index_refined_before];
                        }
                    }
                    else {
                        alpha = alpha / (rKnotsURefined[k + l - 1] - knots_u_old[i + l - degree_u - 1]);
                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u_refined, nb_cp_v, index, m);
                            IndexType cp_index_refined_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u_refined, nb_cp_v, index - 1, m);

                            const array_1d<double, 3> cp_coordinates = alpha * rPointsRefined[cp_index_refined_after] + (1.0 - alpha) * rPointsRefined[cp_index_refined_before];

                            rPointsRefined(cp_index_refined_after) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            rWeightsRefined[cp_index_refined_after] = rWeightsRefined[cp_index_refined_after] * alpha + rWeightsRefined[cp_index_refined_before] * (1 - alpha);
                        }
                    }
                }
                rKnotsURefined[k - 1] = rKnotsUToInsert[j];

                k -= 1;
            }

            for (IndexType i = 0; i < nb_cp_u_refined * nb_cp_v; ++i) {
                rPointsRefined(i)->Coordinates() = rPointsRefined(i)->Coordinates() / rWeightsRefined[i];
            }
        }
        else {
            KRATOS_WARNING("::NurbsSurfaceRefinementUtilities::KnotRefinementU") << "No refinement applied as knot vector is empty. "
                << "Possible errors: knots to insert are not within the range of knot_vector_u: " << knots_u_old << std::endl;
        }
    }

    void NurbsSurfaceRefinementUtilities::KnotRefinementV(
        NurbsSurfaceGeometryType& rGeometry,
        std::vector<double>& rKnotsVToInsert,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsVRefined,
        Vector& rWeightsRefined)
    {
        const Vector& knots_v_old = rGeometry.KnotsV();

        SortAndFilter(rKnotsVToInsert, knots_v_old);
        if (rKnotsVToInsert.size() > 0) {
            const SizeType nb_knots_v_to_insert = rKnotsVToInsert.size();

            const SizeType degree_v = rGeometry.PolynomialDegree(1);

            const SizeType nb_cp_u = rGeometry.PointsNumberInDirection(0);
            const SizeType nb_cp_v_old = rGeometry.PointsNumberInDirection(1);

            Vector weights_old = rGeometry.Weights();
            if (weights_old.size() != rGeometry.size()) {
                weights_old.resize(rGeometry.size());
                std::fill(weights_old.begin(), weights_old.end(), 1.0);
            }

            const SizeType nb_knots_v_old = knots_v_old.size();

            const SizeType a = NurbsUtilities::GetUpperSpan(degree_v, knots_v_old, knots_v_old[0]);
            const SizeType b = NurbsUtilities::GetUpperSpan(degree_v, knots_v_old, knots_v_old[knots_v_old.size() - 1]);

            const SizeType nb_cp_v_refined = nb_cp_v_old + nb_knots_v_to_insert;
            const SizeType nb_knots_v_refined = nb_knots_v_old + nb_knots_v_to_insert;

            // Resize reference variables
            if (rKnotsVRefined.size() != nb_knots_v_refined) {
                rKnotsVRefined.resize(nb_knots_v_refined);
            }
            rKnotsVRefined = ZeroVector(nb_knots_v_refined);
            if (rWeightsRefined.size() != nb_cp_v_refined * nb_cp_u) {
                rWeightsRefined.resize(nb_cp_v_refined * nb_cp_u);
            }
            rWeightsRefined = ZeroVector(nb_cp_v_refined * nb_cp_u);
            if (rPointsRefined.size() != nb_cp_v_refined * nb_cp_u) {
                rPointsRefined.resize(nb_cp_v_refined * nb_cp_u);
            }

            for (IndexType i = 0; i < a - degree_v + 2; ++i) {
                for (IndexType m = 0; m < nb_cp_u; ++m) {
                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_refined, m, i);
                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_old, m, i);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                }
            }

            for (IndexType i = b; i < nb_cp_v_old; ++i) {
                for (IndexType m = 0; m < nb_cp_u; ++m) {
                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_refined, m, i + nb_knots_v_to_insert);
                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_old, m, i);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                }
            }

            // Create new knot span vector
            for (IndexType i = 0; i < a + 1; i++) {
                rKnotsVRefined[i] = knots_v_old[i];
            }
            for (IndexType i = b + degree_v - 1; i < nb_knots_v_old; i++) {
                rKnotsVRefined[i + nb_knots_v_to_insert] = knots_v_old[i];
            }

            const IndexType r = nb_knots_v_to_insert - 1;

            IndexType i = b + 2 + degree_v - 1;
            IndexType k = b + 2 + degree_v + r;

            for (int j = r; j >= 0; --j) {
                while ((rKnotsVToInsert[j] <= knots_v_old[i - 1]) && (i > a + 1)) {
                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                        IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u, nb_cp_v_refined, m, k - degree_v - 1);
                        IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u, nb_cp_v_old, m, i - degree_v - 1);

                        const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                        rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                        rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                    }
                    rKnotsVRefined[k - 1] = knots_v_old[i - 1];

                    k -= 1;
                    i -= 1;
                }

                for (IndexType m = 0; m < nb_cp_u; m++) {
                    IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_refined, m, k - degree_v);
                    IndexType cp_index_refined_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_refined, m, k - degree_v - 1);

                    rPointsRefined(cp_index_refined_after) = rPointsRefined(cp_index_refined_before);
                    rWeightsRefined[cp_index_refined_after] = rWeightsRefined[cp_index_refined_before];
                }

                for (IndexType l = 1; l < degree_v + 1; l++) {
                    const IndexType index = k - degree_v + l;
                    auto alpha = rKnotsVRefined[k + l - 1] - rKnotsVToInsert[j];

                    if (std::abs(alpha) < 1e-7) {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v_refined, m, index);
                            IndexType cp_index_refined_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v_refined, m, index - 1);

                            rPointsRefined(cp_index_refined_after) = rPointsRefined(cp_index_refined_before);
                            rWeightsRefined[cp_index_refined_after] = rWeightsRefined[cp_index_refined_before];
                        }
                    }
                    else {
                        alpha = alpha / (rKnotsVRefined[k + l - 1] - knots_v_old[i + l - degree_v - 1]);
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v_refined, m, index);
                            IndexType cp_index_refined_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v_refined, m, index - 1);

                            const array_1d<double, 3> cp_coordinates = alpha * rPointsRefined[cp_index_refined_after] + (1.0 - alpha) * rPointsRefined[cp_index_refined_before];

                            rPointsRefined(cp_index_refined_after) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            rWeightsRefined[cp_index_refined_after] = rWeightsRefined[cp_index_refined_after] * alpha + rWeightsRefined[cp_index_refined_before] * (1 - alpha);
                        }
                    }
                }
                rKnotsVRefined[k - 1] = rKnotsVToInsert[j];

                k -= 1;
            }

            for (IndexType i = 0; i < nb_cp_v_refined * nb_cp_u; ++i) {
                rPointsRefined(i)->Coordinates() = rPointsRefined(i)->Coordinates() / rWeightsRefined[i];
            }
        }
        else {
            KRATOS_WARNING("::NurbsSurfaceRefinementUtilities::KnotRefinementV") << "No refinement applied as knot vector is empty. "
                << "Possible errors: knots to insert are not within the range of knot_vector_v: " << knots_v_old << std::endl;
        }
    }

    void NurbsSurfaceRefinementUtilities::DegreeElevationU(
        NurbsSurfaceGeometryType& rGeometry,
        SizeType& rDegreeUToElevate,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsURefined,
        Vector& rWeightsRefined)
    {
        const Vector& knots_u_old = rGeometry.KnotsU();

        if (rDegreeUToElevate > 0) {
            const SizeType degree_u_old = rGeometry.PolynomialDegree(0);
            const SizeType degree_u = degree_u_old + rDegreeUToElevate;

            const SizeType nb_cp_u_old = rGeometry.PointsNumberInDirection(0);
            const SizeType nb_cp_v = rGeometry.PointsNumberInDirection(1);

            Vector weights_old = rGeometry.Weights();
            if (weights_old.size() != rGeometry.size()) {
                weights_old.resize(rGeometry.size());
                std::fill(weights_old.begin(), weights_old.end(), 1.0);
            }

            const SizeType nb_knots_u_old = knots_u_old.size();
            SizeType number_of_non_zero_spans_u = 0;
            for (IndexType i = 0; i < nb_knots_u_old - 1; i++)
            {
                if (knots_u_old[i] != knots_u_old[i + 1])
                {
                    number_of_non_zero_spans_u = number_of_non_zero_spans_u + 1; 
                } 
            }

            const SizeType nb_knots_u_refined = nb_knots_u_old + rDegreeUToElevate * (number_of_non_zero_spans_u + 1);
            const SizeType nb_cp_u_refined = nb_knots_u_refined - degree_u + 1;

            // Resize reference variables
            if (rKnotsURefined.size() != nb_knots_u_refined) {
                rKnotsURefined.resize(nb_knots_u_refined);
            }
            rKnotsURefined = ZeroVector(nb_knots_u_refined);
            if (rWeightsRefined.size() != nb_cp_u_refined * nb_cp_v) {
                rWeightsRefined.resize(nb_cp_u_refined * nb_cp_v);
            }
            rWeightsRefined = ZeroVector(nb_cp_u_refined * nb_cp_v);
            if (rPointsRefined.size() != nb_cp_u_refined * nb_cp_v) {
                rPointsRefined.resize(nb_cp_u_refined * nb_cp_v);
            }

            // Initialize additional variables
            Matrix bezier_alphas = ZeroMatrix(degree_u + 1, degree_u_old + 1);
            PointerVector<NodeType> bezier_points;
            if (bezier_points.size() != (degree_u_old + 1) * nb_cp_v) {
                bezier_points.resize((degree_u_old + 1) * nb_cp_v);
            }
            Vector bezier_points_weights = ZeroVector((degree_u_old + 1) * nb_cp_v);
            PointerVector<NodeType> elevate_bezier_points;
            if (elevate_bezier_points.size() != (degree_u + 1) * nb_cp_v) {
                elevate_bezier_points.resize((degree_u + 1) * nb_cp_v);
            }
            Vector elevate_bezier_points_weights = ZeroVector((degree_u + 1) * nb_cp_v);
            PointerVector<NodeType> next_bezier_points;
            if (next_bezier_points.size() != (degree_u_old - 1) * nb_cp_v) {
                next_bezier_points.resize((degree_u_old - 1) * nb_cp_v);
            }
            Vector next_bezier_points_weights = ZeroVector((degree_u_old - 1) * nb_cp_v);

            // Compute Bézier degree elevation coefficients
            bezier_alphas(0, 0) = 1.0;
            bezier_alphas(degree_u, degree_u_old) = 1.0;

            for(IndexType i = 1; i < (degree_u / 2) + 1; ++i)
            {
                double inv = 1.0 / NurbsUtilities::GetBinomCoefficient(degree_u, i);
                SizeType min_degree = std::min(degree_u_old, i);
                int index = std::max(0, static_cast<int>(i - rDegreeUToElevate));

                for(IndexType j = index; j < min_degree + 1; ++j)
                {
                    bezier_alphas(i, j) = inv * NurbsUtilities::GetBinomCoefficient(degree_u_old, j) * NurbsUtilities::GetBinomCoefficient(rDegreeUToElevate, i - j);
                    bezier_alphas(degree_u - i, degree_u_old - j) = bezier_alphas(i, j);
                }
            }

            // Initialize coefficients
            SizeType a = degree_u_old - 1;
            SizeType b = degree_u_old;
            int r = -1;
            IndexType knot_index = degree_u + 1;
            IndexType control_point_index = 1;
            double knot_u_old_a = knots_u_old[a];

            // Create new knot span vector
            for (IndexType i = 0; i < degree_u; ++i) {
                rKnotsURefined[i] = knot_u_old_a;
            }

            // Create new node vector
            for (IndexType m = 0; m < nb_cp_v; ++m) {
                IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                    nb_cp_u_refined, nb_cp_v, 0, m);
                IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                    nb_cp_u_old, nb_cp_v, 0, m);

                const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
            }

            // Initilaize first Bézier segment
            for (IndexType i = 0; i < degree_u_old + 1; ++i) {
                for (IndexType m = 0; m < nb_cp_v; ++m) {
                    IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        (degree_u_old + 1), nb_cp_v, i, m);
                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_old, nb_cp_v, i, m);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                }
            }

            while (b < nb_knots_u_old)
            {
                IndexType i = b;
                while (b < (nb_knots_u_old - 1) && knots_u_old[b] == knots_u_old[b + 1])
                {
                    b++;
                }
                if (b + 1 == nb_knots_u_old)
                {
                    b++;
                } 
                IndexType mult = b - i + 1;
                double knot_u_old_b = knots_u_old[i];
                int r_old = r;
                r = degree_u_old - mult;

                IndexType left_bezier, right_bezier;
                if (r_old > 0)
                {
                    left_bezier = (r_old + 2) / 2; 
                }
                else
                {
                    left_bezier = 1;
                }
                if (r > 0)
                {
                    right_bezier = degree_u - ((r + 1) / 2);  
                }
                else
                {
                    right_bezier = degree_u;
                }
                
                // Insert knot to get Bezier segment
                if (r > 0)
                {
                    Vector alphas = ZeroVector(degree_u_old - 1);
                    for (IndexType k = degree_u_old; k > mult; --k)
                    {
                        alphas[k - mult - 1] = (knot_u_old_b - knot_u_old_a)/(knots_u_old[a + k] - knot_u_old_a); 
                    }

                    for (IndexType j = 1; j < static_cast<IndexType>(r + 1); ++j) 
                    {
                        int save = r - j;
                        IndexType s = mult + j;

                        for (IndexType k = degree_u_old; k >= s; --k)
                        {
                            for (IndexType m = 0; m < nb_cp_v; ++m) {
                                IndexType cp_index_bpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u_old + 1), nb_cp_v, k, m);
                                IndexType cp_index_bpts_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u_old + 1), nb_cp_v, k - 1, m);

                                const array_1d<double, 3> cp_coordinates = alphas[k - s] * bezier_points[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points[cp_index_bpts_before];

                                bezier_points(cp_index_bpts_after) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                bezier_points_weights[cp_index_bpts_after] = alphas[k - s] * bezier_points_weights[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points_weights[cp_index_bpts_before];
                            }
                        }

                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u_old + 1), nb_cp_v, degree_u_old, m);
                            IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u_old - 1), nb_cp_v, save, m);

                            next_bezier_points(cp_index_nextbpts) = bezier_points(cp_index_bpts);
                            next_bezier_points_weights[cp_index_nextbpts] = bezier_points_weights[cp_index_bpts];
                        }
                    }
                }

                // Degree elevate Bezier
                for (IndexType i = left_bezier; i < degree_u + 1; ++i)
                {
                    for (IndexType m = 0; m < nb_cp_v; ++m) {
                        IndexType cp_index_ebpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            (degree_u + 1), nb_cp_v, i, m);

                        const array_1d<double, 3> cp_coordinates = ZeroVector(3);

                        elevate_bezier_points(cp_index_ebpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                        elevate_bezier_points_weights[cp_index_ebpts] = 0.0;
                    }

                    SizeType min_degree = std::min(degree_u_old, i);
                    int index = std::max(0, static_cast<int>(i - rDegreeUToElevate));

                    for (IndexType j = index; j < min_degree + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u_old + 1), nb_cp_v, j, m);
                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u + 1), nb_cp_v, i, m);

                            const array_1d<double, 3> cp_coordinates = elevate_bezier_points[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points[cp_index_bpts];

                            elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            elevate_bezier_points_weights[cp_index_etbpts] = elevate_bezier_points_weights[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points_weights[cp_index_bpts];
                        }
                    }
                }

                // Knot removal 
                if (r_old > 1)
                {
                    IndexType first = knot_index - 2;
                    IndexType last = knot_index;

                    double beta = (knot_u_old_b - rKnotsURefined[knot_index - 2]) / (knot_u_old_b - knot_u_old_a);

                    for (IndexType t = 1; t < static_cast<IndexType>(r_old); ++t)
                    {
                        IndexType i = first;
                        IndexType j = last;
                        IndexType k = j - knot_index + 1;

                        while (j - i > t)
                        {
                            if (i < control_point_index)
                            {
                                double alpha = (knot_u_old_b - rKnotsURefined[i - 1]) / (knot_u_old_a - rKnotsURefined[i - 1]);
                                for (IndexType m = 0; m < nb_cp_v; ++m) {
                                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        nb_cp_u_refined, nb_cp_v, i, m);
                                    IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        nb_cp_u_refined, nb_cp_v, i - 1, m);

                                    const array_1d<double, 3> cp_coordinates = alpha * rPointsRefined[cp_index_refined] + (1.0 - alpha) * rPointsRefined[cp_index_refined_before];

                                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                    rWeightsRefined[cp_index_refined] = rWeightsRefined[cp_index_refined] * alpha + rWeightsRefined[cp_index_refined_before] * (1 - alpha);
                                }
                            }
                            if (j >= left_bezier)
                            {
                                if ((j - t) <= (knot_index - degree_u + r_old))
                                {
                                    double gamma = (knot_u_old_b - rKnotsURefined[j - t - 1]) / (knot_u_old_b - knot_u_old_a);
                                    for (IndexType m = 0; m < nb_cp_v; ++m) {
                                        IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            (degree_u + 1), nb_cp_v, k, m);
                                        IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            (degree_u + 1), nb_cp_v, k + 1, m);

                                        const array_1d<double, 3> cp_coordinates = gamma * elevate_bezier_points[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points[cp_index_etbpts_after];

                                        elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                        elevate_bezier_points_weights[cp_index_etbpts] = gamma * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                    }

                                }
                                else
                                {
                                    for (IndexType m = 0; m < nb_cp_v; ++m) {
                                        IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            (degree_u + 1), nb_cp_v, k, m);
                                        IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            (degree_u + 1), nb_cp_v, k + 1, m);

                                        const array_1d<double, 3> cp_coordinates = beta * elevate_bezier_points[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points[cp_index_etbpts_after];

                                        elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                        elevate_bezier_points_weights[cp_index_etbpts] = beta * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                    }
                                }
                            }
                            ++i;
                            j = j -1;
                            k = k -1;
                        }
                        --first;
                        ++last;
                    }
                } 

                // Create next knot span vector
                if ((a + 1) != degree_u_old) 
                {
                    for (IndexType i = 0; i < (degree_u - r_old); ++i)
                    {
                        rKnotsURefined[knot_index - 1] = knot_u_old_a;
                        knot_index = knot_index + 1;
                    }
                }

                // Create next node vector
                for (IndexType j = left_bezier; j < right_bezier + 1; ++j)
                {
                    for (IndexType m = 0; m < nb_cp_v; ++m) {
                        IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u_refined, nb_cp_v, control_point_index, m);
                        IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            (degree_u + 1), nb_cp_v, j, m);

                        rPointsRefined(cp_index) = elevate_bezier_points(cp_index_etbpts);
                        rWeightsRefined[cp_index] = elevate_bezier_points_weights[cp_index_etbpts];
                    }
                    control_point_index = control_point_index + 1;
                }

                // Setup for the next pass through loop
                if (b < nb_knots_u_old)
                {
                    for (IndexType j = 0; j < static_cast<IndexType>(r); ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            IndexType cp_index_btps = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u_old + 1), nb_cp_v, j, m);
                            IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u_old - 1), nb_cp_v, j, m);

                            bezier_points(cp_index_btps) = next_bezier_points(cp_index_nextbpts);
                            bezier_points_weights[cp_index_btps] = next_bezier_points_weights[cp_index_nextbpts];
                        }
                    }
                    for (IndexType j = r; j < degree_u_old + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u_old + 1), nb_cp_v, j, m);
                            IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u_old, nb_cp_v, (b - degree_u_old + j + 1), m);

                            const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                            bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                        }
                    }

                    a = b;
                    b = b + 1;
                    knot_u_old_a = knot_u_old_b;
                }
                else
                {
                    for (IndexType i = 0; i < degree_u; ++i)
                    {
                        rKnotsURefined[knot_index - 1 + i] = knot_u_old_b;
                    }
                }
            }
            
            for (IndexType i = 0; i < nb_cp_u_refined * nb_cp_v; ++i) {
                rPointsRefined(i)->Coordinates() = rPointsRefined(i)->Coordinates() / rWeightsRefined[i];
            }
        }
        else {
            KRATOS_WARNING("::NurbsSurfaceRefinementUtilities::DegreeElevationU") << "No elevation applied as degree is zero. "
                << std::endl;
        }
    }

    void NurbsSurfaceRefinementUtilities::DegreeElevationV(
        NurbsSurfaceGeometryType& rGeometry,
        SizeType& rDegreeVToElevate,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsVRefined,
        Vector& rWeightsRefined)
    {
        const Vector& knots_v_old = rGeometry.KnotsV();

        if (rDegreeVToElevate > 0) {
            const SizeType degree_v_old = rGeometry.PolynomialDegree(1);
            const SizeType degree_v = degree_v_old + rDegreeVToElevate;

            const SizeType nb_cp_u = rGeometry.PointsNumberInDirection(0);
            const SizeType nb_cp_v_old = rGeometry.PointsNumberInDirection(1);

            Vector weights_old = rGeometry.Weights();
            if (weights_old.size() != rGeometry.size()) {
                weights_old.resize(rGeometry.size());
                std::fill(weights_old.begin(), weights_old.end(), 1.0);
            }

            const SizeType nb_knots_v_old = knots_v_old.size();
            SizeType number_of_non_zero_spans_v = 0;
            for (IndexType i = 0; i < nb_knots_v_old - 1; i++)
            {
                if (knots_v_old[i] != knots_v_old[i + 1])
                {
                    number_of_non_zero_spans_v = number_of_non_zero_spans_v + 1; 
                } 
            }

            const SizeType nb_knots_v_refined = nb_knots_v_old + rDegreeVToElevate * (number_of_non_zero_spans_v + 1);
            const SizeType nb_cp_v_refined = nb_knots_v_refined - degree_v + 1;

            // Resize reference variables
            if (rKnotsVRefined.size() != nb_knots_v_refined) {
                rKnotsVRefined.resize(nb_knots_v_refined);
            }
            rKnotsVRefined = ZeroVector(nb_knots_v_refined);
            if (rWeightsRefined.size() != nb_cp_v_refined * nb_cp_u) {
                rWeightsRefined.resize(nb_cp_v_refined * nb_cp_u);
            }
            rWeightsRefined = ZeroVector(nb_cp_v_refined * nb_cp_u);
            if (rPointsRefined.size() != nb_cp_v_refined * nb_cp_u) {
                rPointsRefined.resize(nb_cp_v_refined * nb_cp_u);
            }

            // Initialize additional variables
            Matrix bezier_alphas = ZeroMatrix(degree_v + 1, degree_v_old + 1);
            PointerVector<NodeType> bezier_points;
            if (bezier_points.size() != (degree_v_old + 1) * nb_cp_u) {
                bezier_points.resize((degree_v_old + 1) * nb_cp_u);
            }
            Vector bezier_points_weights = ZeroVector((degree_v_old + 1) * nb_cp_u);
            PointerVector<NodeType> elevate_bezier_points;
            if (elevate_bezier_points.size() != (degree_v + 1) * nb_cp_u) {
                elevate_bezier_points.resize((degree_v + 1) * nb_cp_u);
            }
            Vector elevate_bezier_points_weights = ZeroVector((degree_v + 1) * nb_cp_u);
            PointerVector<NodeType> next_bezier_points;
            if (next_bezier_points.size() != (degree_v_old - 1) * nb_cp_u) {
                next_bezier_points.resize((degree_v_old - 1) * nb_cp_u);
            }
            Vector next_bezier_points_weights = ZeroVector((degree_v_old - 1) * nb_cp_u);

            // Compute Bézier degree elevation coefficients
            bezier_alphas(0, 0) = 1.0;
            bezier_alphas(degree_v, degree_v_old) = 1.0;

            for(IndexType i = 1; i < (degree_v / 2) + 1; ++i)
            {
                double inv = 1.0 / NurbsUtilities::GetBinomCoefficient(degree_v, i);
                SizeType min_degree = std::min(degree_v_old, i);
                int index = std::max(0, static_cast<int>(i - rDegreeVToElevate));

                for(IndexType j = index; j < min_degree + 1; ++j)
                {
                    bezier_alphas(i, j) = inv * NurbsUtilities::GetBinomCoefficient(degree_v_old, j) * NurbsUtilities::GetBinomCoefficient(rDegreeVToElevate, i - j);
                    bezier_alphas(degree_v - i, degree_v_old - j) = bezier_alphas(i, j);
                }
            }

            // Initialize coefficients
            SizeType a = degree_v_old - 1;
            SizeType b = degree_v_old;
            int r = -1;
            IndexType knot_index = degree_v + 1;
            IndexType control_point_index = 1;
            double knot_v_old_a = knots_v_old[a];

            // Create new knot span vector
            for (IndexType i = 0; i < degree_v; ++i) {
                rKnotsVRefined[i] = knot_v_old_a;
            }

            // Create new node vector
            for (IndexType m = 0; m < nb_cp_u; ++m) {
                IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                    nb_cp_u, nb_cp_v_refined, m, 0);
                IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                    nb_cp_u, nb_cp_v_old, m, 0);

                const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
            }

            // Initilaize first Bézier segment
            for (IndexType i = 0; i < degree_v_old + 1; ++i) {
                for (IndexType m = 0; m < nb_cp_u; ++m) {
                    IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, (degree_v_old + 1), m, i);
                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_old, m, i);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                }
            }

            while (b < nb_knots_v_old)
            {
                IndexType i = b;
                while (b < (nb_knots_v_old - 1) && knots_v_old[b] == knots_v_old[b + 1])
                {
                    b++;
                }
                if (b + 1 == nb_knots_v_old)
                {
                    b++;
                } 
                IndexType mult = b - i + 1;
                double knot_v_old_b = knots_v_old[i];
                int r_old = r;
                r = degree_v_old - mult;

                IndexType left_bezier, right_bezier;
                if (r_old > 0)
                {
                    left_bezier = (r_old + 2) / 2; 
                }
                else
                {
                    left_bezier = 1;
                }
                if (r > 0)
                {
                    right_bezier = degree_v - ((r + 1) / 2);  
                }
                else
                {
                    right_bezier = degree_v;
                }
                
                // Insert knot to get Bezier segment
                if (r > 0)
                {
                    Vector alphas = ZeroVector(degree_v_old - 1);
                    for (IndexType k = degree_v_old; k > mult; --k)
                    {
                        alphas[k - mult - 1] = (knot_v_old_b - knot_v_old_a)/(knots_v_old[a + k] - knot_v_old_a); 
                    }

                    for (IndexType j = 1; j < static_cast<IndexType>(r + 1); ++j) 
                    {
                        int save = r - j;
                        IndexType s = mult + j;

                        for (IndexType k = degree_v_old; k >= s; --k)
                        {
                            for (IndexType m = 0; m < nb_cp_u; ++m) {
                                IndexType cp_index_bpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    nb_cp_u, (degree_v_old + 1), m, k);
                                IndexType cp_index_bpts_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    nb_cp_u, (degree_v_old + 1), m, k - 1);

                                const array_1d<double, 3> cp_coordinates = alphas[k - s] * bezier_points[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points[cp_index_bpts_before];

                                bezier_points(cp_index_bpts_after) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                bezier_points_weights[cp_index_bpts_after] = alphas[k - s] * bezier_points_weights[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points_weights[cp_index_bpts_before];
                            }
                        }

                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v_old + 1), m, degree_v_old);
                            IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v_old - 1), m, save);

                            next_bezier_points(cp_index_nextbpts) = bezier_points(cp_index_bpts);
                            next_bezier_points_weights[cp_index_nextbpts] = bezier_points_weights[cp_index_bpts];
                        }
                    }
                }

                // Degree elevate Bezier
                for (IndexType i = left_bezier; i < degree_v + 1; ++i)
                {
                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                        IndexType cp_index_ebpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u, (degree_v + 1), m, i);

                        const array_1d<double, 3> cp_coordinates = ZeroVector(3);

                        elevate_bezier_points(cp_index_ebpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                        elevate_bezier_points_weights[cp_index_ebpts] = 0.0;
                    }

                    SizeType min_degree = std::min(degree_v_old, i);
                    int index = std::max(0, static_cast<int>(i - rDegreeVToElevate));

                    for (IndexType j = index; j < min_degree + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v_old + 1), m, j);
                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v + 1), m, i);

                            const array_1d<double, 3> cp_coordinates = elevate_bezier_points[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points[cp_index_bpts];

                            elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            elevate_bezier_points_weights[cp_index_etbpts] = elevate_bezier_points_weights[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points_weights[cp_index_bpts];
                        }
                    }
                }

                // Knot removal 
                if (r_old > 1)
                {
                    IndexType first = knot_index - 2;
                    IndexType last = knot_index;

                    double beta = (knot_v_old_b - rKnotsVRefined[knot_index - 2]) / (knot_v_old_b - knot_v_old_a);

                    for (IndexType t = 1; t < static_cast<IndexType>(r_old); ++t)
                    {
                        IndexType i = first;
                        IndexType j = last;
                        IndexType k = j - knot_index + 1;

                        while (j - i > t)
                        {
                            if (i < control_point_index)
                            {
                                double alpha = (knot_v_old_b - rKnotsVRefined[i - 1]) / (knot_v_old_a - rKnotsVRefined[i - 1]);
                                for (IndexType m = 0; m < nb_cp_u; ++m) {
                                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        nb_cp_u, nb_cp_v_refined, m, i);
                                    IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        nb_cp_u, nb_cp_v_refined, m, i - 1);

                                    const array_1d<double, 3> cp_coordinates = alpha * rPointsRefined[cp_index_refined] + (1.0 - alpha) * rPointsRefined[cp_index_refined_before];

                                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                    rWeightsRefined[cp_index_refined] = rWeightsRefined[cp_index_refined] * alpha + rWeightsRefined[cp_index_refined_before] * (1 - alpha);
                                }
                            }
                            if (j >= left_bezier)
                            {
                                if ((j - t) <= (knot_index - degree_v + r_old))
                                {
                                    double gamma = (knot_v_old_b - rKnotsVRefined[j - t - 1]) / (knot_v_old_b - knot_v_old_a);
                                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                                        IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u, (degree_v + 1), m, k);
                                        IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u, (degree_v + 1), m, k + 1);

                                        const array_1d<double, 3> cp_coordinates = gamma * elevate_bezier_points[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points[cp_index_etbpts_after];

                                        elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                        elevate_bezier_points_weights[cp_index_etbpts] = gamma * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                    }

                                }
                                else
                                {
                                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                                        IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u, (degree_v + 1), m, k);
                                        IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u, (degree_v + 1), m, k + 1);

                                        const array_1d<double, 3> cp_coordinates = beta * elevate_bezier_points[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points[cp_index_etbpts_after];

                                        elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                        elevate_bezier_points_weights[cp_index_etbpts] = beta * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                    }
                                }
                            }
                            ++i;
                            j = j -1;
                            k = k -1;
                        }
                        --first;
                        ++last;
                    }
                } 

                // Create next knot span vector
                if ((a + 1) != degree_v_old) 
                {
                    for (IndexType i = 0; i < (degree_v - r_old); ++i)
                    {
                        rKnotsVRefined[knot_index - 1] = knot_v_old_a;
                        knot_index = knot_index + 1;
                    }
                }

                // Create next node vector
                for (IndexType j = left_bezier; j < right_bezier + 1; ++j)
                {
                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                        IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u, nb_cp_v_refined, m, control_point_index);
                        IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u, (degree_v + 1), m, j);

                        rPointsRefined(cp_index) = elevate_bezier_points(cp_index_etbpts);
                        rWeightsRefined[cp_index] = elevate_bezier_points_weights[cp_index_etbpts];
                    }
                    control_point_index = control_point_index + 1;
                }

                // Setup for the next pass through loop
                if (b < nb_knots_v_old)
                {
                    for (IndexType j = 0; j < static_cast<IndexType>(r); ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            IndexType cp_index_btps = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v_old + 1), m, j);
                            IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v_old - 1), m, j);

                            bezier_points(cp_index_btps) = next_bezier_points(cp_index_nextbpts);
                            bezier_points_weights[cp_index_btps] = next_bezier_points_weights[cp_index_nextbpts];
                        }
                    }
                    for (IndexType j = r; j < degree_v_old + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v_old + 1), m, j);
                            IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v_old, m, (b - degree_v_old + j + 1));

                            const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                            bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                        }
                    }

                    a = b;
                    b = b + 1;
                    knot_v_old_a = knot_v_old_b;
                }
                else
                {
                    for (IndexType i = 0; i < degree_v; ++i)
                    {
                        rKnotsVRefined[knot_index - 1 + i] = knot_v_old_b;
                    }
                }
            }
            
            for (IndexType i = 0; i < nb_cp_v_refined * nb_cp_u; ++i) {
                rPointsRefined(i)->Coordinates() = rPointsRefined(i)->Coordinates() / rWeightsRefined[i];
            }
        }
        else {
            KRATOS_WARNING("::NurbsSurfaceRefinementUtilities::DegreeElevationV") << "No elevation applied as degree is zero. "
                << std::endl;
        }
    }

    ///@}

    void NurbsSurfaceRefinementUtilities::SortAndFilter(
        std::vector<double>& rKnotsToInsert,
        const Vector& rKnotsOld)
    {
        std::sort(rKnotsToInsert.begin(), rKnotsToInsert.end());

        const double lower_bound = std::min(rKnotsOld[0], rKnotsOld[rKnotsOld.size() - 1]);
        const double upper_bound = std::max(rKnotsOld[0], rKnotsOld[rKnotsOld.size() - 1]);

        auto start = std::lower_bound(rKnotsToInsert.begin(), rKnotsToInsert.end(), lower_bound);
        auto end = std::upper_bound(rKnotsToInsert.begin(), rKnotsToInsert.end(), upper_bound);

        rKnotsToInsert = std::vector<double>(start, end);
    }

} // End namespace kratos