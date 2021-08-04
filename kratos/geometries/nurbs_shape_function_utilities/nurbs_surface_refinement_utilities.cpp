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
        KRATOS_WATCH(rKnotsUToInsert)
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
        KRATOS_WATCH(knots_v_old)
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