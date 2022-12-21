//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_NURBS_SURFACE_REFINEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_NURBS_SURFACE_REFINEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/nurbs_surface_geometry.h"
#include "nurbs_utilities.h"
#include "includes/node.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) NurbsSurfaceRefinementUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node<3> NodeType;

    typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceGeometryType;
    typedef typename NurbsSurfaceGeometryType::Pointer NurbsSurfaceGeometryPointerType;

    ///@}
    ///@name Operations
    ///@{

    /*
     * @brief Refines the u-knot vector of a NurbsSurfaceGeometry.
     * @details This function adopts the surface knot-refinement algorithm from
     *          Piegl 1995, Algorithm A5.5.
     *
     * Portet from carat++ (https://www.bgu.tum.de/st/software/forschung/carat/)
     * and ANurbs (https://github.com/oberbichler/ANurbs)
     *
     * @param rGeometry surface to be refined.
     * @param rKnotsUToInsert Knots to be inserted.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsURefined the new knot vector.
     * @param rWeightsRefined the new weight vector.
     */
    template<class TNurbsSurfaceGeometryType>
    static void KnotRefinementU(
        TNurbsSurfaceGeometryType & rGeometry,
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

    /*
     * @brief Refines the v-knot vector of a NurbsSurfaceGeometry.
     * @details This function adopts the surface knot-refinement algorithm from
     *          Piegl 1995, Algorithm A5.5.
     *
     * Portet from carat++ (https://www.bgu.tum.de/st/software/forschung/carat/)
     * and ANurbs (https://github.com/oberbichler/ANurbs)
     *
     * @param rGeometry surface to be refined.
     * @param rKnotsVToInsert Knots to be inserted.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsVRefined the new knot vector.
     * @param rWeightsRefined the new weight vector.
     */
    template<class TNurbsSurfaceGeometryType>
    static void KnotRefinementV(
        TNurbsSurfaceGeometryType& rGeometry,
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

    /*
     * @brief Elevates the degree-u of a NurbsSurfaceGeometry.
     * @details This function adopts the surface degree elevation algorithm from
     *          Piegl 1995, Algorithm A5.10.
     *
     * Portet from carat++ (https://www.bgu.tum.de/st/software/forschung/carat/)
     *
     * @param rGeometry surface to be refined.
     * @param rDegreeUToElevate Degree to be elevated.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsURefined the new knot vector.
     * @param rWeightsRefined the new weight vector.
     */
    static void DegreeElevationU(
        NurbsSurfaceGeometryType& rGeometry,
        SizeType& rDegreeUToElevate,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsURefined,
        Vector& rWeightsRefined);

    /*
     * @brief Elevates the degree-v of a NurbsSurfaceGeometry.
     * @details This function adopts the surface degree elevation algorithm from
     *          Piegl 1995, Algorithm A5.10.
     *
     * Portet from carat++ (https://www.bgu.tum.de/st/software/forschung/carat/)
     *
     * @param rGeometry surface to be refined.
     * @param rDegreeVToElevate Degree to be elevated.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsVRefined the new knot vector.
     * @param rWeightsRefined the new weight vector.
     */
    static void DegreeElevationV(
        NurbsSurfaceGeometryType& rGeometry,
        SizeType& rDegreeVToElevate,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsVRefined,
        Vector& rWeightsRefined);

    ///@}
    ///@name Utilities
    ///@{

    /*
     * @brief Sorts the knot vector, which is required for the refinement.
     *        Filters the knots to insert which are within the former knot span
     *        and accordingly coordinates of the former geometry.
     *
     * @param rKnotsToInsert the knot vector to be sorted and filtered.
     * @param rKnotsOld the old knot vector, defining the limits.
     */
    static void SortAndFilter(
        std::vector<double>&rKnotsToInsert,
        const Vector& rKnotsOld);

};

} // End namespace kratos

#endif // KRATOS_NURBS_SURFACE_REFINEMENT_UTILITIES_H_INCLUDED