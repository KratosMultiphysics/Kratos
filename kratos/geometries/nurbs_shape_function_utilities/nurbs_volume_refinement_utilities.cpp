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
#include "nurbs_volume_refinement_utilities.h"

namespace Kratos {

    /**
     * @brief Refines the u-knotvector of a NurbsVolumeGeometry.
     * @details This function adopts the Volume knot-refinement algorithm from Piegl1995 (p.164 Algorithm A5.5).
     *          The algorithm is modified to suit trivariant B-Spline volumes.
     * @param rGeometry Volume to be refined.
     * @param rKnotsUToInsert Knots to be inserted.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsURefined the new knot vector.
     * @note This function does not consider weights, thus only B-Spline-Volumes can be refined.
     **/
    void NurbsVolumeRefinementUtilities::KnotRefinementU(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsUToInsert,
                            PointerVector<NodeType>& rPointsRefined, Vector& rKnotsURefined ){

        // Sort the knots which are to be inserted!
        std::sort(rKnotsUToInsert.begin(),rKnotsUToInsert.end());
        // Get current order
        const SizeType polynomial_degree_u = rGeometry.PolynomialDegreeU();

        // Get current knot information
        const Kratos::Vector& old_knots_u = rGeometry.KnotsU();

        const SizeType old_num_of_knots_u = rGeometry.NumberOfKnotsU();
        // Get current cp's information
        const SizeType old_num_of_cp_u = rGeometry.NumberOfControlPointsU();
        const SizeType old_num_of_cp_v = rGeometry.NumberOfControlPointsV();
        const SizeType old_num_of_cp_w = rGeometry.NumberOfControlPointsW();

        // Get current span's
        SizeType a = NurbsUtilities::GetLowerSpan(polynomial_degree_u, old_knots_u, rKnotsUToInsert.front());
        SizeType b = NurbsUtilities::GetLowerSpan(polynomial_degree_u, old_knots_u, rKnotsUToInsert.back()) + 1;

        SizeType r = rKnotsUToInsert.size();
        // Initialize new containers
        SizeType new_num_of_knots_u = old_num_of_knots_u + r;
        rKnotsURefined.resize(new_num_of_knots_u);

        SizeType new_num_of_cp_u = old_num_of_cp_u + r;
        rPointsRefined.resize(new_num_of_cp_u*old_num_of_cp_v*old_num_of_cp_w);

        r = r - 1;
        // Create new ordered knot vector.
        for( IndexType i = 0; i <= a; i++)
            rKnotsURefined[i] = old_knots_u[i];
        for( IndexType i = b+polynomial_degree_u; i < old_num_of_knots_u; ++i)
            rKnotsURefined[i+r+1] = old_knots_u[i];

        // Copy unaltered cp's.
        for( IndexType column=0; column < old_num_of_cp_v; ++column){
            for( IndexType depth=0; depth < old_num_of_cp_w; ++depth){
                // Copy unaltered control points.
                // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                for( IndexType i = 0; i <= a - polynomial_degree_u + 1; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                }
                // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                for( IndexType i = b; i < old_num_of_cp_u; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i+r+1, column, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                }
            }
        }

        // Find new CP's
        const SizeType p = polynomial_degree_u;
        int i = b + p - 1;
        int k = b + p + r;
        for( int j = r; j >= 0; j--){
            while( rKnotsUToInsert[j] <= old_knots_u[i] && i > static_cast<int>(a) ){
                for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                    for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                        // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                        // Also keep attention to the first argument. The left index is mapped to new_num_of_cp_u, but the right index is mapped
                        // to old_num_of_cp_u.
                        IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, k-p-1+1, column, depth);
                        IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i-p-1+1, column, depth);
                        rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                    }
                }
                rKnotsURefined[k] = old_knots_u[i];
                k--;
                i--;
            }
            for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                   // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, k-p-1+1, column, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, k-p+1, column, depth);
                    rPointsRefined(cp_index_left) =  rPointsRefined(cp_index_right);
                }
            }
            for( IndexType l=1; l <= p; ++l){
                IndexType index = k - p + l;
                double alpha = rKnotsURefined[k+l] - rKnotsUToInsert[j];
                if( std::abs(alpha) < 1e-10){
                    for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                        for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                            // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                            IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index-1+1, column, depth);
                            IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index+1, column, depth);
                            rPointsRefined(cp_index_left) = rPointsRefined(cp_index_right);
                        }
                    }
                }
                else {
                    alpha = alpha / ( rKnotsURefined[k+l] - old_knots_u[i-p+l]);
                    for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                        for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                            // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                            IndexType cp_index_minus = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index-1+1, column, depth);
                            IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index+1, column, depth);

                            double x = alpha * rPointsRefined[cp_index_minus][0] + (1.0-alpha) * rPointsRefined[cp_index][0];
                            double y = alpha * rPointsRefined[cp_index_minus][1] + (1.0-alpha) * rPointsRefined[cp_index][1];
                            double z = alpha * rPointsRefined[cp_index_minus][2] + (1.0-alpha) * rPointsRefined[cp_index][2];

                            rPointsRefined(cp_index_minus) = Kratos::make_intrusive<NodeType>(0, x, y, z);
                        }
                    }
                }
            }
            rKnotsURefined[k] = rKnotsUToInsert[j];
            k -= 1;
        }
    }


    void NurbsVolumeRefinementUtilities::KnotRefinementV(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsVToInsert,
                                        PointerVector<NodeType>& rPointsRefined, Vector& rKnotsVRefined ){

        // Sort the knots which are to be inserted!
        std::sort(rKnotsVToInsert.begin(),rKnotsVToInsert.end());
        // Get current order
        const SizeType polynomial_degree_v = rGeometry.PolynomialDegreeV();
        // Get current knot information
        const Kratos::Vector& old_knots_v = rGeometry.KnotsV();

        const SizeType old_num_of_knots_v = rGeometry.NumberOfKnotsV();
        // Get current cp's information
        const SizeType old_num_of_cp_u = rGeometry.NumberOfControlPointsU();
        const SizeType old_num_of_cp_v = rGeometry.NumberOfControlPointsV();
        const SizeType old_num_of_cp_w = rGeometry.NumberOfControlPointsW();

        // Get current span's
        SizeType a = NurbsUtilities::GetLowerSpan(polynomial_degree_v, old_knots_v, rKnotsVToInsert.front());
        SizeType b = NurbsUtilities::GetLowerSpan(polynomial_degree_v, old_knots_v, rKnotsVToInsert.back()) + 1;

        SizeType r = rKnotsVToInsert.size();
        // Initialize new containers
        SizeType new_num_of_knots_v = old_num_of_knots_v + r;
        rKnotsVRefined.resize(new_num_of_knots_v);

        SizeType new_num_of_cp_v = old_num_of_cp_v + r;
        rPointsRefined.resize(old_num_of_cp_u*new_num_of_cp_v*old_num_of_cp_w);

        r = r - 1;
        // Create new ordered knot vector.
        for( IndexType i = 0; i <= a; i++)
            rKnotsVRefined[i] = old_knots_v[i];
        for( IndexType i = b+polynomial_degree_v; i < old_num_of_knots_v; ++i)
            rKnotsVRefined[i+r+1] = old_knots_v[i];

        // Copy unaltered cp's.
        for( IndexType row=0; row < old_num_of_cp_u; ++row){
            for( IndexType depth=0; depth < old_num_of_cp_w; ++depth){
                // Copy unaltered control points.
                // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                for( IndexType i = 0; i <= a - polynomial_degree_v + 1; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, i, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, row, i, depth);
                    rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                }
                // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                for( IndexType i = b; i < old_num_of_cp_v; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, i+r+1, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, row, i, depth);
                    rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                }
            }
        }

        // Find new CP's
        const SizeType p = polynomial_degree_v;
        int i = b + p - 1;
        int k = b + p + r;
        for( int j = r; j >= 0; j--){
            while( rKnotsVToInsert[j] <= old_knots_v[i] && i > static_cast<int>(a)){
                for( IndexType row=0; row < old_num_of_cp_u; ++row) {
                    for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                        // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                        // Also keep attention to the second argument. The left index is mapped to new_num_of_cp_v, but the right index is mapped
                        // to old_num_of_cp_v.
                        IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, k-p-1+1, depth);
                        IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, row, i-p-1+1, depth);
                        rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                    }
                }
                rKnotsVRefined[k] = old_knots_v[i];
                k--;
                i--;
            }
            for( IndexType row=0; row < old_num_of_cp_u; ++row) {
                for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                   // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, k-p-1+1, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, k-p+1, depth);
                    rPointsRefined(cp_index_left) =  rPointsRefined(cp_index_right);
                }
            }
            for( IndexType l=1; l <= p; ++l){
                IndexType index = k - p + l;
                double alpha = rKnotsVRefined[k+l] - rKnotsVToInsert[j];
                if( std::abs(alpha) < 1e-10){
                    for( IndexType row=0; row < old_num_of_cp_u; ++row) {
                        for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                            // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                            IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, index-1+1, depth);
                            IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, index+1, depth);
                            rPointsRefined(cp_index_left) = rPointsRefined(cp_index_right);
                        }
                    }
                }
                else {
                    alpha = alpha / ( rKnotsVRefined[k+l] - old_knots_v[i-p+l]);
                    for( IndexType row=0; row < old_num_of_cp_u; ++row) {
                        for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                            // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                            IndexType cp_index_minus = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, index-1+1, depth);
                            IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                old_num_of_cp_u, new_num_of_cp_v, old_num_of_cp_w, row, index+1, depth);

                            double x = alpha * rPointsRefined[cp_index_minus][0] + (1.0-alpha) * rPointsRefined[cp_index][0];
                            double y = alpha * rPointsRefined[cp_index_minus][1] + (1.0-alpha) * rPointsRefined[cp_index][1];
                            double z = alpha * rPointsRefined[cp_index_minus][2] + (1.0-alpha) * rPointsRefined[cp_index][2];

                            rPointsRefined(cp_index_minus) = Kratos::make_intrusive<NodeType>(0, x, y, z);
                        }
                    }
                }

            }
            rKnotsVRefined[k] = rKnotsVToInsert[j];
            k -= 1;
        }
    }

    void NurbsVolumeRefinementUtilities::KnotRefinementW(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsWToInsert,
                                         PointerVector<NodeType>& rPointsRefined, Vector& rKnotsWRefined ){

        // Sort the knots which are to be inserted!
        std::sort(rKnotsWToInsert.begin(),rKnotsWToInsert.end());
        // Get current order
        const SizeType polynomial_degree_w = rGeometry.PolynomialDegreeW();
        // Get current knot information
        const Kratos::Vector& old_knots_w = rGeometry.KnotsW();

        const SizeType old_num_of_knots_w = rGeometry.NumberOfKnotsW();
        // Get current cp's information
        const SizeType old_num_of_cp_u = rGeometry.NumberOfControlPointsU();
        const SizeType old_num_of_cp_v = rGeometry.NumberOfControlPointsV();
        const SizeType old_num_of_cp_w = rGeometry.NumberOfControlPointsW();

        // Get current span's
        SizeType a = NurbsUtilities::GetLowerSpan(polynomial_degree_w, old_knots_w, rKnotsWToInsert.front());
        SizeType b = NurbsUtilities::GetLowerSpan(polynomial_degree_w, old_knots_w, rKnotsWToInsert.back()) + 1;

        SizeType r = rKnotsWToInsert.size();
        // Initialize new containers
        SizeType new_num_of_knots_w = old_num_of_knots_w + r;
        rKnotsWRefined.resize(new_num_of_knots_w);

        SizeType new_num_of_cp_w = old_num_of_cp_w + r;
        rPointsRefined.resize(old_num_of_cp_u*old_num_of_cp_v*new_num_of_cp_w);

        r = r - 1;
        // Create new ordered knot vector.
        for( IndexType i = 0; i <= a; i++)
            rKnotsWRefined[i] = old_knots_w[i];
        for( IndexType i = b+polynomial_degree_w; i < old_num_of_knots_w; ++i)
            rKnotsWRefined[i+r+1] = old_knots_w[i];

        // Copy unaltered cp's.
        for( IndexType row=0; row < old_num_of_cp_u; ++row){
            for( IndexType column=0; column < old_num_of_cp_v; ++column){
                // Copy unaltered control points.
                // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                for( IndexType i = 0; i <= a - polynomial_degree_w + 1; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, i);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, row, column, i);
                    rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                }
                // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                for( IndexType i = b; i < old_num_of_cp_w; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, i+r+1);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, row, column, i);
                    rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                }
            }
        }

        // Find new CP's
        const SizeType p = polynomial_degree_w;
        int i = b + p - 1;
        int k = b + p + r;
        for( int j = r; j >= 0; j--){
            while( rKnotsWToInsert[j] <= old_knots_w[i] && i > static_cast<int>(a)){
                for( IndexType row=0; row < old_num_of_cp_u; ++row) {
                    for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                        // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                        // Also keep attention to the second argument. The left index is mapped to new_num_of_cp_w, but the right index is mapped
                        // to old_num_of_cp_w.
                        IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, k-p-1+1);
                        IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, row, column, i-p-1+1);
                        rPointsRefined(cp_index_left) = Kratos::make_intrusive<NodeType>(0, rGeometry[cp_index_right]);
                    }
                }
                rKnotsWRefined[k] = old_knots_w[i];
                k--;
                i--;
            }
            for( IndexType row=0; row < old_num_of_cp_u; ++row) {
                for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                   // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, k-p-1+1);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, k-p+1);
                    rPointsRefined(cp_index_left) =  rPointsRefined(cp_index_right);
                }
            }
            for( IndexType l=1; l <= p; ++l){
                IndexType index = k - p + l;
                double alpha = rKnotsWRefined[k+l] - rKnotsWToInsert[j];
                if( std::abs(alpha) < 1e-10){
                    for( IndexType row=0; row < old_num_of_cp_u; ++row) {
                        for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                            // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                            IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, index-1+1);
                            IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, index+1);
                            rPointsRefined(cp_index_left) = rPointsRefined(cp_index_right);
                        }
                    }
                }
                else {
                    alpha = alpha / ( rKnotsWRefined[k+l] - old_knots_w[i-p+l]);
                    for( IndexType row=0; row < old_num_of_cp_u; ++row) {
                        for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                            // Note: The "+1" accounts for the fact that first and last knots only appear "p"-times (and not "p+1"-times).
                            IndexType cp_index_minus = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, index-1+1);
                            IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                old_num_of_cp_u, old_num_of_cp_v, new_num_of_cp_w, row, column, index+1);

                            double x = alpha * rPointsRefined[cp_index_minus][0] + (1.0-alpha) * rPointsRefined[cp_index][0];
                            double y = alpha * rPointsRefined[cp_index_minus][1] + (1.0-alpha) * rPointsRefined[cp_index][1];
                            double z = alpha * rPointsRefined[cp_index_minus][2] + (1.0-alpha) * rPointsRefined[cp_index][2];

                            rPointsRefined(cp_index_minus) = Kratos::make_intrusive<NodeType>(0, x, y, z);
                        }
                    }
                }

            }
            rKnotsWRefined[k] = rKnotsWToInsert[j];
            k -= 1;
        }
        
    }

    //########################################################################################
    // OLD STUFF
    //########################################################################################

    void NurbsVolumeRefinementUtilities::DegreeElevationU(
                                        NurbsVolumeGeometryType& rGeometry,
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
            const SizeType nb_cp_w = rGeometry.PointsNumberInDirection(2);

            // Attention: Weights are not yet implemented.
            // Vector weights_old = rGeometry.Weights();
            // if (weights_old.size() != rGeometry.size()) {
            //     weights_old.resize(rGeometry.size());
            //     std::fill(weights_old.begin(), weights_old.end(), 1.0);
            // }
            Vector weights_old; // = rGeometry.Weights();
            if (weights_old.size() != rGeometry.size()) {
                weights_old.resize(rGeometry.size());
                std::fill(weights_old.begin(), weights_old.end(), 1.0);
            }
            //---------------------------

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
            if (rWeightsRefined.size() != nb_cp_u_refined * nb_cp_v * nb_cp_w) {
                rWeightsRefined.resize(nb_cp_u_refined * nb_cp_v * nb_cp_w);
            }
            rWeightsRefined = ZeroVector(nb_cp_u_refined * nb_cp_v * nb_cp_w);
            if (rPointsRefined.size() != nb_cp_u_refined * nb_cp_v * nb_cp_w) {
                rPointsRefined.resize(nb_cp_u_refined * nb_cp_v * nb_cp_w);
            }

            // Initialize additional variables
            Matrix bezier_alphas = ZeroMatrix(degree_u + 1, degree_u_old + 1);
            PointerVector<NodeType> bezier_points;
            if (bezier_points.size() != (degree_u_old + 1) * nb_cp_v * nb_cp_w) {
                bezier_points.resize((degree_u_old + 1) * nb_cp_v * nb_cp_w);
            }
            Vector bezier_points_weights = ZeroVector((degree_u_old + 1) * nb_cp_v * nb_cp_w);
            PointerVector<NodeType> elevate_bezier_points;
            if (elevate_bezier_points.size() != (degree_u + 1) * nb_cp_v * nb_cp_w) {
                elevate_bezier_points.resize((degree_u + 1) * nb_cp_v * nb_cp_w);
            }
            Vector elevate_bezier_points_weights = ZeroVector((degree_u + 1) * nb_cp_v * nb_cp_w);
            PointerVector<NodeType> next_bezier_points;
            if (next_bezier_points.size() != (degree_u_old - 1) * nb_cp_v * nb_cp_w) {
                next_bezier_points.resize((degree_u_old - 1) * nb_cp_v * nb_cp_w);
            }
            Vector next_bezier_points_weights = ZeroVector((degree_u_old - 1) * nb_cp_v * nb_cp_w);

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
                for (IndexType n = 0; n < nb_cp_w; ++n) {
                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_refined, nb_cp_v, nb_cp_w, 0, m, n);

                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u_old, nb_cp_v, nb_cp_w, 0, m, n);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                }
            }

            // Initilaize first Bézier segment
            for (IndexType i = 0; i < degree_u_old + 1; ++i) {
                for (IndexType m = 0; m < nb_cp_v; ++m) {
                    for (IndexType n = 0; n < nb_cp_w; ++n) {
                        IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            (degree_u_old + 1), nb_cp_v, nb_cp_w, i, m, n);
                        IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u_old, nb_cp_v, nb_cp_w, i, m, n);

                        const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                        bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                        bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                    }
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
                                for (IndexType n = 0; n < nb_cp_w; ++n) {
                                    IndexType cp_index_bpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        (degree_u_old + 1), nb_cp_v, nb_cp_w, k, m, n);
                                    IndexType cp_index_bpts_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        (degree_u_old + 1), nb_cp_v, nb_cp_w, k - 1, m, n);

                                    const array_1d<double, 3> cp_coordinates = alphas[k - s] * bezier_points[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points[cp_index_bpts_before];

                                    bezier_points(cp_index_bpts_after) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                    bezier_points_weights[cp_index_bpts_after] = alphas[k - s] * bezier_points_weights[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points_weights[cp_index_bpts_before];
                                }
                            }
                        }

                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            for (IndexType n = 0; n < nb_cp_w; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u_old + 1), nb_cp_v, nb_cp_w, degree_u_old, m, n);
                                IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u_old - 1), nb_cp_v, nb_cp_w, save, m, n);

                                next_bezier_points(cp_index_nextbpts) = bezier_points(cp_index_bpts);
                                next_bezier_points_weights[cp_index_nextbpts] = bezier_points_weights[cp_index_bpts];
                            }
                        }
                    }
                }

                // Degree elevate Bezier
                for (IndexType i = left_bezier; i < degree_u + 1; ++i)
                {
                    for (IndexType m = 0; m < nb_cp_v; ++m) {
                        for (IndexType n = 0; n < nb_cp_w; ++n) {
                            IndexType cp_index_ebpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u + 1), nb_cp_v, nb_cp_w, i, m, n);

                            const array_1d<double, 3> cp_coordinates = ZeroVector(3);

                            elevate_bezier_points(cp_index_ebpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            elevate_bezier_points_weights[cp_index_ebpts] = 0.0;
                        }
                    }

                    SizeType min_degree = std::min(degree_u_old, i);
                    int index = std::max(0, static_cast<int>(i - rDegreeUToElevate));

                    for (IndexType j = index; j < min_degree + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            for (IndexType n = 0; n < nb_cp_w; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u_old + 1), nb_cp_v, nb_cp_w, j, m, n);
                                IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u + 1), nb_cp_v, nb_cp_w, i, m, n);

                                const array_1d<double, 3> cp_coordinates = elevate_bezier_points[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points[cp_index_bpts];

                                elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                elevate_bezier_points_weights[cp_index_etbpts] = elevate_bezier_points_weights[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points_weights[cp_index_bpts];
                            }
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
                                    for (IndexType n = 0; n < nb_cp_w; ++n) {
                                        IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u_refined, nb_cp_v, nb_cp_w, i, m, n);
                                        IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u_refined, nb_cp_v, nb_cp_w, i - 1, m, n);

                                        const array_1d<double, 3> cp_coordinates = alpha * rPointsRefined[cp_index_refined] + (1.0 - alpha) * rPointsRefined[cp_index_refined_before];

                                        rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                        rWeightsRefined[cp_index_refined] = rWeightsRefined[cp_index_refined] * alpha + rWeightsRefined[cp_index_refined_before] * (1 - alpha);
                                    }
                                }
                            }
                            if (j >= left_bezier)
                            {
                                if ((j - t) <= (knot_index - degree_u + r_old))
                                {
                                    double gamma = (knot_u_old_b - rKnotsURefined[j - t - 1]) / (knot_u_old_b - knot_u_old_a);
                                    for (IndexType m = 0; m < nb_cp_v; ++m) {
                                        for (IndexType n = 0; n < nb_cp_w; ++n) {
                                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                (degree_u + 1), nb_cp_v, nb_cp_w, k, m, n);
                                            IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                (degree_u + 1), nb_cp_v, nb_cp_w, k + 1, m, n);

                                            const array_1d<double, 3> cp_coordinates = gamma * elevate_bezier_points[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points[cp_index_etbpts_after];

                                            elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                            elevate_bezier_points_weights[cp_index_etbpts] = gamma * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                        }
                                    }

                                }
                                else
                                {
                                    for (IndexType m = 0; m < nb_cp_v; ++m) {
                                        for (IndexType n = 0; n < nb_cp_w; ++n) {
                                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                (degree_u + 1), nb_cp_v, nb_cp_w, k, m, n);
                                            IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                (degree_u + 1), nb_cp_v, nb_cp_w, k + 1, m, n);

                                            const array_1d<double, 3> cp_coordinates = beta * elevate_bezier_points[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points[cp_index_etbpts_after];

                                            elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                            elevate_bezier_points_weights[cp_index_etbpts] = beta * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                        }
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
                        for (IndexType n = 0; n < nb_cp_w; ++n) {
                            IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u_refined, nb_cp_v, nb_cp_w, control_point_index, m, n);
                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                (degree_u + 1), nb_cp_v, nb_cp_w, j, m, n);

                            rPointsRefined(cp_index) = elevate_bezier_points(cp_index_etbpts);
                            rWeightsRefined[cp_index] = elevate_bezier_points_weights[cp_index_etbpts];
                        }
                    }
                    control_point_index = control_point_index + 1;
                }

                // Setup for the next pass through loop
                if (b < nb_knots_u_old)
                {
                    for (IndexType j = 0; j < static_cast<IndexType>(r); ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            for (IndexType n = 0; n < nb_cp_w; ++n) {
                                IndexType cp_index_btps = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u_old + 1), nb_cp_v, nb_cp_w, j, m, n);
                                IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u_old - 1), nb_cp_v, nb_cp_w, j, m, n);

                                bezier_points(cp_index_btps) = next_bezier_points(cp_index_nextbpts);
                                bezier_points_weights[cp_index_btps] = next_bezier_points_weights[cp_index_nextbpts];
                            }
                        }
                    }
                    for (IndexType j = r; j < degree_u_old + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_v; ++m) {
                            for (IndexType n = 0; n < nb_cp_w; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    (degree_u_old + 1), nb_cp_v, nb_cp_w, j, m, n);
                                IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    nb_cp_u_old, nb_cp_v, nb_cp_w, (b - degree_u_old + j + 1), m, n);

                                const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                                bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                            }
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
            
            for (IndexType i = 0; i < nb_cp_u_refined * nb_cp_v * nb_cp_w; ++i) {
                rPointsRefined(i)->Coordinates() = rPointsRefined(i)->Coordinates() / rWeightsRefined[i];
            }
        }
        else {
            KRATOS_WARNING("::NurbsVolumeRefinementUtilities::DegreeElevationU") << "No elevation applied as degree is zero. "
                << std::endl;
        }
    }

    void NurbsVolumeRefinementUtilities::DegreeElevationV(
                                        NurbsVolumeGeometryType& rGeometry,
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
            const SizeType nb_cp_w = rGeometry.PointsNumberInDirection(2);

            // Attention: Weights are not yet implemented.
            // Vector weights_old = rGeometry.Weights();
            // if (weights_old.size() != rGeometry.size()) {
            //     weights_old.resize(rGeometry.size());
            //     std::fill(weights_old.begin(), weights_old.end(), 1.0);
            // }
            Vector weights_old; // = rGeometry.Weights();
            if (weights_old.size() != rGeometry.size()) {
                weights_old.resize(rGeometry.size());
                std::fill(weights_old.begin(), weights_old.end(), 1.0);
            }
            //---------------------------

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
            if (rWeightsRefined.size() != nb_cp_v_refined * nb_cp_u * nb_cp_w) {
                rWeightsRefined.resize(nb_cp_v_refined * nb_cp_u * nb_cp_w);
            }
            rWeightsRefined = ZeroVector(nb_cp_v_refined * nb_cp_u * nb_cp_w);
            if (rPointsRefined.size() != nb_cp_v_refined * nb_cp_u * nb_cp_w) {
                rPointsRefined.resize(nb_cp_v_refined * nb_cp_u * nb_cp_w);
            }

            // Initialize additional variables
            Matrix bezier_alphas = ZeroMatrix(degree_v + 1, degree_v_old + 1);
            PointerVector<NodeType> bezier_points;
            if (bezier_points.size() != (degree_v_old + 1) * nb_cp_u * nb_cp_w) {
                bezier_points.resize((degree_v_old + 1) * nb_cp_u * nb_cp_w);
            }
            Vector bezier_points_weights = ZeroVector((degree_v_old + 1) * nb_cp_u * nb_cp_w);
            PointerVector<NodeType> elevate_bezier_points;
            if (elevate_bezier_points.size() != (degree_v + 1) * nb_cp_u * nb_cp_w) {
                elevate_bezier_points.resize((degree_v + 1) * nb_cp_u * nb_cp_w);
            }
            Vector elevate_bezier_points_weights = ZeroVector((degree_v + 1) * nb_cp_u * nb_cp_w);
            PointerVector<NodeType> next_bezier_points;
            if (next_bezier_points.size() != (degree_v_old - 1) * nb_cp_u * nb_cp_w) {
                next_bezier_points.resize((degree_v_old - 1) * nb_cp_u * nb_cp_w);
            }
            Vector next_bezier_points_weights = ZeroVector((degree_v_old - 1) * nb_cp_u * nb_cp_w);

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
                for (IndexType n = 0; n < nb_cp_w; ++n) {
                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_refined, nb_cp_w, m, 0, n);

                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v_old, nb_cp_w, m, 0, n);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                }
            }

            // Initilaize first Bézier segment
            for (IndexType i = 0; i < degree_v_old + 1; ++i) {
                for (IndexType m = 0; m < nb_cp_u; ++m) {
                    for (IndexType n = 0; n < nb_cp_w; ++n) {
                        IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u,  (degree_v_old + 1),  nb_cp_w, m, i, n);
                        IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u, nb_cp_v_old, nb_cp_w, m, i, n);

                        const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                        bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                        bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                    }
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
                                for (IndexType n = 0; n < nb_cp_w; ++n) {
                                    IndexType cp_index_bpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        nb_cp_u, (degree_v_old + 1), nb_cp_w, m, k, n);
                                    IndexType cp_index_bpts_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        nb_cp_u, (degree_v_old + 1), nb_cp_w, m, k-1, n);

                                    const array_1d<double, 3> cp_coordinates = alphas[k - s] * bezier_points[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points[cp_index_bpts_before];

                                    bezier_points(cp_index_bpts_after) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                    bezier_points_weights[cp_index_bpts_after] = alphas[k - s] * bezier_points_weights[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points_weights[cp_index_bpts_before];
                                }
                            }
                        }

                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            for (IndexType n = 0; n < nb_cp_w; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, (degree_v_old + 1), nb_cp_w, m, degree_v_old, n);
                                IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, (degree_v_old - 1), nb_cp_w, m, save, n);

                                next_bezier_points(cp_index_nextbpts) = bezier_points(cp_index_bpts);
                                next_bezier_points_weights[cp_index_nextbpts] = bezier_points_weights[cp_index_bpts];
                            }
                        }
                    }
                }

                // Degree elevate Bezier
                for (IndexType i = left_bezier; i < degree_v + 1; ++i)
                {
                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                        for (IndexType n = 0; n < nb_cp_w; ++n) {
                            IndexType cp_index_ebpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v + 1), nb_cp_w, m, i, n);

                            const array_1d<double, 3> cp_coordinates = ZeroVector(3);

                            elevate_bezier_points(cp_index_ebpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            elevate_bezier_points_weights[cp_index_ebpts] = 0.0;
                        }
                    }

                    SizeType min_degree = std::min(degree_v_old, i);
                    int index = std::max(0, static_cast<int>(i - rDegreeVToElevate));

                    for (IndexType j = index; j < min_degree + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            for (IndexType n = 0; n < nb_cp_w; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    nb_cp_u,(degree_v_old + 1), nb_cp_w, m, j, n);
                                IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    nb_cp_u, (degree_v + 1), nb_cp_w, m, i, n);

                                const array_1d<double, 3> cp_coordinates = elevate_bezier_points[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points[cp_index_bpts];

                                elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                elevate_bezier_points_weights[cp_index_etbpts] = elevate_bezier_points_weights[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points_weights[cp_index_bpts];
                            }
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
                                    for (IndexType n = 0; n < nb_cp_w; ++n) {
                                        IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u, nb_cp_v_refined, nb_cp_w, m, i, n);
                                        IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u, nb_cp_v_refined, nb_cp_w, m, i-1, n);

                                        const array_1d<double, 3> cp_coordinates = alpha * rPointsRefined[cp_index_refined] + (1.0 - alpha) * rPointsRefined[cp_index_refined_before];

                                        rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                        rWeightsRefined[cp_index_refined] = rWeightsRefined[cp_index_refined] * alpha + rWeightsRefined[cp_index_refined_before] * (1 - alpha);
                                    }
                                }
                            }
                            if (j >= left_bezier)
                            {
                                if ((j - t) <= (knot_index - degree_v + r_old))
                                {
                                    double gamma = (knot_v_old_b - rKnotsVRefined[j - t - 1]) / (knot_v_old_b - knot_v_old_a);
                                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                                        for (IndexType n = 0; n < nb_cp_w; ++n) {
                                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                 nb_cp_u, (degree_v + 1), nb_cp_w, m, k, n);
                                            IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                 nb_cp_u, (degree_v + 1),nb_cp_w, m, k+1, n);

                                            const array_1d<double, 3> cp_coordinates = gamma * elevate_bezier_points[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points[cp_index_etbpts_after];

                                            elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                            elevate_bezier_points_weights[cp_index_etbpts] = gamma * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                        }
                                    }

                                }
                                else
                                {
                                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                                        for (IndexType n = 0; n < nb_cp_w; ++n) {
                                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                nb_cp_u, (degree_v + 1), nb_cp_w, m, k, n);
                                            IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                nb_cp_u, (degree_v + 1), nb_cp_w, m, k+1, n);

                                            const array_1d<double, 3> cp_coordinates = beta * elevate_bezier_points[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points[cp_index_etbpts_after];

                                            elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                            elevate_bezier_points_weights[cp_index_etbpts] = beta * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                        }
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
                        for (IndexType n = 0; n < nb_cp_w; ++n) {
                            IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v_refined, nb_cp_w, m, control_point_index, n);
                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, (degree_v + 1), nb_cp_w, m, j, n);

                            rPointsRefined(cp_index) = elevate_bezier_points(cp_index_etbpts);
                            rWeightsRefined[cp_index] = elevate_bezier_points_weights[cp_index_etbpts];
                        }
                    }
                    control_point_index = control_point_index + 1;
                }

                // Setup for the next pass through loop
                if (b < nb_knots_v_old)
                {
                    for (IndexType j = 0; j < static_cast<IndexType>(r); ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            for (IndexType n = 0; n < nb_cp_w; ++n) {
                                IndexType cp_index_btps = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, (degree_v_old + 1), nb_cp_w, m, j, n);
                                IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, (degree_v_old - 1), nb_cp_w, m, j, n);

                                bezier_points(cp_index_btps) = next_bezier_points(cp_index_nextbpts);
                                bezier_points_weights[cp_index_btps] = next_bezier_points_weights[cp_index_nextbpts];
                            }
                        }
                    }
                    for (IndexType j = r; j < degree_v_old + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            for (IndexType n = 0; n < nb_cp_w; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, (degree_v_old + 1), nb_cp_w, m, j, n);
                                IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, nb_cp_v_old, nb_cp_w, m, (b - degree_v_old + j + 1), n);

                                const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                                bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                            }
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
            
            for (IndexType i = 0; i < nb_cp_v_refined * nb_cp_u * nb_cp_w; ++i) {
                rPointsRefined(i)->Coordinates() = rPointsRefined(i)->Coordinates() / rWeightsRefined[i];
            }
        }
        else {
            KRATOS_WARNING("::NurbsVolumeRefinementUtilities::DegreeElevationU") << "No elevation applied as degree is zero. "
                << std::endl;
        }
    }


    void NurbsVolumeRefinementUtilities::DegreeElevationW(
                                        NurbsVolumeGeometryType& rGeometry,
                                        SizeType& rDegreeWToElevate,
                                        PointerVector<NodeType>& rPointsRefined,
                                        Vector& rKnotsWRefined,
                                        Vector& rWeightsRefined)
    {
        const Vector& knots_w_old = rGeometry.KnotsW();

        if (rDegreeWToElevate > 0) {
            const SizeType degree_w_old = rGeometry.PolynomialDegree(2);
            const SizeType degree_w = degree_w_old + rDegreeWToElevate;

            const SizeType nb_cp_u = rGeometry.PointsNumberInDirection(0);
            const SizeType nb_cp_v = rGeometry.PointsNumberInDirection(1);
            const SizeType nb_cp_w_old = rGeometry.PointsNumberInDirection(2);

            // Attention: Weights are not yet implemented.
            // Vector weights_old = rGeometry.Weights();
            // if (weights_old.size() != rGeometry.size()) {
            //     weights_old.resize(rGeometry.size());
            //     std::fill(weights_old.begin(), weights_old.end(), 1.0);
            // }
            Vector weights_old; // = rGeometry.Weights();
            if (weights_old.size() != rGeometry.size()) {
                weights_old.resize(rGeometry.size());
                std::fill(weights_old.begin(), weights_old.end(), 1.0);
            }
            //---------------------------

            const SizeType nb_knots_w_old = knots_w_old.size();
            SizeType number_of_non_zero_spans_w = 0;
            for (IndexType i = 0; i < nb_knots_w_old - 1; i++)
            {
                if (knots_w_old[i] != knots_w_old[i + 1])
                {
                    number_of_non_zero_spans_w = number_of_non_zero_spans_w + 1; 
                } 
            }

            const SizeType nb_knots_w_refined = nb_knots_w_old + rDegreeWToElevate * (number_of_non_zero_spans_w + 1);
            const SizeType nb_cp_w_refined = nb_knots_w_refined - degree_w + 1;

            // Resize reference variables
            if (rKnotsWRefined.size() != nb_knots_w_refined) {
                rKnotsWRefined.resize(nb_knots_w_refined);
            }
            rKnotsWRefined = ZeroVector(nb_knots_w_refined);
            if (rWeightsRefined.size() != nb_cp_w_refined * nb_cp_u * nb_cp_v) {
                rWeightsRefined.resize(nb_cp_w_refined * nb_cp_u * nb_cp_v);
            }
            rWeightsRefined = ZeroVector(nb_cp_w_refined * nb_cp_u * nb_cp_v);
            if (rPointsRefined.size() != nb_cp_w_refined * nb_cp_u * nb_cp_v) {
                rPointsRefined.resize(nb_cp_w_refined * nb_cp_u * nb_cp_v);
            }

            // Initialize additional variables
            Matrix bezier_alphas = ZeroMatrix(degree_w + 1, degree_w_old + 1);
            PointerVector<NodeType> bezier_points;
            if (bezier_points.size() != (degree_w_old + 1) * nb_cp_u * nb_cp_v) {
                bezier_points.resize((degree_w_old + 1) * nb_cp_u * nb_cp_v);
            }
            Vector bezier_points_weights = ZeroVector((degree_w_old + 1) * nb_cp_u * nb_cp_v);
            PointerVector<NodeType> elevate_bezier_points;
            if (elevate_bezier_points.size() != (degree_w + 1) * nb_cp_u * nb_cp_v) {
                elevate_bezier_points.resize((degree_w + 1) * nb_cp_u * nb_cp_v);
            }
            Vector elevate_bezier_points_weights = ZeroVector((degree_w + 1) * nb_cp_u * nb_cp_v);
            PointerVector<NodeType> next_bezier_points;
            if (next_bezier_points.size() != (degree_w_old - 1) * nb_cp_u * nb_cp_v) {
                next_bezier_points.resize((degree_w_old - 1) * nb_cp_u * nb_cp_v);
            }
            Vector next_bezier_points_weights = ZeroVector((degree_w_old - 1) * nb_cp_u * nb_cp_v);

            // Compute Bézier degree elevation coefficients
            bezier_alphas(0, 0) = 1.0;
            bezier_alphas(degree_w, degree_w_old) = 1.0;

            for(IndexType i = 1; i < (degree_w / 2) + 1; ++i)
            {
                double inv = 1.0 / NurbsUtilities::GetBinomCoefficient(degree_w, i);
                SizeType min_degree = std::min(degree_w_old, i);
                int index = std::max(0, static_cast<int>(i - rDegreeWToElevate));

                for(IndexType j = index; j < min_degree + 1; ++j)
                {
                    bezier_alphas(i, j) = inv * NurbsUtilities::GetBinomCoefficient(degree_w_old, j) * NurbsUtilities::GetBinomCoefficient(rDegreeWToElevate, i - j);
                    bezier_alphas(degree_w - i, degree_w_old - j) = bezier_alphas(i, j);
                }
            }

            // Initialize coefficients
            SizeType a = degree_w_old - 1;
            SizeType b = degree_w_old;
            int r = -1;
            IndexType knot_index = degree_w + 1;
            IndexType control_point_index = 1;
            double knot_w_old_a = knots_w_old[a];

            // Create new knot span vector
            for (IndexType i = 0; i < degree_w; ++i) {
                rKnotsWRefined[i] = knot_w_old_a;
            }

            // Create new node vector
            for (IndexType m = 0; m < nb_cp_u; ++m) {
                for (IndexType n = 0; n < nb_cp_v; ++n) {
                    IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v, nb_cp_w_refined,  m, n, 0);

                    IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        nb_cp_u, nb_cp_v, nb_cp_w_old, m, n, 0);

                    const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                    rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                    rWeightsRefined[cp_index_refined] = weights_old[cp_index_old];
                }
            }

            // Initilaize first Bézier segment
            for (IndexType i = 0; i < degree_w_old + 1; ++i) {
                for (IndexType m = 0; m < nb_cp_u; ++m) {
                    for (IndexType n = 0; n < nb_cp_v; ++n) {
                        IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u,  nb_cp_v, (degree_w_old + 1), m, n, i);
                        IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            nb_cp_u, nb_cp_v, nb_cp_w_old, m, n, i);

                        const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                        bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                        bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                    }
                }
            }

            while (b < nb_knots_w_old)
            {
                IndexType i = b;
                while (b < (nb_knots_w_old - 1) && knots_w_old[b] == knots_w_old[b + 1])
                {
                    b++;
                }
                if (b + 1 == nb_knots_w_old)
                {
                    b++;
                } 
                IndexType mult = b - i + 1;
                double knot_w_old_b = knots_w_old[i];
                int r_old = r;
                r = degree_w_old - mult;

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
                    right_bezier = degree_w - ((r + 1) / 2);  
                }
                else
                {
                    right_bezier = degree_w;
                }
                
                // Insert knot to get Bezier segment
                if (r > 0)
                {
                    Vector alphas = ZeroVector(degree_w_old - 1);
                    for (IndexType k = degree_w_old; k > mult; --k)
                    {
                        alphas[k - mult - 1] = (knot_w_old_b - knot_w_old_a)/(knots_w_old[a + k] - knot_w_old_a); 
                    }

                    for (IndexType j = 1; j < static_cast<IndexType>(r + 1); ++j) 
                    {
                        int save = r - j;
                        IndexType s = mult + j;

                        for (IndexType k = degree_w_old; k >= s; --k)
                        {
                            for (IndexType m = 0; m < nb_cp_u; ++m) {
                                for (IndexType n = 0; n < nb_cp_v; ++n) {
                                    IndexType cp_index_bpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        nb_cp_u, nb_cp_v, (degree_w_old + 1),m, n, k);
                                    IndexType cp_index_bpts_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                        nb_cp_u, nb_cp_v, (degree_w_old + 1), m, n, k-1);

                                    const array_1d<double, 3> cp_coordinates = alphas[k - s] * bezier_points[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points[cp_index_bpts_before];

                                    bezier_points(cp_index_bpts_after) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                    bezier_points_weights[cp_index_bpts_after] = alphas[k - s] * bezier_points_weights[cp_index_bpts_after] + (1 - alphas[k - s]) * bezier_points_weights[cp_index_bpts_before];
                                }
                            }
                        }

                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            for (IndexType n = 0; n < nb_cp_v; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, nb_cp_v, (degree_w_old + 1), m, n, degree_w_old);
                                IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, nb_cp_v, (degree_w_old - 1), m, n, save);

                                next_bezier_points(cp_index_nextbpts) = bezier_points(cp_index_bpts);
                                next_bezier_points_weights[cp_index_nextbpts] = bezier_points_weights[cp_index_bpts];
                            }
                        }
                    }
                }

                // Degree elevate Bezier
                for (IndexType i = left_bezier; i < degree_w + 1; ++i)
                {
                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                        for (IndexType n = 0; n < nb_cp_v; ++n) {
                            IndexType cp_index_ebpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v, (degree_w + 1), m, n, i);

                            const array_1d<double, 3> cp_coordinates = ZeroVector(3);

                            elevate_bezier_points(cp_index_ebpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                            elevate_bezier_points_weights[cp_index_ebpts] = 0.0;
                        }
                    }

                    SizeType min_degree = std::min(degree_w_old, i);
                    int index = std::max(0, static_cast<int>(i - rDegreeWToElevate));

                    for (IndexType j = index; j < min_degree + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            for (IndexType n = 0; n < nb_cp_v; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    nb_cp_u, nb_cp_v, (degree_w_old + 1), m, n, j);
                                IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                    nb_cp_u, nb_cp_v, (degree_w + 1), m, n, i);

                                const array_1d<double, 3> cp_coordinates = elevate_bezier_points[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points[cp_index_bpts];

                                elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                elevate_bezier_points_weights[cp_index_etbpts] = elevate_bezier_points_weights[cp_index_etbpts] + bezier_alphas(i, j) * bezier_points_weights[cp_index_bpts];
                            }
                        }
                    }
                }

                // Knot removal 
                if (r_old > 1)
                {
                    IndexType first = knot_index - 2;
                    IndexType last = knot_index;

                    double beta = (knot_w_old_b - rKnotsWRefined[knot_index - 2]) / (knot_w_old_b - knot_w_old_a);

                    for (IndexType t = 1; t < static_cast<IndexType>(r_old); ++t)
                    {
                        IndexType i = first;
                        IndexType j = last;
                        IndexType k = j - knot_index + 1;

                        while (j - i > t)
                        {
                            if (i < control_point_index)
                            {
                                double alpha = (knot_w_old_b - rKnotsWRefined[i - 1]) / (knot_w_old_a - rKnotsWRefined[i - 1]);
                                for (IndexType m = 0; m < nb_cp_u; ++m) {
                                    for (IndexType n = 0; n < nb_cp_v; ++n) {
                                        IndexType cp_index_refined = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u, nb_cp_v, nb_cp_w_refined, m, n, i);
                                        IndexType cp_index_refined_before = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                            nb_cp_u, nb_cp_v, nb_cp_w_refined, m, n, i-1);

                                        const array_1d<double, 3> cp_coordinates = alpha * rPointsRefined[cp_index_refined] + (1.0 - alpha) * rPointsRefined[cp_index_refined_before];

                                        rPointsRefined(cp_index_refined) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                        rWeightsRefined[cp_index_refined] = rWeightsRefined[cp_index_refined] * alpha + rWeightsRefined[cp_index_refined_before] * (1 - alpha);
                                    }
                                }
                            }
                            if (j >= left_bezier)
                            {
                                if ((j - t) <= (knot_index - degree_w + r_old))
                                {
                                    double gamma = (knot_w_old_b - rKnotsWRefined[j - t - 1]) / (knot_w_old_b - knot_w_old_a);
                                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                                        for (IndexType n = 0; n < nb_cp_v; ++n) {
                                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                 nb_cp_u, nb_cp_v, (degree_w + 1), m, n, k);
                                            IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                 nb_cp_u, nb_cp_v, (degree_w + 1), m, n, k+1);

                                            const array_1d<double, 3> cp_coordinates = gamma * elevate_bezier_points[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points[cp_index_etbpts_after];

                                            elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                            elevate_bezier_points_weights[cp_index_etbpts] = gamma * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - gamma) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                        }
                                    }

                                }
                                else
                                {
                                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                                        for (IndexType n = 0; n < nb_cp_v; ++n) {
                                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                nb_cp_u,  nb_cp_v, (degree_w + 1), m, n, k);
                                            IndexType cp_index_etbpts_after = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                                nb_cp_u, nb_cp_v, (degree_w + 1), m, n, k+1);

                                            const array_1d<double, 3> cp_coordinates = beta * elevate_bezier_points[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points[cp_index_etbpts_after];

                                            elevate_bezier_points(cp_index_etbpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                            elevate_bezier_points_weights[cp_index_etbpts] = beta * elevate_bezier_points_weights[cp_index_etbpts] + (1.0 - beta) * elevate_bezier_points_weights[cp_index_etbpts_after];
                                        }
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
                if ((a + 1) != degree_w_old) 
                {
                    for (IndexType i = 0; i < (degree_w - r_old); ++i)
                    {
                        rKnotsWRefined[knot_index - 1] = knot_w_old_a;
                        knot_index = knot_index + 1;
                    }
                }

                // Create next node vector
                for (IndexType j = left_bezier; j < right_bezier + 1; ++j)
                {
                    for (IndexType m = 0; m < nb_cp_u; ++m) {
                        for (IndexType n = 0; n < nb_cp_v; ++n) {
                            IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v, nb_cp_w_refined, m, n, control_point_index);
                            IndexType cp_index_etbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                nb_cp_u, nb_cp_v, (degree_w + 1), m, n, j);

                            rPointsRefined(cp_index) = elevate_bezier_points(cp_index_etbpts);
                            rWeightsRefined[cp_index] = elevate_bezier_points_weights[cp_index_etbpts];
                        }
                    }
                    control_point_index = control_point_index + 1;
                }

                // Setup for the next pass through loop
                if (b < nb_knots_w_old)
                {
                    for (IndexType j = 0; j < static_cast<IndexType>(r); ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            for (IndexType n = 0; n < nb_cp_v; ++n) {
                                IndexType cp_index_btps = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, nb_cp_v, (degree_w_old + 1), m, n, j);
                                IndexType cp_index_nextbpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, nb_cp_v, (degree_w_old - 1), m, n, j);

                                bezier_points(cp_index_btps) = next_bezier_points(cp_index_nextbpts);
                                bezier_points_weights[cp_index_btps] = next_bezier_points_weights[cp_index_nextbpts];
                            }
                        }
                    }
                    for (IndexType j = r; j < degree_w_old + 1; ++j)
                    {
                        for (IndexType m = 0; m < nb_cp_u; ++m) {
                            for (IndexType n = 0; n < nb_cp_v; ++n) {
                                IndexType cp_index_bpts = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, nb_cp_v, (degree_w_old + 1), m, n, j);
                                IndexType cp_index_old = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                     nb_cp_u, nb_cp_v, nb_cp_w_old, m, n, (b - degree_w_old + j + 1));

                                const array_1d<double, 3> cp_coordinates = rGeometry[cp_index_old] * weights_old[cp_index_old];

                                bezier_points(cp_index_bpts) = Kratos::make_intrusive<NodeType>(0, cp_coordinates);
                                bezier_points_weights[cp_index_bpts] = weights_old[cp_index_old];
                            }
                        }
                    }

                    a = b;
                    b = b + 1;
                    knot_w_old_a = knot_w_old_b;
                }
                else
                {
                    for (IndexType i = 0; i < degree_w; ++i)
                    {
                        rKnotsWRefined[knot_index - 1 + i] = knot_w_old_b;
                    }
                }
            }
            
            for (IndexType i = 0; i < nb_cp_w_refined * nb_cp_u * nb_cp_v; ++i) {
                rPointsRefined(i)->Coordinates() = rPointsRefined(i)->Coordinates() / rWeightsRefined[i];
            }
        }
        else {
            KRATOS_WARNING("::NurbsVolumeRefinementUtilities::DegreeElevationU") << "No elevation applied as degree is zero. "
                << std::endl;
        }
    }





} // End namespace kratos