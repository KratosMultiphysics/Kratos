//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_NURBS_KNOT_INSERTION_H_INCLUDED )
#define  KRATOS_NURBS_KNOT_INSERTION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/nurbs_volume_geometry.h"
#include "nurbs_utilities.h"
#include "includes/node.h"

namespace Kratos {

template<class TContainerPointType>
class KnotInsertion
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node<3> NodeType;

    typedef NurbsVolumeGeometry<TContainerPointType> NurbsVolumeGeometryType;

    ///@}
    ///@name Operations
    ///@{
    static void KnotRefinement(NurbsVolumeGeometryType& rGeometry, std::vector<double> u ){

        std::sort(u.begin(),u.end());
        // Get current order
        SizeType polynomial_degree_u = rGeometry.PolynomialDegreeU();
        SizeType polynomial_degree_v = rGeometry.PolynomialDegreeV();
        SizeType polynomial_degree_w = rGeometry.PolynomialDegreeW();
        // Get current knot information
        const Kratos::Vector& old_knots_u = rGeometry.KnotsU();
        const Kratos::Vector& old_knots_v = rGeometry.KnotsV();
        const Kratos::Vector& old_knots_w = rGeometry.KnotsW();

        SizeType old_num_of_knots_u = rGeometry.NumberOfKnotsU();
        // Get current cp's information
        SizeType old_num_of_cp_u = rGeometry.NumberOfControlPointsU();
        SizeType old_num_of_cp_v = rGeometry.NumberOfControlPointsV();
        SizeType old_num_of_cp_w = rGeometry.NumberOfControlPointsW();

        SizeType a = NurbsUtilities::GetLowerSpan(polynomial_degree_u, old_knots_u, u[0]);
        std::cout << "a: " << a << std::endl;
        // Check the +1 again
        SizeType b = NurbsUtilities::GetLowerSpan(polynomial_degree_u, old_knots_u, u[u.size()-1]) + 1;
        std::cout << "b: " << b << std::endl;

        SizeType r = u.size();
        // Initialize new containers
        SizeType new_num_of_knots_u = old_num_of_knots_u + r;
        Kratos::Vector new_knots_u(new_num_of_knots_u);

        SizeType new_num_of_cp_u = old_num_of_cp_u + r;
        PointerVector<NodeType> new_points(new_num_of_cp_u*old_num_of_cp_v*old_num_of_cp_w);

        r = r - 1;

        // Create new ordered knot vector.
        for( IndexType i = 0; i <= a; i++)
            new_knots_u[i] = old_knots_u[i];
        for( IndexType i = b+polynomial_degree_u; i < old_num_of_knots_u; ++i)
            new_knots_u[i+r+1] = old_knots_u[i];

        //Copy unaltered cp's.
        SizeType id = 1;
        for( IndexType column=0; column < old_num_of_cp_v; ++column){
            for( IndexType depth=0; depth < old_num_of_cp_w; ++depth){
                //Copy unaltered control points.
                //Note: The +1 accounts for the fact that we have only p-multiplicity at the knot ends.
                for( IndexType i = 0; i <= a - polynomial_degree_u + 1; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    new_points(cp_index_left) = NodeType::Pointer(new NodeType(id, rGeometry[cp_index_right][0], rGeometry[cp_index_right][1], rGeometry[cp_index_right][2]));
                    id += 1;
                }
                // The +1 accounts for the fact that we have only p-multiplicity at the knot ends.
                for( IndexType i = b; i < old_num_of_cp_u; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i+r+1, column, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    new_points(cp_index_left) = NodeType::Pointer(new NodeType(id, rGeometry[cp_index_right][0], rGeometry[cp_index_right][1], rGeometry[cp_index_right][2]));
                    id += 1;
                }
            }
        }

        int p = polynomial_degree_u;
        int i = b + p - 1;
        int k = b + p + r;
        for( int j = r; j >= 0; j--){ // Could be whole loop as well
            while( u[j] <= old_knots_u[i] && i > a){
                for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                    for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                        IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, k-p-1 +1, column, depth);
                        IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                            old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i - p -1 + 1, column, depth);
                        new_points(cp_index_left) = NodeType::Pointer(new NodeType(id, rGeometry[cp_index_right][0], rGeometry[cp_index_right][1], rGeometry[cp_index_right][2]));
                        id++;
                    }
                }
                new_knots_u[k] = old_knots_u[i];
                k--;
                i--;
                std::cout << "llop very uno" << std::endl;
            }
            for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, k-p-1 +1, column, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, k-p+1, column, depth);
                    new_points(cp_index_left) = NodeType::Pointer(new NodeType(id, new_points[cp_index_right][0], new_points[cp_index_right][1], new_points[cp_index_right][2]));
                    id++;
                }
            }
            for( IndexType l=1; l <= p; ++l){
                IndexType index = k - p + l;
                double alpha = new_knots_u[k+l] - u[j];
                if( std::abs(alpha) < 1e-10){
                    for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                        for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                            IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index-1+1, column, depth);
                            IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index+1, column, depth);
                            new_points(cp_index_left) = NodeType::Pointer(new NodeType(id, new_points[cp_index_right][0], new_points[cp_index_right][1], new_points[cp_index_right][2]));
                            id++;
                        }
                    }
                    std::cout << "llop uno" << std::endl;
                }
                else {
                    alpha = alpha / ( new_knots_u[k+l] - old_knots_u[i-p+l]);
                    for( IndexType column=0; column < old_num_of_cp_v; ++column) {
                        for( IndexType depth=0; depth < old_num_of_cp_w; ++depth) {
                            IndexType cp_index_minus = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index-1+1, column, depth);
                            IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                                new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index+1, column, depth);

                            double x = alpha * new_points[cp_index_minus][0] + (1.0-alpha) * new_points[cp_index][0];
                            double y = alpha * new_points[cp_index_minus][1] + (1.0-alpha) * new_points[cp_index][1];
                            double z = alpha * new_points[cp_index_minus][2] + (1.0-alpha) * new_points[cp_index][2];

                            new_points(cp_index_minus) = NodeType::Pointer(new NodeType(id, x, y, z));
                            id++;
                        }
                    }
                    std::cout << "llop doce" << std::endl;
                }

            }
            new_knots_u[k] = u[j];
            k -= 1;
        }

        rGeometry = NurbsVolumeGeometry<PointerVector<NodeType>>(
            new_points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
            new_knots_u, old_knots_v, old_knots_w);

    }

    static void Insert_knot_u(NurbsVolumeGeometryType& rGeometry, double u, SizeType r){

        // Get current order
        SizeType polynomial_degree_u = rGeometry.PolynomialDegreeU();
        SizeType polynomial_degree_v = rGeometry.PolynomialDegreeV();
        SizeType polynomial_degree_w = rGeometry.PolynomialDegreeW();
        // Get current knot information
        const Kratos::Vector& old_knots_u = rGeometry.KnotsU();
        const Kratos::Vector& old_knots_v = rGeometry.KnotsV();
        const Kratos::Vector& old_knots_w = rGeometry.KnotsW();

        SizeType old_num_of_knots_u = rGeometry.NumberOfKnotsU();
        // Get current cp's information
        SizeType old_num_of_cp_u = rGeometry.NumberOfControlPointsU();
        SizeType old_num_of_cp_v = rGeometry.NumberOfControlPointsV();
        SizeType old_num_of_cp_w = rGeometry.NumberOfControlPointsW();
        // Initial multiplicity
        // Todo: Add the capability to insert knots that already exist.
        SizeType s = 0;

        // Initialize new containers
        SizeType new_num_of_knots_u = old_num_of_knots_u + r;
        Kratos::Vector new_knots_u(new_num_of_knots_u);

        SizeType new_num_of_cp_u = old_num_of_cp_u + r;
        PointerVector<NodeType> new_points(new_num_of_cp_u*old_num_of_cp_v*old_num_of_cp_w);

        // Get Index of the lower knotspan [pole_u ... u ... span_end]
        SizeType pole_u = NurbsUtilities::GetLowerSpan(polynomial_degree_u, old_knots_u, u);

        // Create new ordered knot vector. Insert u r-times after pole_u.
        for( IndexType i = 0; i <= pole_u; i++)
            new_knots_u[i] = old_knots_u[i];
        for( IndexType i = 1; i <= r; ++i)
            new_knots_u[pole_u+i] = u;
        for( IndexType i = pole_u+1; i < old_num_of_knots_u; ++i)
            new_knots_u[i+r] = old_knots_u[i];

        // Store the alphas
        SizeType L;
        Matrix alpha;
        alpha.resize(polynomial_degree_u - 1 -s, r);
        for( IndexType j = 1; j <= r; j++){
            L = pole_u - polynomial_degree_u + j;
            for( IndexType i = 0; i <= polynomial_degree_u - j -s; i++){
                alpha(i,j-1) = (u - old_knots_u[L+i])/(old_knots_u[i+pole_u+1] - old_knots_u[L+i]);
            }
        }
        std::cout << "Hallo" << std::endl;
        for( IndexType column=0; column < old_num_of_cp_v; ++column){
            for( IndexType depth=0; depth < old_num_of_cp_w; ++depth){
                //Copy unaltered control points.
                //Note: The +1 accounts for the fact that we have only p-multiplicity at the knot ends.
                std::cout << "Hallo22" << std::endl;
                SizeType id = 1;
                for( IndexType i = 0; i <= pole_u - polynomial_degree_u + 1; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    new_points(cp_index_left) = NodeType::Pointer(new NodeType(id, rGeometry[cp_index_right][0], rGeometry[cp_index_right][1], rGeometry[cp_index_right][2]));
                    id += 1;
                }
                // The +1 accounts for the fact that we have only p-multiplicity at the knot ends.
                for( IndexType i = pole_u + 1; i < old_num_of_cp_u; ++i){
                    IndexType cp_index_left = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i+r, column, depth);
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i, column, depth);
                    new_points(cp_index_left) = NodeType::Pointer(new NodeType(id, rGeometry[cp_index_right][0], rGeometry[cp_index_right][1], rGeometry[cp_index_right][2]));
                    id += 1;
                }

                // Temporary store cp's which will be replaces
                // P - s - 1 will be replaced
                // Starting at k - p + 1 (k:pole_u)
                // Note: The +1 accounts for the fact that we have only p-multiplicity at the knot ends.
                PointerVector<NodeType> tmp_points(old_num_of_cp_u*old_num_of_cp_v*old_num_of_cp_w);
                for( IndexType i = 0; i <= polynomial_degree_u - s; ++i){
                    int index = pole_u - polynomial_degree_u + i + 1;
                    IndexType cp_index_right = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        old_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, index, column, depth);
                    tmp_points(i) = rGeometry(cp_index_right);
                }

                // Find new cp's
                IndexType L = 100;
                for( IndexType j = 1; j <= r; ++j){ // Insert knot r-times
                    L = pole_u - polynomial_degree_u + j;
                    for( IndexType i = 0; i <= polynomial_degree_u - j - s; i++){
                        tmp_points[i][0] = alpha(i,j-1) * tmp_points[i+1][0] + (1.0 - alpha(i,j-1)) * tmp_points[i][0];
                        tmp_points[i][1] = alpha(i,j-1) * tmp_points[i+1][1] + (1.0 - alpha(i,j-1)) * tmp_points[i][1];
                        tmp_points[i][2] = alpha(i,j-1) * tmp_points[i+1][2] + (1.0 - alpha(i,j-1)) * tmp_points[i][2];
                    }
                    IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, L+1, column, depth);
                    new_points(cp_index) = NodeType::Pointer(new NodeType(id, tmp_points[0][0], tmp_points[0][1], tmp_points[0][2]));
                    id++;

                    int index = polynomial_degree_u - j;
                    // std::cout << "index2: " << pole_u+r-j << std::endl;
                    // std::cout << "value: " << tmp_points[index][0] << std::endl;
                    cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, pole_u+r-j+1, column, depth);
                    new_points(cp_index) = NodeType::Pointer(new NodeType(id, tmp_points[index][0], tmp_points[index][1], tmp_points[index][2]));
                    id++;
                }
                for( IndexType i = L + 1; i < pole_u; ++i){
                    IndexType cp_index = NurbsUtilities::GetVectorIndexFromMatrixIndices(
                        new_num_of_cp_u, old_num_of_cp_v, old_num_of_cp_w, i+1, column, depth);
                    new_points(cp_index) = NodeType::Pointer(new NodeType(id, tmp_points[i-L][0], tmp_points[i-L][1], tmp_points[i-L][2]));
                    id++;
                }
            }
        }

        rGeometry = NurbsVolumeGeometry<PointerVector<NodeType>>(
            new_points, polynomial_degree_u, polynomial_degree_v, polynomial_degree_w,
            new_knots_u, old_knots_v, old_knots_w);
    }

    static void insert_knot_u(NurbsVolumeGeometryType& rGeometry, double u, SizeType r){

        // Get current order
        SizeType polynomial_degree_u = rGeometry.PolynomialDegreeU();
        // Get current knot information
        const Kratos::Vector& old_knots_u = rGeometry.KnotsU();
        SizeType old_num_of_knots_u = rGeometry.NumberOfKnotsU();
        // Get current cp's information
        SizeType old_num_of_cp_u = rGeometry.NumberOfControlPointsU();

        // Initial multiplicity
        // Todo: Add the capability to insert knots that already exist.
        SizeType s = 0;

        // Initialize new containers
        SizeType new_num_of_knots_u = old_num_of_knots_u + r;
        Kratos::Vector new_knots_u(new_num_of_knots_u);

        SizeType new_num_of_cp_u = old_num_of_cp_u + r;
        PointerVector<NodeType> new_points(new_num_of_cp_u);

        // Get Index of the lower knotspan [pole_u ... u ... span_end]
        SizeType pole_u = NurbsUtilities::GetLowerSpan(polynomial_degree_u, old_knots_u, u);

        // Create new ordered knot vector. Insert u r-times after pole_u.
        for( IndexType i = 0; i <= pole_u; i++)
            new_knots_u[i] = old_knots_u[i];
        for( IndexType i = 1; i <= r; ++i)
            new_knots_u[pole_u+i] = u;
        for( IndexType i = pole_u+1; i < old_num_of_knots_u; ++i)
            new_knots_u[i+r] = old_knots_u[i];

        // Copy unaltered control points.
        // Note: The +1 accounts for the fact that we have only p-multiplicity at the knot ends.
        SizeType id = 1;
        for( IndexType i = 0; i <= pole_u - polynomial_degree_u + 1; ++i){
            new_points(i) = NodeType::Pointer(new NodeType(id, rGeometry[i][0], rGeometry[i][1], rGeometry[i][2]));
            id += 1;
        }
        // The +1 accounts for the fact that we have only p-multiplicity at the knot ends.
        for( IndexType i = pole_u + 1; i < old_num_of_cp_u; ++i){
            new_points(i+r) = NodeType::Pointer(new NodeType(id, rGeometry[i][0], rGeometry[i][1], rGeometry[i][2]));
            id += 1;
        }

        // Temporary store cp's which will be replaces
        // P - s - 1 will be replaced
        // Starting at k - p + 1 (k:pole_u)
        // Note: The +1 accounts for the fact that we have only p-multiplicity at the knot ends.
        PointerVector<NodeType> tmp_points(new_num_of_cp_u);
        for( IndexType i = 0; i <= polynomial_degree_u - s; ++i){
            int index = pole_u - polynomial_degree_u + i + 1;
            tmp_points(i) = rGeometry(index);
        }

        // Find new cp's
        IndexType L;
        for( IndexType j = 1; j <= r; ++j){ // Insert knot r-times
            L = pole_u - polynomial_degree_u + j;
            for( IndexType i = 0; i <= polynomial_degree_u - j - s; i++){
                double alpha = (u - old_knots_u[L+i]) / (old_knots_u[i+pole_u+1] - old_knots_u[L+i]);
                tmp_points[i][0] = alpha * tmp_points[i+1][0] + (1.0 - alpha) * tmp_points[i][0];
                tmp_points[i][1] = alpha * tmp_points[i+1][1] + (1.0 - alpha) * tmp_points[i][1];
                tmp_points[i][2] = alpha * tmp_points[i+1][2] + (1.0 - alpha) * tmp_points[i][2];
            }
            std::cout << "index1: " << L << std::endl;
            std::cout << "value: " << tmp_points[0][0] << std::endl;
            new_points(L+1) = NodeType::Pointer(new NodeType(id, tmp_points[0][0], tmp_points[0][1], tmp_points[0][2]));
            id++;
            int index = polynomial_degree_u - j;
            std::cout << "index2: " << pole_u+r-j << std::endl;
            std::cout << "value: " << tmp_points[index][0] << std::endl;
            new_points(pole_u+r-j+1) = NodeType::Pointer(new NodeType(id, tmp_points[index][0], tmp_points[index][1], tmp_points[index][2]));
            id++;
        }
        for( IndexType i = L + 1; i < pole_u; ++i){
            std::cout << "index3: " << i << std::endl;
            std::cout << "value: " << tmp_points[i-L][0] << std::endl;
            new_points(i+1) = NodeType::Pointer(new NodeType(id, tmp_points[i-L][0], tmp_points[i-L][1], tmp_points[i-L][2]));
            id++;
        }
        for( int i = 0; i < new_num_of_cp_u; ++i){
            std::cout << new_points[i] << std::endl;
        }

        std::cout << "new knots u: \n" << new_knots_u << std::endl;
    }
    ///@}

};

} // End Namespace

#endif