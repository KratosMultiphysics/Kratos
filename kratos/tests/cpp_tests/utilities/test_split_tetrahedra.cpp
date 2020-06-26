//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// Project includes
#include "includes/node.h"
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/split_tetrahedra.h"

namespace Kratos {
namespace Testing {

void SearchExternalFaces(
    const std::vector<Geometry<Node<3>>> &rSubTetrasFaces,
    std::vector<Geometry<Node<3>>> &rSubTetrasExtFaces)
{
    for (unsigned int i_face = 0; i_face < rSubTetrasFaces.size(); ++i_face) {
        unsigned int i_face_times = 0;
        const Geometry<Node<3>> &r_i_face = rSubTetrasFaces[i_face];
        const unsigned int i_0 = r_i_face[0].Id();
        const unsigned int i_1 = r_i_face[1].Id();
        const unsigned int i_2 = r_i_face[2].Id();

        // Check the current face against all faces
        for (unsigned int j_face = 0; j_face < rSubTetrasFaces.size(); ++j_face) {
            const Geometry<Node<3>> &r_j_face = rSubTetrasFaces[j_face];
            const unsigned int j_0 = r_j_face[0].Id();
            const unsigned int j_1 = r_j_face[1].Id();
            const unsigned int j_2 = r_j_face[2].Id();

            if (i_0 == j_0 || i_0 == j_1 || i_0 == j_2) {
                if (i_1 == j_0 || i_1 == j_1 || i_1 == j_2) {
                    if (i_2 == j_0 || i_2 == j_1 || i_2 == j_2) {
                        i_face_times++;
                    }
                }
            }
        }

        // If the face it's external, there must be a unique occurrence
        if (i_face_times == 1) {
            rSubTetrasExtFaces.push_back(r_i_face);
        }
    }
}

void PrintEdgesAndExternalFaces(
    const unsigned int Permutation,
    const std::vector<int> &rEdgeIntVect,
    const std::vector<Geometry<Node<3>>> rSubTetrasExtFaces)
{
    std::cout << "(("
        << rEdgeIntVect[0] << ","
        << rEdgeIntVect[1] << ","
        << rEdgeIntVect[2] << ","
        << rEdgeIntVect[3] << ","
        << rEdgeIntVect[4] << ","
        << rEdgeIntVect[5] << ") , [";

    for (unsigned int i_ext_face = 0; i_ext_face < rSubTetrasExtFaces.size(); i_ext_face++) {
        const Geometry<Node<3>> &r_i_ext_face = rSubTetrasExtFaces[i_ext_face];
        const unsigned int i_ext_0 = r_i_ext_face[0].Id();
        const unsigned int i_ext_1 = r_i_ext_face[1].Id();
        const unsigned int i_ext_2 = r_i_ext_face[2].Id();
        std::cout << "(" << i_ext_0 << "," << i_ext_1 << "," << i_ext_2 << ")";
        if (i_ext_face != (rSubTetrasExtFaces.size() - 1)) {
            std::cout << ",";
        }
    }

    if (Permutation != 729) {
        std::cout << "])," << std::endl;
    } else {
        std::cout << "])"  << std::endl;
    }
}

KRATOS_TEST_CASE_IN_SUITE(TetrahedraSplitModes, KratosCoreFastSuite)
{
    std::vector<int> aux_ids(11);
    std::vector<int> edge_ids(6);

    aux_ids[0] = 2; aux_ids[1] = 1; aux_ids[2] = 4; aux_ids[3] = 3; // Vertices
    aux_ids[4] = -1; // Edge 01
    aux_ids[5] =  1; // Edge 02 (split)
    aux_ids[6] = -1; // Edge 03
    aux_ids[7] =  1; // Edge 12 (split)
    aux_ids[8] = -1; // Edge 13
    aux_ids[9] =  1; // Edge 23 (split)
    aux_ids[10] = -1; // Internal (Steiner) node not required

    TetrahedraSplit::TetrahedraSplitMode(aux_ids.data(), edge_ids.data());

    const std::vector<int> expected_results{0, 5, 3, 7, 3, 9};
    for (std::size_t i = 0; i < 6; ++i) {
        KRATOS_CHECK_EQUAL(edge_ids[i], expected_results[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TetrahedraSplitEdgesPatterns, KratosCoreFastSuite)
{
    std::vector<Node<3>::Pointer> nodes_pointer_vect(11);

    // Tetrahedra nodes
    nodes_pointer_vect[0] = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
    nodes_pointer_vect[1] = Kratos::make_intrusive<Node<3>>(2, 1.0, 0.0, 0.0);
    nodes_pointer_vect[2] = Kratos::make_intrusive<Node<3>>(3, 0.0, 1.0, 0.0);
    nodes_pointer_vect[3] = Kratos::make_intrusive<Node<3>>(4, 0.0, 0.0, 1.0);

    // Edge intersection candidate nodes
    nodes_pointer_vect[4] = Kratos::make_intrusive<Node<3>>(5, 0.5, 0.0, 0.0); // Edge 01
    nodes_pointer_vect[5] = Kratos::make_intrusive<Node<3>>(6, 0.0, 0.5, 0.0); // Edge 02
    nodes_pointer_vect[6] = Kratos::make_intrusive<Node<3>>(7, 0.0, 0.0, 0.5); // Edge 03
    nodes_pointer_vect[7] = Kratos::make_intrusive<Node<3>>(8, 0.5, 0.5, 0.0); // Edge 12
    nodes_pointer_vect[8] = Kratos::make_intrusive<Node<3>>(9, 0.5, 0.0, 0.5); // Edge 13
    nodes_pointer_vect[9] = Kratos::make_intrusive<Node<3>>(10, 0.0, 0.5, 0.5); // Edge 23

    // Baricenter
    nodes_pointer_vect[10] = Kratos::make_intrusive<Node<3>>(11, 0.25, 0.25, 0.25); // Baricenter

    // Loop all the splitting patterns
    std::vector<std::vector<unsigned int>> a(6);
    a[0] = {0,1,4}; // Edge 01 candidate nodes
    a[1] = {0,2,5}; // Edge 02 candidate nodes
    a[2] = {0,3,6}; // Edge 03 candidate nodes
    a[3] = {1,2,7}; // Edge 12 candidate nodes
    a[4] = {1,3,8}; // Edge 13 candidate nodes
    a[5] = {2,3,9}; // Edge 23 candidate nodes

    // These are auxiliary bools used for debugging and for generating the patterns
    // Switch them to true to compute the split external faces and to print the permutation info
    const bool compute_faces = false;
    const bool print_edges_and_ext_faces = false;

    // Loop and check all the edges permutations
    int n_elems=0;
    int steiner_node;
    int n_split_edges;
    std::vector<int> t(56);
    std::vector<int> edge_int_vect(6);
    unsigned int permutation = 0;
    for(auto i0 : a[0]) {
        for(auto i1 : a[1]) {
            for(auto i2 : a[2]) {
                for(auto i3 : a[3]) {
                    for(auto i4 : a[4]) {
                        for(auto i5 : a[5]) {

                            permutation++;
                            edge_int_vect[0] = i0;
                            edge_int_vect[1] = i1;
                            edge_int_vect[2] = i2;
                            edge_int_vect[3] = i3;
                            edge_int_vect[4] = i4;
                            edge_int_vect[5] = i5;

                            double sub_vol = 0.0;
                            double tot_vol = 0.0;

                            // Call the split tetrahedra utility
                            TetrahedraSplit::Split_Tetrahedra(edge_int_vect.data(), t.data(), &n_elems, &n_split_edges, &steiner_node);

                            // Declare the faces array (only needed for debugging)
                            std::vector<Geometry<Node<3>>> sub_tetras_faces(n_elems * 4);
                            std::vector<Geometry<Node<3>>> sub_tetras_ext_faces;

                            // Reconstruct each one of the splitting tetras
                            for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
                                Tetrahedra3D4<Node<3>>::Pointer p_sub_tetra = Kratos::make_shared<Tetrahedra3D4<Node<3>>>(
                                    nodes_pointer_vect[t[i_elem*4]],
                                    nodes_pointer_vect[t[i_elem*4 + 1]],
                                    nodes_pointer_vect[t[i_elem*4 + 2]],
                                    nodes_pointer_vect[t[i_elem*4 + 3]]);

                                // Check subtetra volumes
                                sub_vol = p_sub_tetra->Volume();
                                tot_vol += sub_vol;

                                // If required, save sub tetra faces
                                if (compute_faces) {
                                    const auto faces = p_sub_tetra->GenerateFaces();
                                    for (unsigned int i_face = 0; i_face < 4; ++i_face) {
                                       sub_tetras_faces[i_elem * 4 + i_face] = faces[i_face];
                                    }
                                }

                                KRATOS_ERROR_IF(sub_vol < 1e-12) << "Negative subdivision " <<
                                    i_elem << " volume for edges: " <<
                                    edge_int_vect[0] << " " <<
                                    edge_int_vect[1] << " " <<
                                    edge_int_vect[2] << " " <<
                                    edge_int_vect[3] << " " <<
                                    edge_int_vect[4] << " " <<
                                    edge_int_vect[5] << std::endl;
                            }

                            KRATOS_CHECK_NEAR(tot_vol, 1.0/6.0, 1e-12);

                            // Search for the external faces
                            if (compute_faces) {
                                SearchExternalFaces(sub_tetras_faces, sub_tetras_ext_faces);
                            }

                            // Print splitting pattern and exterior faces
                            if (print_edges_and_ext_faces) {
                                PrintEdgesAndExternalFaces(permutation, edge_int_vect, sub_tetras_ext_faces);
                            }

                        }
                    }
                }
            }
        }
    }


}

}   // namespace Testing
}  // namespace Kratos.


