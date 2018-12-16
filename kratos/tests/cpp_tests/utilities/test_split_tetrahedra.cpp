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

KRATOS_TEST_CASE_IN_SUITE(SplitTetrahedra, KratosCoreFastSuite)
{
    std::vector<Node<3>::Pointer> nodes_pointer_vect(11);

    // Tetrahedra nodes
    nodes_pointer_vect[0] = Kratos::make_shared<Node<3>>(1, 0.0, 0.0, 0.0);
    nodes_pointer_vect[1] = Kratos::make_shared<Node<3>>(2, 1.0, 0.0, 0.0);
    nodes_pointer_vect[2] = Kratos::make_shared<Node<3>>(3, 0.0, 1.0, 0.0);
    nodes_pointer_vect[3] = Kratos::make_shared<Node<3>>(4, 0.0, 0.0, 1.0);

    // Edge intersection candidate nodes
    nodes_pointer_vect[4] = Kratos::make_shared<Node<3>>(5, 0.5, 0.0, 0.0); // Edge 01
    nodes_pointer_vect[5] = Kratos::make_shared<Node<3>>(6, 0.0, 0.5, 0.0); // Edge 02
    nodes_pointer_vect[6] = Kratos::make_shared<Node<3>>(7, 0.0, 0.0, 0.5); // Edge 03
    nodes_pointer_vect[7] = Kratos::make_shared<Node<3>>(8, 0.5, 0.5, 0.0); // Edge 12
    nodes_pointer_vect[8] = Kratos::make_shared<Node<3>>(9, 0.5, 0.0, 0.5); // Edge 13
    nodes_pointer_vect[9] = Kratos::make_shared<Node<3>>(10, 0.0, 0.5, 0.5); // Edge 23

    // Baricenter
    nodes_pointer_vect[10] = Kratos::make_shared<Node<3>>(11, 1.0/3.0, 1.0/3.0, 1.0/3.0); // Baricenter

    // Loop all the splitting patterns
    std::vector<std::vector<unsigned int>> a(6);
    a[0] = {0,1,4}; // Edge 01 candidate nodes
    a[1] = {0,2,5}; // Edge 02 candidate nodes 
    a[2] = {0,3,6}; // Edge 03 candidate nodes
    a[3] = {1,2,7}; // Edge 12 candidate nodes
    a[4] = {1,3,8}; // Edge 13 candidate nodes
    a[5] = {2,3,9}; // Edge 23 candidate nodes

    int n_elems;
    int steiner_node;
    int n_split_edges;
    std::vector<int> t(56);
    std::vector<int> edge_int_vect(6);

    unsigned int counter = 0;
    for(auto i0 : a[0]) {
        for(auto i1 : a[1]) {
            for(auto i2 : a[2]) {
                for(auto i3 : a[3]) {
                    for(auto i4 : a[4]) {
                        for(auto i5 : a[5]) {
                            counter++;
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
                            // std::cout << n_elems << " ";
                            std::cout << "Permutation: " << counter << " edges: ";
                            for(auto k : edge_int_vect)
                                std::cout << k << " ";
                            std::cout << std::endl;
                            std::cout << " n_elems: " << n_elems << std::endl;
                            std::cout << " n_split_edges: " << n_split_edges << std::endl;

                            std::cout << " t: ";
                            for (unsigned int i = 0; i < t.size(); ++i) {
                                std::cout << t[i] << " ";
                            }
                            std::cout << std::endl;

                            // Reconstruct each one of the splitting tetras
                            for (unsigned int i_elem = 0; i_elem < n_elems; ++i_elem) {
                                std::cout << " i_elem: " << i_elem 
                                << " i0: " << (*(nodes_pointer_vect[t[i_elem*4]])).Id()
                                << " i1: " << (*(nodes_pointer_vect[t[i_elem*4 + 1]])).Id()
                                << " i2: " << (*(nodes_pointer_vect[t[i_elem*4 + 2]])).Id()
                                << " i3: " << (*(nodes_pointer_vect[t[i_elem*4 + 3]])).Id() << std::endl;
                                Tetrahedra3D4<Node<3>>::Pointer p_sub_tetra = Kratos::make_shared<Tetrahedra3D4<Node<3>>>(
                                    nodes_pointer_vect[t[i_elem*4]],
                                    nodes_pointer_vect[t[i_elem*4 + 1]],
                                    nodes_pointer_vect[t[i_elem*4 + 2]],
                                    nodes_pointer_vect[t[i_elem*4 + 3]]);

                                // Check subtetra volumes
                                sub_vol = p_sub_tetra->Volume();
                                tot_vol += sub_vol;

                                std::cout << " Sub. vol " << i_elem << " : " << sub_vol << std::endl;

                                KRATOS_ERROR_IF(sub_vol < 0.0) << "Negative sub tetra volume for edges: " << 
                                    edge_int_vect[0] << " " <<
                                    edge_int_vect[1] << " " <<
                                    edge_int_vect[2] << " " <<
                                    edge_int_vect[3] << " " <<
                                    edge_int_vect[4] << " " <<
                                    edge_int_vect[5] << " " << std::endl;
                            }

                            std::cout << " Tot. vol: " << tot_vol << std::endl;

                            KRATOS_CHECK_NEAR(tot_vol, 1.0/6.0, 1e-4);
                        }
                    }
                }
            }
        }
    }


}

}   // namespace Testing
}  // namespace Kratos.


