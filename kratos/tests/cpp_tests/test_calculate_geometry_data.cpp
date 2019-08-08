//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/geometry_utilities.h"

#include "utilities/builtin_timer.h"

/* Utilities */

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;
        typedef std::size_t IndexSize;

        KRATOS_TEST_CASE_IN_SUITE(TestCalculateGeometryData, KratosCoreFastSuite)
        {
            //auto node_1(Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0));
            //auto node_2(Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0));
            //auto node_3(Kratos::make_intrusive<NodeType>(3, 0.5, 1.0, 0.0));
            //auto node_4(Kratos::make_intrusive<NodeType>(4, 0.5, 0.3, 1.0));

            //const Geometry<NodeType>::Pointer p_geom(Kratos::make_shared<Tetrahedra3D4<NodeType>>(node_1, node_2, node_3, node_4));
            NodeType::Pointer p_node_1 = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
            NodeType::Pointer p_node_2 = Kratos::make_intrusive<Node<3>>(2, 1.0, 0.0, 0.0);
            NodeType::Pointer p_node_3 = Kratos::make_intrusive<Node<3>>(3, 1.0, 1.0, 0.0);

            // Now we create the geometry
            std::vector<NodeType::Pointer> triangle_nodes(3);
            triangle_nodes[0] = p_node_1;
            triangle_nodes[1] = p_node_2;
            triangle_nodes[2] = p_node_3;
            Triangle2D3 <NodeType> triangle(PointerVector<NodeType>{triangle_nodes});

            const std::size_t num_evaluations = 1e8;

            BuiltinTimer timer;

            array_1d<double, 3> N_all;
            for (IndexType i = 0; i<num_evaluations; ++i) {

                // Computing the info
                BoundedMatrix<double, 3, 2> DN_DX;
                array_1d<double, 3> N;
                double area;

                GeometryUtils::CalculateGeometryData(triangle, DN_DX, N, area);

                N_all += N;
            }

            KRATOS_WATCH(N_all)

            std::cout << std::endl << num_evaluations << " evaluations took " << timer.ElapsedSeconds() << std::endl;

            NodeType::Pointer p_node_5 = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
            NodeType::Pointer p_node_6 = Kratos::make_intrusive<Node<3>>(2, 1.0, 0.0, 0.0);
            NodeType::Pointer p_node_7 = Kratos::make_intrusive<Node<3>>(3, 1.0, 1.0, 0.0);
            NodeType::Pointer p_node_8 = Kratos::make_intrusive<Node<3>>(4, 1.0, 1.0, 1.0);

            // Now we create the geometry
            std::vector<NodeType::Pointer> tetrahedra_nodes(4);
            tetrahedra_nodes[0] = p_node_5;
            tetrahedra_nodes[1] = p_node_6;
            tetrahedra_nodes[2] = p_node_7;
            tetrahedra_nodes[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra(PointerVector<NodeType>{tetrahedra_nodes});

            BuiltinTimer timer2;

            array_1d<double, 4> N_all2;
            for (IndexType i = 0; i < num_evaluations; ++i) {

                // Computing the info
                BoundedMatrix<double, 4, 3> DN_DX;
                array_1d<double, 4> N;
                double volume;

                GeometryUtils::CalculateGeometryData(tetrahedra, DN_DX, N, volume);

                N_all2 += N;
            }
            KRATOS_WATCH(N_all2)

            std::cout << std::endl << num_evaluations << " evaluations took " << timer2.ElapsedSeconds() << std::endl;

        }

    } // namespace Testing
}  // namespace Kratos.