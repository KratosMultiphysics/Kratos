//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_IGA_TEST_CREATION_UTILITY )
#define  KRATOS_IGA_TEST_CREATION_UTILITY

#include "geometries/geometry.h"
#include "geometries/nurbs_surface_geometry.h"

namespace Kratos
{
namespace TestCreationUtility
{
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node<3> NodeType;
    typedef PointerVector<NodeType> NodeVector;
    typedef Geometry<NodeType> GeometryType;
    typedef NurbsSurfaceGeometry<3, NodeVector> NurbsSurfaceType;

    ///@}
    ///@name Operations
    ///@{

    void AddDisplacementDofs(ModelPart& rModelPart) {
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }
    }

    NurbsSurfaceType GenerateNurbsSurface(ModelPart& rModelPart) {
        NodeVector points(6);

        points(0) = rModelPart.CreateNewNode(1,  0, 5,  0);
        points(1) = rModelPart.CreateNewNode(2,  5, 5,  0);
        points(2) = rModelPart.CreateNewNode(3, 10, 5, -4);
        points(3) = rModelPart.CreateNewNode(4,  0, 0,  0);
        points(4) = rModelPart.CreateNewNode(5,  5, 0,  0);
        points(5) = rModelPart.CreateNewNode(6, 10, 0, -4);

        Vector knot_u = ZeroVector(4);
        knot_u[0] = 0.0;
        knot_u[1] = 0.0;
        knot_u[2] = 10.0;
        knot_u[3] = 10.0;
        Vector knot_v = ZeroVector(2);
        knot_v[0] = 0.0;
        knot_v[1] = 5.0;

        int p = 2;
        int q = 1;

        return NurbsSurfaceType(
            points, p, q, knot_u, knot_v);
    }

    typename Geometry<NodeType>::Pointer GetQuadraturePointGeometry(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        typename GeometryType::IntegrationPointsArrayType integration_points(1);
        integration_points[0] = IntegrationPoint;
        typename GeometryType::GeometriesArrayType result_geometries;
        if (PolynomialDegree == SizeType(2)) {
            GenerateNurbsSurface(rModelPart).CreateQuadraturePointGeometries(
                result_geometries, 3, integration_points);
        }
        return result_geometries(0);
    }

}; /// namespace TestCreationUtility
} /// namespace Kratos.

#endif /// KRATOS_IGA_TEST_CREATION_UTILITY
