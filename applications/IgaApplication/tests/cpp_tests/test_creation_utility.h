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

    NurbsSurfaceType GenerateNurbsSurface(ModelPart& rModelPart, SizeType PolynomialDegree) {
        
        SizeType p = PolynomialDegree;
        SizeType q = 1;

        Vector knot_u = ZeroVector(2*(p+1));
        for(IndexType index = 0; index < knot_u.size()/2; ++index) {
            knot_u[index] = 0.0;
        }
        for(IndexType index = knot_u.size()/2; index < knot_u.size(); ++index) {
            knot_u[index] = 1.0;
        }

        Vector knot_v = ZeroVector(2*(q+1));
        for(IndexType index = 0; index < knot_v.size()/2; ++index) {
            knot_v[index] = 0.0;
        }
        for(IndexType index = knot_v.size()/2; index < knot_v.size(); ++index) {
            knot_v[index] = 1.0;
        }

        NodeVector points((p+1)*(q+1));
        for(IndexType index = 0; index < points.size()/2; ++index) {
            double x = 1.0*index/PolynomialDegree;
            double y = -0.05;
            double z = 0.0;
            points(index) = rModelPart.CreateNewNode(index + 1, x, y, z);
        }
        for(IndexType index = points.size()/2; index < points.size(); ++index) {
            double x = 1.0*(index - PolynomialDegree -1)/PolynomialDegree;
            double y = 0.05;
            double z = 0.0;
            points(index) = rModelPart.CreateNewNode(index + 1, x, y, z);
        }

        return NurbsSurfaceType(
            points, p, q, knot_u, knot_v);
    }

    typename Geometry<NodeType>::Pointer GetQuadraturePointGeometry(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        typename GeometryType::IntegrationPointsArrayType integration_points(1);
        integration_points[0] = IntegrationPoint;
        typename GeometryType::GeometriesArrayType result_geometries;
        
        GenerateNurbsSurface(rModelPart,PolynomialDegree).CreateQuadraturePointGeometries(
                result_geometries, 3, integration_points);
        
        return result_geometries(0);
    }

    ///@}
}; /// namespace TestCreationUtility
} /// namespace Kratos.

#endif /// KRATOS_IGA_TEST_CREATION_UTILITY
