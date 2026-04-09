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
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"

#include "iga_application_variables.h"

namespace Kratos
{
namespace TestCreationUtility
{
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node NodeType;
    typedef PointerVector<NodeType> NodeVector;
    typedef Geometry<NodeType> GeometryType;
    typedef NurbsSurfaceGeometry<3, NodeVector> NurbsSurfaceType;
    typedef NurbsVolumeGeometry<NodeVector> NurbsVolumeType;

    ///@}
    ///@name Operations
    ///@{

    inline void AddDisplacementDofs(ModelPart& rModelPart) {
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }
    }

    inline void AddDirectorInc2DDofs(ModelPart& rModelPart) {
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DIRECTORINC_X);
            r_node.AddDof(DIRECTORINC_Y);
        }
    }

    inline typename NurbsSurfaceType::Pointer GenerateNurbsSurface(
        ModelPart& rModelPart, 
        SizeType PolynomialDegree)
    {
        SizeType p = PolynomialDegree;
        SizeType q = 1;

        Vector knot_u = ZeroVector(2*(p+1));
        for(IndexType index = knot_u.size()/2; index < knot_u.size(); ++index) {
            knot_u[index] = 1.0;
        }

        array_1d<double, 4> knot_v;

        knot_v[0] = 0.0;
        knot_v[1] = 0.0;
        knot_v[2] = 1.0;
        knot_v[3] = 1.0;
        
        NodeVector points((p+1)*(q+1));

        if (p == 3)
        {
            points(0) = rModelPart.CreateNewNode(1, 0.0, -0.05, 0.0);
            points(1) = rModelPart.CreateNewNode(2, 0.333333333333333, -0.05, 0.0);
            points(2) = rModelPart.CreateNewNode(3, 0.666666666666667, -0.05, 0.0);
            points(3) = rModelPart.CreateNewNode(4, 1.0, -0.05, 0.0);
            points(4) = rModelPart.CreateNewNode(5, 0.0, 0.05, 0.0);
            points(5) = rModelPart.CreateNewNode(6, 0.333333333333333, 0.05, 0.0);
            points(6) = rModelPart.CreateNewNode(7, 0.666666666666667, 0.05, 0.0);
            points(7) = rModelPart.CreateNewNode(8, 1.0, 0.05, 0.0);
        }
        else if (p == 4)
        {
            points(0) = rModelPart.CreateNewNode(1, 0.0, -0.05, 0.0);
            points(1) = rModelPart.CreateNewNode(2, 0.25, -0.05, 0.0);
            points(2) = rModelPart.CreateNewNode(3, 0.5, -0.05, 0.0);
            points(3) = rModelPart.CreateNewNode(4, 0.75, -0.05, 0.0);
            points(4) = rModelPart.CreateNewNode(5, 1.0, -0.05, 0.0);
            points(5) = rModelPart.CreateNewNode(6, 0.0, 0.05, 0.0);
            points(6) = rModelPart.CreateNewNode(7, 0.25, 0.05, 0.0);
            points(7) = rModelPart.CreateNewNode(8, 0.5, 0.05, 0.0);
            points(8) = rModelPart.CreateNewNode(9, 0.75, 0.05, 0.0);
            points(9) = rModelPart.CreateNewNode(10, 1.0, 0.05, 0.0);
        }
        else if (p == 5)
        {
            points(0) = rModelPart.CreateNewNode(1, 0.0, -0.05, 0.0);
            points(1) = rModelPart.CreateNewNode(2, 0.2, -0.05, 0.0);
            points(2) = rModelPart.CreateNewNode(3, 0.4, -0.05, 0.0);
            points(3) = rModelPart.CreateNewNode(4, 0.6, -0.05, 0.0);
            points(4) = rModelPart.CreateNewNode(5, 0.8, -0.05, 0.0);
            points(5) = rModelPart.CreateNewNode(6, 1.0, -0.05, 0.0);
            points(6) = rModelPart.CreateNewNode(7, 0.0, 0.05, 0.0);
            points(7) = rModelPart.CreateNewNode(8, 0.2, 0.05, 0.0);
            points(8) = rModelPart.CreateNewNode(9, 0.4, 0.05, 0.0);
            points(9) = rModelPart.CreateNewNode(10, 0.6, 0.05, 0.0);
            points(10) = rModelPart.CreateNewNode(11, 0.8, 0.05, 0.0);
            points(11) = rModelPart.CreateNewNode(12, 1.0, 0.05, 0.0);
        }

        return Kratos::make_shared<NurbsSurfaceType>(
            points, p, q, knot_u, knot_v);
    }

    inline typename Geometry<NodeType>::Pointer GetQuadraturePointGeometry(
        ModelPart& rModelPart, 
        SizeType PolynomialDegree, 
        IntegrationPoint<3> IntegrationPoint)
    {
        typename GeometryType::IntegrationPointsArrayType integration_points(1);
        integration_points[0] = IntegrationPoint;
        typename GeometryType::GeometriesArrayType result_geometries;

        auto p_nurbs_surface = GenerateNurbsSurface(rModelPart, PolynomialDegree);
        p_nurbs_surface->SetId(1);
        IntegrationInfo integration_info = p_nurbs_surface->GetDefaultIntegrationInfo();
        p_nurbs_surface->CreateQuadraturePointGeometries(
                result_geometries, 3, integration_points, integration_info);
        rModelPart.AddGeometry(p_nurbs_surface);

        return result_geometries(0);
    }


    inline typename Geometry<NodeType>::Pointer GetQuadraturePointGeometryOnCurve(
        ModelPart& rModelPart,
        SizeType PolynomialDegree,
        IntegrationPoint<3> IntegrationPoint)
    {
        typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
        IntegrationPointsArrayType integration_points(1);
        integration_points[0] = IntegrationPoint;
        typedef typename Geometry<NodeType>::GeometriesArrayType GeometriesArrayType;

        GeometriesArrayType result_geometries;

        // Create the BrepCurve in Surface
        // Assign the points belonging to the curve
        PointerVector<Node> points_curve;
        points_curve.push_back(Node::Pointer(new Node(1, 0.0, 0.05)));
        points_curve.push_back(Node::Pointer(new Node(2, 1.0, 0.05)));
        // Assign the curve's knot vector
        Vector knot_vector_curve = ZeroVector(4);
        knot_vector_curve[0] = 0.0;
        knot_vector_curve[1] = 0.0;
        knot_vector_curve[2] = 1.0;
        knot_vector_curve[3] = 1.0;
        // Polynomial degree of the curve
        int p_curve = 1;
        // Create the 2D embedded curve
        auto curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Node>>>(points_curve, p_curve, knot_vector_curve);

        auto p_nurbs_surface = GenerateNurbsSurface(rModelPart, 3);
        p_nurbs_surface->SetId(1);
        IntegrationInfo integration_info_surface = p_nurbs_surface->GetDefaultIntegrationInfo();
        p_nurbs_surface->CreateQuadraturePointGeometries(
                result_geometries, 3, integration_points, integration_info_surface);

        auto p_curve_nurbs_surface = Kratos::make_shared<NurbsCurveOnSurfaceGeometry<3, PointerVector<NodeType>, PointerVector<NodeType>>>(p_nurbs_surface, curve);

        p_curve_nurbs_surface->SetId(2);
        IntegrationInfo integration_info = p_curve_nurbs_surface->GetDefaultIntegrationInfo();
        
        p_curve_nurbs_surface->CreateQuadraturePointGeometries(result_geometries, 3, integration_points, integration_info);

        rModelPart.AddGeometry(p_nurbs_surface);
        rModelPart.AddGeometry(p_curve_nurbs_surface);

        return result_geometries(0);
    }

    inline typename NurbsVolumeType::Pointer GenerateRectangularNurbsVolumeP2(
        ModelPart& rModelPart,
        const double LengthX = 2.0,
        const double LengthY = 1.0,
        const double LengthZ = 0.5)
    {
        constexpr SizeType polynomial_degree = 2;
        Vector knot_u = ZeroVector(2 * (polynomial_degree + 1));
        Vector knot_v = ZeroVector(2 * (polynomial_degree + 1));
        Vector knot_w = ZeroVector(2 * (polynomial_degree + 1));

        for (IndexType i = polynomial_degree + 1; i < knot_u.size(); ++i) {
            knot_u[i] = 1.0;
            knot_v[i] = 1.0;
            knot_w[i] = 1.0;
        }

        NodeVector points((polynomial_degree + 1) * (polynomial_degree + 1) * (polynomial_degree + 1));

        IndexType node_id = 1;
        for (IndexType k = 0; k < polynomial_degree + 1; ++k) {
            for (IndexType j = 0; j < polynomial_degree + 1; ++j) {
                for (IndexType i = 0; i < polynomial_degree + 1; ++i) {
                    const double x = LengthX * static_cast<double>(i) / static_cast<double>(polynomial_degree);
                    const double y = LengthY * static_cast<double>(j) / static_cast<double>(polynomial_degree);
                    const double z = LengthZ * static_cast<double>(k) / static_cast<double>(polynomial_degree);
                    points(node_id - 1) = rModelPart.CreateNewNode(node_id, x, y, z);
                    ++node_id;
                }
            }
        }

        return Kratos::make_shared<NurbsVolumeType>(
            points, polynomial_degree, polynomial_degree, polynomial_degree, knot_u, knot_v, knot_w);
    }

    inline typename Geometry<NodeType>::Pointer GetQuadraturePointGeometryFromRectangularVolumeP2(
        ModelPart& rModelPart,
        IntegrationPoint<3> IntegrationPoint)
    {
        typename GeometryType::IntegrationPointsArrayType integration_points(1);
        integration_points[0] = IntegrationPoint;
        typename GeometryType::GeometriesArrayType result_geometries;

        auto p_nurbs_volume = GenerateRectangularNurbsVolumeP2(rModelPart);
        p_nurbs_volume->SetId(1);
        IntegrationInfo integration_info = p_nurbs_volume->GetDefaultIntegrationInfo();
        p_nurbs_volume->CreateQuadraturePointGeometries(
            result_geometries, 2, integration_points, integration_info);
        rModelPart.AddGeometry(p_nurbs_volume);

        return result_geometries(0);
    }

    ///@}
}; /// namespace TestCreationUtility
} /// namespace Kratos.

#endif /// KRATOS_IGA_TEST_CREATION_UTILITY
