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
<<<<<<< HEAD
        ModelPart& rModelPart,
        SizeType PolynomialDegree,
        IntegrationPoint<3> IntegrationPoint)
=======
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
>>>>>>> 3D_sbm
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

<<<<<<< HEAD
=======
        // // Assign the points belonging to the surface
        // PointerVector<Node> points_surface;
        // points_surface.push_back(Node::Pointer(new Node( 0,  0,  0)));
        // points_surface.push_back(Node::Pointer(new Node( 1,  0,  0)));
        // points_surface.push_back(Node::Pointer(new Node( 2,  0,  0)));
        // points_surface.push_back(Node::Pointer(new Node( 0,  1,  0)));
        // points_surface.push_back(Node::Pointer(new Node( 1,  1,  0)));
        // points_surface.push_back(Node::Pointer(new Node( 2,  2,  0)));
        // points_surface.push_back(Node::Pointer(new Node( 0,  2,  0)));
        // points_surface.push_back(Node::Pointer(new Node( 1,  2,  0)));
        // points_surface.push_back(Node::Pointer(new Node( 2,  2,  0)));

        // // Assign the surface's knot vectors
        // Vector knot_vector_u_surface = ZeroVector(4);
        // knot_vector_u_surface[0] = 0.0;
        // knot_vector_u_surface[1] = 0.0;
        // knot_vector_u_surface[2] = 2.0;
        // knot_vector_u_surface[3] = 2.0;

        // Vector knot_vector_v_surface = ZeroVector(4);
        // knot_vector_v_surface[0] = 0.0;
        // knot_vector_v_surface[1] = 0.0;
        // knot_vector_v_surface[2] = 1.0;
        // knot_vector_v_surface[3] = 1.0;

        // // Polynomial degrees
        // int p_surface = 2;
        // int q_surface = 2;

        // // Create a 3D surface
        // auto p_nurbs_surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_surface, p_surface,
        //     q_surface, knot_vector_u_surface, knot_vector_v_surface);

>>>>>>> 3D_sbm
        auto p_nurbs_surface = GenerateNurbsSurface(rModelPart, 3);
        p_nurbs_surface->SetId(1);
        IntegrationInfo integration_info_surface = p_nurbs_surface->GetDefaultIntegrationInfo();
        p_nurbs_surface->CreateQuadraturePointGeometries(
                result_geometries, 3, integration_points, integration_info_surface);

        auto p_curve_nurbs_surface = Kratos::make_shared<NurbsCurveOnSurfaceGeometry<3, PointerVector<NodeType>, PointerVector<NodeType>>>(p_nurbs_surface, curve);

        p_curve_nurbs_surface->SetId(2);
        IntegrationInfo integration_info = p_curve_nurbs_surface->GetDefaultIntegrationInfo();
<<<<<<< HEAD
=======

        // p_curve_nurbs_surface->CreateIntegrationPoints(integration_points, integration_info);
>>>>>>> 3D_sbm
        
        p_curve_nurbs_surface->CreateQuadraturePointGeometries(result_geometries, 3, integration_points, integration_info);

        rModelPart.AddGeometry(p_nurbs_surface);
        rModelPart.AddGeometry(p_curve_nurbs_surface);

        return result_geometries(0);
    }

    ///@}
}; /// namespace TestCreationUtility
} /// namespace Kratos.

#endif /// KRATOS_IGA_TEST_CREATION_UTILITY
