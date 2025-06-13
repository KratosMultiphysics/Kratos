import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA

class TestCreationUtility:

    @staticmethod
    def GenerateNurbsSurface(model_part, polynomial_degree):
        p = polynomial_degree
        q = 1  
        
        knot_u = KM.Vector(2 * (p + 1))
        for i in range((p + 1), 2 * (p + 1)):
            knot_u[i] = 1.0

        knot_v = KM.Vector(4)
        knot_v[0] = 0.0
        knot_v[1] = 0.0
        knot_v[2] = 1.0
        knot_v[3] = 1.0

        points = KM.NodesVector()

        if p == 3:
            points.append(model_part.CreateNewNode(1, 0.0, -0.05, 0.0))
            points.append(model_part.CreateNewNode(2, 0.333333333333333, -0.05, 0.0))
            points.append(model_part.CreateNewNode(3, 0.666666666666667, -0.05, 0.0))
            points.append(model_part.CreateNewNode(4, 1.0, -0.05, 0.0))
            points.append(model_part.CreateNewNode(5, 0.0, 0.05, 0.0))
            points.append(model_part.CreateNewNode(6, 0.333333333333333, 0.05, 0.0))
            points.append(model_part.CreateNewNode(7, 0.666666666666667, 0.05, 0.0))
            points.append(model_part.CreateNewNode(8, 1.0, 0.05, 0.0))

        elif p == 4:
            points.append(model_part.CreateNewNode(1, 0.0, -0.05, 0.0))
            points.append(model_part.CreateNewNode(2, 0.25, -0.05, 0.0))
            points.append(model_part.CreateNewNode(3, 0.5, -0.05, 0.0))
            points.append(model_part.CreateNewNode(4, 0.75, -0.05, 0.0))
            points.append(model_part.CreateNewNode(5, 1.0, -0.05, 0.0))
            points.append(model_part.CreateNewNode(6, 0.0, 0.05, 0.0))
            points.append(model_part.CreateNewNode(7, 0.25, 0.05, 0.0))
            points.append(model_part.CreateNewNode(8, 0.5, 0.05, 0.0))
            points.append(model_part.CreateNewNode(9, 0.75, 0.05, 0.0))
            points.append(model_part.CreateNewNode(10, 1.0, 0.05, 0.0))

        elif p == 5:
            points.append(model_part.CreateNewNode(1, 0.0, -0.05, 0.0))
            points.append(model_part.CreateNewNode(2, 0.2, -0.05, 0.0))
            points.append(model_part.CreateNewNode(3, 0.4, -0.05, 0.0))
            points.append(model_part.CreateNewNode(4, 0.6, -0.05, 0.0))
            points.append(model_part.CreateNewNode(5, 0.8, -0.05, 0.0))
            points.append(model_part.CreateNewNode(6, 1.0, -0.05, 0.0))
            points.append(model_part.CreateNewNode(7, 0.0, 0.05, 0.0))
            points.append(model_part.CreateNewNode(8, 0.2, 0.05, 0.0))
            points.append(model_part.CreateNewNode(9, 0.4, 0.05, 0.0))
            points.append(model_part.CreateNewNode(10, 0.6, 0.05, 0.0))
            points.append(model_part.CreateNewNode(11, 0.8, 0.05, 0.0))
            points.append(model_part.CreateNewNode(12, 1.0, 0.05, 0.0))

        else:
            raise ValueError(f"Polynomial degree {p} not supported in this test utility.")

        surface = KM.NurbsSurfaceGeometry3D(points, p, q, knot_u, knot_v)
        return surface

    @staticmethod
    def GetQuadraturePointGeometry(model_part, polynomial_degree, integration_point):
        surface = TestCreationUtility.GenerateNurbsSurface(model_part, polynomial_degree)
        surface.SetId(1)
        geom_vector = KM.GeometriesVector()
        surface.CreateQuadraturePointGeometries(geom_vector, 3, [integration_point])
        model_part.AddGeometry(surface)
        
        return geom_vector[0]
    
    @staticmethod
    def GetQuadraturePointGeometryOnCurve(model_part, polynomial_degree, integration_point):
        # Create the embedded curve (parametric in 2D, embedded in surface)
        points_curve = KM.NodesVector()
        points_curve.append(KM.Node(1, 0.0, 0.05, 0.0))
        points_curve.append(KM.Node(2, 1.0, 0.05, 0.0))

        knot_vector_curve = KM.Vector(4)
        knot_vector_curve[0] = 0.0
        knot_vector_curve[1] = 0.0
        knot_vector_curve[2] = 1.0
        knot_vector_curve[3] = 1.0

        degree_curve = 1

        curve = KM.NurbsCurveGeometry2D(points_curve, degree_curve, knot_vector_curve)

        # Generate the NURBS surface
        surface = TestCreationUtility.GenerateNurbsSurface(model_part, 3)
        surface.SetId(1)
        model_part.AddGeometry(surface)

        # Create the curve on surface geometry
        curve_on_surface = KM.NurbsCurveOnSurfaceGeometry3D(surface, curve)
        curve_on_surface.SetId(2)
        model_part.AddGeometry(curve_on_surface)

        # Convert integration point to expected input (list of lists with [xi, eta, zeta, weight])

        geom_vector = KM.GeometriesVector()

        # Create quadrature point geometry (Python binding version does not take integration info)
        curve_on_surface.CreateQuadraturePointGeometries(geom_vector, 3, [integration_point])

        return geom_vector[0]