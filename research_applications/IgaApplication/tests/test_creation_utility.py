import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA

class TestCreationUtility:

    def GenerateNurbsSurfaceP2(model_part):
        p = 2
        q = 2
        
        knot_u = KM.Vector(6)
        knot_u[0] = 0.0
        knot_u[1] = 0.0
        knot_u[2] = 0.0
        knot_u[3] = 1.0
        knot_u[4] = 1.0
        knot_u[5] = 1.0

        knot_v = KM.Vector(6)
        knot_v[0] = 0.0
        knot_v[1] = 0.0
        knot_v[2] = 0.0
        knot_v[3] = 1.0
        knot_v[4] = 1.0
        knot_v[5] = 1.0

        points = KM.NodesVector()

        # 4x4 control points grid (u in rows, v in cols)
        id = 1
        for i in range(5):  # v direction
            for j in range(5):  # u direction
                x = j * 1.0
                y = i * 1.0
                z = 0.0
                points.append(model_part.CreateNewNode(id, x, y, z))
                id += 1

        surface = KM.NurbsSurfaceGeometry3D(points, p, q, knot_u, knot_v)
        return surface


    @staticmethod
    def GetQuadraturePointGeometryP2(model_part, integration_point):
        surface = TestCreationUtility.GenerateNurbsSurfaceP2(model_part)
        surface.SetId(1)
        geom_vector = KM.GeometriesVector()
        surface.CreateQuadraturePointGeometries(geom_vector, 3, [integration_point])
        model_part.AddGeometry(surface)
        
        return geom_vector[0]
    
    @staticmethod
    def GetQuadraturePointGeometryOnCurveP2(model_part, integration_point):
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
        surface = TestCreationUtility.GenerateNurbsSurfaceP2(model_part)
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