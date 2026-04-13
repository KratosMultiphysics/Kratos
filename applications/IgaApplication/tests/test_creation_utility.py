import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA

class TestCreationUtility:

    @staticmethod
    def GenerateNurbsVolumeP2Rectangular(model_part, length_x=2.0, length_y=1.0, length_z=0.5):
        p = 2
        q = 2
        r = 2

        knot_u = KM.Vector(6)
        knot_v = KM.Vector(6)
        knot_w = KM.Vector(6)
        for knot_vector in (knot_u, knot_v, knot_w):
            knot_vector[0] = 0.0
            knot_vector[1] = 0.0
            knot_vector[2] = 0.0
            knot_vector[3] = 1.0
            knot_vector[4] = 1.0
            knot_vector[5] = 1.0

        points = KM.NodesVector()

        node_id = 1
        for k in range(3):
            for j in range(3):
                for i in range(3):
                    x = length_x * i / 2.0
                    y = length_y * j / 2.0
                    z = length_z * k / 2.0
                    points.append(model_part.CreateNewNode(node_id, x, y, z))
                    node_id += 1

        return KM.NurbsVolumeGeometry(points, p, q, r, knot_u, knot_v, knot_w)

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

    @staticmethod
    def GetQuadraturePointGeometryFromRectangularVolumeP2(model_part, integration_point):
        volume = TestCreationUtility.GenerateNurbsVolumeP2Rectangular(model_part)
        volume.SetId(1)
        geom_vector = KM.GeometriesVector()
        # Python bindings do not expose the CUSTOM IntegrationInfo path needed to
        # request a single user-defined quadrature point with second derivatives.
        # Use the default quadrature-point geometries and return the first one.
        volume.CreateQuadraturePointGeometries(geom_vector, 3)
        model_part.AddGeometry(volume)

        return geom_vector[0]

    @staticmethod
    def GetQuadraturePointGeometryOnVolumeSurfaceP2(model_part, integration_point):
        volume = TestCreationUtility.GenerateNurbsVolumeP2Rectangular(model_part)
        volume.SetId(1)
        model_part.AddGeometry(volume)

        p1 = KM.Node(1001, 0.0, 0.0, 0.0)
        p2 = KM.Node(1002, 1.0, 0.0, 0.0)
        p3 = KM.Node(1003, 1.0, 1.0, 0.0)
        p4 = KM.Node(1004, 0.0, 1.0, 0.0)
        face_in_volume = KM.Quadrilateral3D4(p1, p2, p3, p4)

        surface_in_volume = KM.SurfaceInNurbsVolumeGeometry(volume, face_in_volume)
        surface_in_volume.SetId(2)
        model_part.AddGeometry(surface_in_volume)

        geom_vector = KM.GeometriesVector()
        surface_in_volume.CreateQuadraturePointGeometries(geom_vector, 3, [integration_point])

        return geom_vector[0]
