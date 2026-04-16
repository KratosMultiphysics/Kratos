import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication
from KratosMultiphysics import KratosUnittest


class TestPythonBindingsBrepSurface(KratosUnittest.TestCase):
    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("ModelPart")

        self.degree_u = 1
        self.degree_v = 1

        self.knot_vector_u = [0.0, 0.0, 1.0, 1.0]
        self.knot_vector_v = [0.0, 0.0, 1.0, 1.0]

        # 2x2 bilinear patch
        self.control_points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0]
        ]

        MakeNodesVectorFromCtrlPts(self.model_part, self.control_points)

    def test_create_brep_surface_from_python_bindings(self):
        weights = MakeVector([1.0, 1.0, 1.0, 1.0])

        brep_surface = CreateBrepSurfaceFromPythonIGABindings(
            self.model_part,
            self.degree_u,
            self.degree_v,
            self.knot_vector_u,
            self.knot_vector_v,
            weights
        )

        self.assertEqual(self.model_part.NumberOfNodes(), 4)
        self.assertEqual(CountGeometries(self.model_part), 1)

        self.assertEqual(brep_surface.Id, 1002)
        self.assertEqual(brep_surface.PolynomialDegree(0), self.degree_u)
        self.assertEqual(brep_surface.PolynomialDegree(1), self.degree_v)

    def test_create_brep_surface_from_python_bindings_with_default_weights(self):
        brep_surface = CreateBrepSurfaceFromPythonIGABindings(
            self.model_part,
            self.degree_u,
            self.degree_v,
            self.knot_vector_u,
            self.knot_vector_v
        )

        self.assertEqual(self.model_part.NumberOfNodes(), 4)
        self.assertEqual(CountGeometries(self.model_part), 1)
        self.assertEqual(brep_surface.Id, 1002)


def MakeVector(values):
    v = KM.Vector(len(values))
    for i, val in enumerate(values):
        v[i] = float(val)
    return v


def MakePointsList2D(p0, p1):
    return [
        KM.Point(float(p0[0]), float(p0[1]), 0.0),
        KM.Point(float(p1[0]), float(p1[1]), 0.0)
    ]


def MakeBoundaryBrepCurve(surface, p0, p1, same_curve_direction=True):
    knot_vector_curve = KM.Vector(4)
    knot_vector_curve[0] = 0.0
    knot_vector_curve[1] = 0.0
    knot_vector_curve[2] = 1.0
    knot_vector_curve[3] = 1.0

    degree_curve = 1
    points_curve = MakePointsList2D(p0, p1)

    curve_2d = KM.NurbsCurveGeometry2DPoint(points_curve, degree_curve, knot_vector_curve)
    return KM.BrepCurveOnSurface(surface, curve_2d, same_curve_direction)


def GetParametricBounds(knot_vector_u, knot_vector_v, degree_u, degree_v):
    u_min = float(knot_vector_u[degree_u])
    u_max = float(knot_vector_u[-degree_u - 1])
    v_min = float(knot_vector_v[degree_v])
    v_max = float(knot_vector_v[-degree_v - 1])
    return u_min, u_max, v_min, v_max


def MakeNodesVectorFromCtrlPts(model_part, control_points):
    points = KM.NodesVector()
    for i, cp in enumerate(control_points, start=1):
        node = model_part.CreateNewNode(i, float(cp[0]), float(cp[1]), float(cp[2]))
        points.append(node)
    return points


def CreateBrepSurfaceFromPythonIGABindings(
    model_part,
    degree_u,
    degree_v,
    knot_vector_u,
    knot_vector_v,
    weights=None
):
    surface_points = KM.NodesVector()
    for node_id in sorted(node.Id for node in model_part.Nodes):
        surface_points.append(model_part.GetNode(node_id))

    knots_u = MakeVector(knot_vector_u)
    knots_v = MakeVector(knot_vector_v)

    if weights is None:
        weights = KM.Vector(model_part.NumberOfNodes())
        for i in range(model_part.NumberOfNodes()):
            weights[i] = 1.0

    nurbs_surface = KM.NurbsSurfaceGeometry3D(
        surface_points,
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        weights
    )
    nurbs_surface.SetId(1001)

    u_min, u_max, v_min, v_max = GetParametricBounds(
        knot_vector_u, knot_vector_v, degree_u, degree_v
    )

    curve_bottom = MakeBoundaryBrepCurve(nurbs_surface, (u_min, v_min), (u_max, v_min), True)
    curve_right  = MakeBoundaryBrepCurve(nurbs_surface, (u_max, v_min), (u_max, v_max), True)
    curve_top    = MakeBoundaryBrepCurve(nurbs_surface, (u_max, v_max), (u_min, v_max), True)
    curve_left   = MakeBoundaryBrepCurve(nurbs_surface, (u_min, v_max), (u_min, v_min), True)

    outer_loops = [[curve_bottom, curve_right, curve_top, curve_left]]
    inner_loops = []

    brep_surface = KM.BrepSurface(nurbs_surface, outer_loops, inner_loops)
    brep_surface.SetId(1002)

    model_part.AddGeometry(brep_surface)
    return brep_surface


def CountGeometries(model_part):
    return sum(1 for _ in model_part.Geometries)

if __name__ == "__main__":
    KratosUnittest.main()