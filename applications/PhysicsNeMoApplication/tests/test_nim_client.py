import io
import os
import struct

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.deployment import nim_client

# The documented output keys of the DoMINO-Automotive-Aero NIM, reproduced by
# the stub so the tests pin the real response layout.
_SURFACE_KEYS = ("surface_coordinates", "pressure_surface", "wall_shear_stress")


class _StubTransport:
    """Records the request and answers with a canned (status, body)."""

    def __init__(self, status=200, body=b""):
        self.status = status
        self.body = body
        self.calls = []

    def __call__(self, url, method, headers, body, timeout):
        self.calls.append({"url": url, "method": method, "headers": headers,
                           "body": body, "timeout": timeout})
        return self.status, self.body


def _NpzBytes(**arrays) -> bytes:
    buffer = io.BytesIO()
    numpy.savez(buffer, **arrays)
    return buffer.getvalue()


class TestNimClient(KratosUnittest.TestCase):
    """The client against a stub transport - the payload contract, pinned."""

    @staticmethod
    def _Tetrahedron():
        points = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                              [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        triangles = numpy.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
        return points, triangles

    def test_BinaryStlLayout(self):
        points, triangles = self._Tetrahedron()
        stl = nim_client.MakeBinaryStlBytes(points, triangles)
        # 80-byte header + uint32 count + 50 bytes per triangle
        self.assertEqual(len(stl), 80 + 4 + 50 * len(triangles))
        self.assertEqual(struct.unpack("<I", stl[80:84])[0], len(triangles))
        # first triangle's vertices, read back from the record layout
        vertex = struct.unpack("<3f", stl[84 + 12:84 + 24])
        self.assertVectorAlmostEqual(list(vertex), list(points[triangles[0][0]]), places=6)
        # the outward-wound facet's unit normal
        normal = struct.unpack("<3f", stl[84:84 + 12])
        self.assertVectorAlmostEqual(list(normal), [0.0, 0.0, -1.0], places=6)

    def test_InferPostsTheDocumentedMultipart(self):
        points, triangles = self._Tetrahedron()
        response = _NpzBytes(pressure_surface=numpy.zeros(3))
        transport = _StubTransport(body=response)
        client = nim_client.NimClient("http://localhost:8000/")
        client.SetTransport(transport)

        stl = nim_client.MakeBinaryStlBytes(points, triangles)
        result = client.Infer(stl, {"stream_velocity": "30.0", "stencil_size": "1",
                                    "point_cloud_size": "500000"})
        self.assertEqual(list(result), ["pressure_surface"])

        call = transport.calls[0]
        self.assertEqual(call["url"], "http://localhost:8000/v1/infer")
        self.assertEqual(call["method"], "POST")
        self.assertIn("multipart/form-data", call["headers"]["Content-Type"])
        body = call["body"]
        # the documented form fields and the design_stl file part, verbatim
        self.assertIn(b'name="stream_velocity"\r\n\r\n30.0', body)
        self.assertIn(b'name="stencil_size"\r\n\r\n1', body)
        self.assertIn(b'name="point_cloud_size"\r\n\r\n500000', body)
        self.assertIn(b'name="design_stl"; filename="design.stl"', body)
        self.assertIn(stl, body)
        # no auth header for a local container
        self.assertNotIn("Authorization", call["headers"])

    def test_ApiKeyBecomesBearerHeader(self):
        transport = _StubTransport(body=_NpzBytes(x=numpy.zeros(1)))
        client = nim_client.NimClient("http://hosted.example", api_key="nvapi-SECRET")
        client.SetTransport(transport)
        client.Infer(b"stl", {})
        self.assertEqual(transport.calls[0]["headers"]["Authorization"], "Bearer nvapi-SECRET")

    def test_EndpointsAndReadiness(self):
        transport = _StubTransport(body=_NpzBytes(x=numpy.zeros(1)))
        client = nim_client.NimClient("http://localhost:8000")
        client.SetTransport(transport)
        client.Infer(b"stl", {}, endpoint="infer/surface")
        self.assertEqual(transport.calls[-1]["url"], "http://localhost:8000/v1/infer/surface")
        with self.assertRaisesRegex(ValueError, "endpoint"):
            client.Infer(b"stl", {}, endpoint="predict")

        self.assertTrue(client.IsReady())
        self.assertEqual(transport.calls[-1]["url"], "http://localhost:8000/v1/health/ready")
        transport.status = 503
        self.assertFalse(client.IsReady())

    def test_ErrorsAreActionable(self):
        client = nim_client.NimClient("http://localhost:8000")
        client.SetTransport(_StubTransport(status=401, body=b"unauthorized"))
        with self.assertRaisesRegex(RuntimeError, "api_key"):
            client.Infer(b"stl", {})
        client.SetTransport(_StubTransport(status=500, body=b"boom"))
        with self.assertRaisesRegex(RuntimeError, "HTTP 500"):
            client.Infer(b"stl", {})
        client.SetTransport(_StubTransport(status=200, body=b"not an npz"))
        with self.assertRaisesRegex(RuntimeError, "npz"):
            client.Infer(b"stl", {})


class TestNimInferenceProcess(KratosUnittest.TestCase):
    """The in-loop process against the stub: gather -> POST -> nearest scatter."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        properties = self.model_part.CreateNewProperties(1)
        corners = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
        for i, (x, y, z) in enumerate(corners):
            self.model_part.CreateNewNode(i + 1, x, y, z)
        for i, triangle in enumerate([[1, 3, 2], [1, 2, 4], [1, 4, 3], [2, 3, 4]]):
            self.model_part.CreateNewCondition(
                "SurfaceCondition3D3N", i + 1, triangle, properties)

    def _Process(self, extra=""):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import nim_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"       : "Main",
                "stream_velocity"       : 25.0,
                "surface_output_fields" : [
                    { "nim_key" : "pressure_surface",  "variable_name" : "PRESSURE" },
                    { "nim_key" : "wall_shear_stress", "variable_name" : "VELOCITY" }
                ]%s
            }
        }""" % extra)
        return nim_inference_process.Factory(settings, self.model)

    def _StubResponse(self):
        # response points AT the node coordinates, so the nearest-neighbour
        # map is the identity and the assertion is exact
        coordinates = numpy.array(
            [[node.X, node.Y, node.Z] for node in self.model_part.Nodes])
        return _NpzBytes(
            surface_coordinates=coordinates,
            pressure_surface=1.0 + coordinates[:, 0] + 2.0 * coordinates[:, 1],
            wall_shear_stress=3.0 * coordinates,
            drag_force=numpy.array(1.25),
            lift_force=numpy.array(-0.5))

    def test_SurfaceFieldsAreMappedAndScalarsCollected(self):
        process = self._Process()
        transport = _StubTransport(body=self._StubResponse())
        process.GetClient().SetTransport(transport)

        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        for node in self.model_part.Nodes:
            self.assertAlmostEqual(
                node.GetSolutionStepValue(Kratos.PRESSURE), 1.0 + node.X + 2.0 * node.Y)
            self.assertVectorAlmostEqual(
                node.GetSolutionStepValue(Kratos.VELOCITY),
                [3.0 * node.X, 3.0 * node.Y, 3.0 * node.Z])
        self.assertAlmostEqual(process.last_scalars["drag_force"], 1.25)
        self.assertAlmostEqual(process.last_scalars["lift_force"], -0.5)

        # the surface that went out: 4 triangles as a binary STL (80-byte
        # header + count + 4 x 50-byte records), with the documented fields
        body = transport.calls[0]["body"]
        self.assertIn(b'name="stream_velocity"\r\n\r\n25.0', body)
        self.assertIn(b'name="design_stl"', body)
        header_at = body.find(b"PhysicsNeMoApplication binary STL")
        self.assertGreater(header_at, -1)
        self.assertEqual(struct.unpack("<I", body[header_at + 80:header_at + 84])[0], 4)
        self.assertEqual(transport.calls[0]["url"], "http://localhost:8000/v1/infer")

    def test_IntervalGating(self):
        process = self._Process(extra=""",
                "output_interval" : 2""")
        transport = _StubTransport(body=self._StubResponse())
        process.GetClient().SetTransport(transport)
        for step in (1, 2, 3):
            self.model_part.ProcessInfo[Kratos.STEP] = step
            process.ExecuteFinalizeSolutionStep()
        self.assertEqual(len(transport.calls), 1)

    def test_MissingResponseKeyIsActionable(self):
        process = self._Process()
        transport = _StubTransport(body=_NpzBytes(coordinates=numpy.zeros((1, 3))))
        process.GetClient().SetTransport(transport)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, "surface_coordinates"):
            process.ExecuteFinalizeSolutionStep()

    def test_NonTriangleSurfaceIsRejected(self):
        part = self.model.CreateModelPart("Quads")
        part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        properties = part.CreateNewProperties(1)
        for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)]):
            part.CreateNewNode(i + 1, float(x), float(y), 0.0)
        part.CreateNewCondition("SurfaceCondition3D4N", 1, [1, 2, 3, 4], properties)
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import nim_inference_process
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name"       : "Quads",
                "surface_output_fields" : [
                    { "nim_key" : "pressure_surface", "variable_name" : "PRESSURE" }
                ]
            }
        }""")
        process = nim_inference_process.Factory(settings, self.model)
        transport = _StubTransport(body=self._StubResponse())
        process.GetClient().SetTransport(transport)
        part.ProcessInfo[Kratos.STEP] = 1
        with self.assertRaisesRegex(ValueError, "TRIANGLE"):
            process.ExecuteFinalizeSolutionStep()


@KratosUnittest.skipUnless(os.environ.get("PHYSICSNEMO_NIM_URL"),
                           "Set PHYSICSNEMO_NIM_URL to a running NIM to verify live.")
class TestNimLive(KratosUnittest.TestCase):
    """Live verification against a real NIM - self-skips without one.

    Running a NIM needs an NGC API key and docker
    (nvcr.io/nim/nvidia/domino-automotive-aero); this environment has
    neither, which is exactly why the payload contract above is pinned
    against the stub instead.
    """

    def test_ReadyAndConfig(self):
        client = nim_client.NimClient(
            os.environ["PHYSICSNEMO_NIM_URL"],
            api_key=os.environ.get("NGC_API_KEY", ""))
        self.assertTrue(client.IsReady())
        self.assertIn("model", client.GetModelConfig().lower())


if __name__ == '__main__':
    KratosUnittest.main()
