"""Client for NVIDIA PhysicsNeMo NIM microservices (DoMINO-Automotive-Aero).

A NIM is a served model behind HTTP: the solver host that calls one needs
neither weights, nor a GPU, nor torch - which is why this module is stdlib
urllib + numpy only. It implements the physics NIMs' documented contract
(nim/physicsnemo/domino-automotive-aero):

- ``POST /v1/infer`` (also ``/v1/infer/surface``, ``/v1/infer/volume``) as
  multipart/form-data: the geometry as an STL file field named
  ``design_stl`` plus string form fields (``stream_velocity``,
  ``stencil_size``, ``point_cloud_size``),
- the response body a numpy ``.npz`` archive (keys such as
  ``surface_coordinates``, ``pressure_surface``, ``wall_shear_stress``,
  ``coordinates``, ``velocity``, ``pressure``, ``drag_force``,
  ``lift_force``),
- ``GET /v1/health/ready`` for readiness, ``GET /v1/model/config`` for the
  served model's configuration.

Honesty note: this contract is implemented from NVIDIA's published NIM
documentation and pinned by payload tests against a stub transport -
running an actual NIM needs an NGC API key and docker, which the reference
environment does not have. ``IsReady()`` is the first thing to call against
a real deployment; a hosted gateway that requires auth is served by the
``api_key`` bearer header (a locally run container needs none - its NGC key
is docker-pull credentials, not an HTTP header).
"""

import io
import uuid

import numpy


def MakeBinaryStlBytes(points, triangles) -> bytes:
    """A binary STL of the given triangle surface, as bytes.

    The physics NIMs take their geometry as an STL upload; this builds one
    from the (N, 3) points / (T, 3) index triangles every mesh-bridge
    surface already provides. Normals are the triangles' own (right-hand
    winding); the attribute byte count is zero, per the format.
    """
    points = numpy.asarray(points, dtype=numpy.float64)
    triangles = numpy.asarray(triangles, dtype=numpy.int64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"Expected (N, 3) points, got shape {points.shape}.")
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError(f"Expected (T, 3) triangles, got shape {triangles.shape}.")

    a, b, c = (points[triangles[:, i]] for i in range(3))
    normals = numpy.cross(b - a, c - a)
    lengths = numpy.linalg.norm(normals, axis=1, keepdims=True)
    normals = numpy.divide(normals, lengths, out=numpy.zeros_like(normals),
                           where=lengths > 0.0)

    record = numpy.zeros(len(triangles), dtype=numpy.dtype([
        ("normal", "<f4", 3), ("vertices", "<f4", (3, 3)), ("attribute", "<u2")]))
    record["normal"] = normals
    record["vertices"] = numpy.stack([a, b, c], axis=1)

    header = b"PhysicsNeMoApplication binary STL".ljust(80, b" ")
    return header + numpy.uint32(len(triangles)).tobytes() + record.tobytes()


def EncodeMultipartFormData(fields, files) -> tuple:
    """(content_type, body) for a multipart/form-data request.

    Args:
        fields: {name: string value} plain form fields.
        files: {name: (filename, bytes)} file fields.
    """
    boundary = uuid.uuid4().hex
    body = io.BytesIO()
    for name, value in fields.items():
        body.write(f"--{boundary}\r\n"
                   f"Content-Disposition: form-data; name=\"{name}\"\r\n\r\n"
                   f"{value}\r\n".encode())
    for name, (filename, content) in files.items():
        body.write(f"--{boundary}\r\n"
                   f"Content-Disposition: form-data; name=\"{name}\"; "
                   f"filename=\"{filename}\"\r\n"
                   f"Content-Type: application/octet-stream\r\n\r\n".encode())
        body.write(content)
        body.write(b"\r\n")
    body.write(f"--{boundary}--\r\n".encode())
    return f"multipart/form-data; boundary={boundary}", body.getvalue()


def _UrllibTransport(url, method, headers, body, timeout):
    """(status, response bytes) through stdlib urllib; the default transport."""
    import urllib.error
    import urllib.request

    request = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status, response.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()
    except urllib.error.URLError as e:
        raise ConnectionError(
            f"Could not reach the NIM at {url} ({e.reason}). Is the container "
            "running? e.g. 'docker run --rm --runtime=nvidia --gpus 1 -p 8000:8000 "
            "-e NGC_API_KEY -t nvcr.io/nim/nvidia/domino-automotive-aero:2.1.0-41313772'"
        ) from e


class NimClient:
    """The HTTP client, with a swappable transport for offline testing.

    Args:
        base_url: e.g. "http://localhost:8000".
        timeout: Per-request timeout in seconds. Inference on a full
            vehicle takes tens of seconds; the documented example uses 120.
        api_key: Optional bearer token, for HOSTED deployments behind an
            authenticating gateway. A locally run NIM container ignores it.
    """

    def __init__(self, base_url: str, timeout: float = 120.0, api_key: str = ""):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.api_key = api_key
        self._transport = _UrllibTransport

    def SetTransport(self, transport) -> None:
        """Replaces the HTTP transport: transport(url, method, headers, body,
        timeout) -> (status, bytes). The test seam, like
        TritonInferenceProcess.SetClient."""
        self._transport = transport

    def _Headers(self, extra=None):
        headers = {"accept": "application/json"}
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"
        headers.update(extra or {})
        return headers

    def _Raise(self, url, status, body):
        snippet = body[:300].decode(errors="replace")
        if status == 401:
            raise RuntimeError(
                f"The NIM at {url} rejected the request as unauthorized (401). A "
                "hosted deployment needs its API key passed as \"api_key\"; a "
                "local container should not require one.")
        raise RuntimeError(f"NIM request to {url} failed with HTTP {status}: {snippet}")

    def IsReady(self) -> bool:
        """True when GET /v1/health/ready answers ready; False when the
        service is unreachable or still starting (models load on boot)."""
        url = f"{self.base_url}/v1/health/ready"
        try:
            status, _ = self._transport(url, "GET", self._Headers(), None, self.timeout)
        except ConnectionError:
            return False
        return status == 200

    def GetModelConfig(self) -> str:
        """The served model's configuration (YAML text) from /v1/model/config."""
        url = f"{self.base_url}/v1/model/config"
        status, body = self._transport(url, "GET", self._Headers(), None, self.timeout)
        if status != 200:
            self._Raise(url, status, body)
        return body.decode()

    def Infer(self, design_stl: bytes, form_fields: dict, endpoint: str = "infer") -> dict:
        """One inference call; the response .npz unpacked to numpy arrays.

        Args:
            design_stl: Binary STL bytes (MakeBinaryStlBytes).
            form_fields: The documented string form fields, e.g.
                {"stream_velocity": "30.0", "stencil_size": "1",
                 "point_cloud_size": "500000"}.
            endpoint: "infer" (volume + surface), "infer/surface" or
                "infer/volume".

        Returns:
            {key: numpy array} with the NIM's documented output keys.
        """
        if endpoint not in ("infer", "infer/surface", "infer/volume"):
            raise ValueError(
                f"Unsupported endpoint \"{endpoint}\". Use \"infer\", "
                "\"infer/surface\" or \"infer/volume\".")
        url = f"{self.base_url}/v1/{endpoint}"
        content_type, body = EncodeMultipartFormData(
            {name: str(value) for name, value in form_fields.items()},
            {"design_stl": ("design.stl", design_stl)})
        status, response = self._transport(
            url, "POST", self._Headers({"Content-Type": content_type}), body, self.timeout)
        if status != 200:
            self._Raise(url, status, response)
        try:
            with numpy.load(io.BytesIO(response)) as archive:
                return {key: numpy.array(archive[key]) for key in archive.files}
        except Exception as e:
            raise RuntimeError(
                f"The NIM at {url} answered HTTP 200 but not with a numpy .npz "
                f"archive ({e}); the first bytes were {response[:40]!r}.") from e
