# Import Kratos Multiphysics and the unittest framework
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import NumPy and standard modules
import numpy as np
import glob
import importlib.util
import os
import tempfile
import xml.etree.ElementTree as ET

# Attempt to import PyVista and the utilities module, set a flag if missing
try:
    import pyvista as pv
    import KratosMultiphysics.pyvista_utilities as pv_utils
    missing_pyvista = False
except ImportError:
    missing_pyvista = True

# The animation output process is a separate module; it is only importable once the
# build has been re-configured, since python_scripts is installed by a glob.
try:
    import KratosMultiphysics.pyvista_animation_output_process as pv_animation_process
    missing_animation_process = False
except ImportError:
    missing_animation_process = True

# GIF encoding needs the optional 'imageio' package, movies also 'imageio-ffmpeg'
missing_imageio = importlib.util.find_spec("imageio") is None
missing_imageio_ffmpeg = missing_imageio or importlib.util.find_spec("imageio_ffmpeg") is None


@KratosUnittest.skipIf(missing_pyvista, "Missing python libraries (pyvista)")
class TestPyVistaUtilities(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()

    # ------------------------------------------------------------------
    # ModelPartToPyVista — basic conversion
    # ------------------------------------------------------------------

    def test_convert_2d_mesh(self):
        mp = self.model.CreateModelPart("Main2D")
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.PRESSURE)

        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        mp.CreateNewNode(4, 1.0, 1.0, 0.0)

        for node in mp.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 10.0 * node.Id)
            node.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.1 * node.Id, 0.2 * node.Id, 0.0])

        # Shift current coordinates to verify undeformed vs deformed selection
        for node in mp.Nodes:
            node.X = node.X0 + 0.5
            node.Y = node.Y0 + 0.5

        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewElement("Element2D3N", 2, [2, 4, 3], prop)

        # Undeformed configuration
        grid = pv_utils.ModelPartToPyVista(
            mp, useDeformedConfiguration=False,
            nodalVariables=[KM.PRESSURE, KM.DISPLACEMENT]
        )

        self.assertIsInstance(grid, pv.UnstructuredGrid)
        self.assertEqual(grid.n_points, 4)
        self.assertEqual(grid.n_cells, 2)

        expected_points = np.array([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]
        ])
        self.assertTrue(np.allclose(grid.points, expected_points))

        self.assertIn("PRESSURE", grid.point_data)
        self.assertIn("DISPLACEMENT", grid.point_data)

        expected_pressure     = np.array([10.0, 20.0, 30.0, 40.0])
        expected_displacement = np.array([
            [0.1, 0.2, 0.0], [0.2, 0.4, 0.0],
            [0.3, 0.6, 0.0], [0.4, 0.8, 0.0]
        ])
        self.assertTrue(np.allclose(grid.point_data["PRESSURE"],     expected_pressure))
        self.assertTrue(np.allclose(grid.point_data["DISPLACEMENT"], expected_displacement))

        # Deformed configuration
        grid_deformed = pv_utils.ModelPartToPyVista(
            mp, useDeformedConfiguration=True, nodalVariables=[KM.PRESSURE]
        )
        expected_deformed = np.array([
            [0.5, 0.5, 0.0], [1.5, 0.5, 0.0],
            [0.5, 1.5, 0.0], [1.5, 1.5, 0.0]
        ])
        self.assertTrue(np.allclose(grid_deformed.points, expected_deformed))

    def test_convert_3d_mesh(self):
        mp = self.model.CreateModelPart("Main3D")
        mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)

        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        mp.CreateNewNode(4, 0.0, 0.0, 1.0)

        prop = mp.GetProperties()[0]
        elem = mp.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], prop)
        elem.SetValue(KM.DENSITY, 7800.0)

        grid = pv_utils.ModelPartToPyVista(mp, elementVariables=[KM.DENSITY])

        self.assertIsInstance(grid, pv.UnstructuredGrid)
        self.assertEqual(grid.n_points, 4)
        self.assertEqual(grid.n_cells, 1)

        self.assertIn("DENSITY", grid.cell_data)
        self.assertAlmostEqual(grid.cell_data["DENSITY"][0], 7800.0)

    def test_multiblock_conversion(self):
        mp = self.model.CreateModelPart("MultiBlockPart")
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)

        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N",   1, [1, 2, 3], prop)
        mp.CreateNewCondition("LineCondition2D2N", 1, [1, 2], prop)

        blocks = pv_utils.ModelPartToPyVista(
            mp, exportElements=True, exportConditions=True
        )

        self.assertIsInstance(blocks, pv.MultiBlock)
        self.assertIn("elements",   blocks.keys())
        self.assertIn("conditions", blocks.keys())

        elem_grid = blocks["elements"]
        cond_grid = blocks["conditions"]
        self.assertEqual(elem_grid.n_cells, 1)
        self.assertEqual(cond_grid.n_cells, 1)
        self.assertEqual(elem_grid.celltypes[0], 5)   # VTK_TRIANGLE
        self.assertEqual(cond_grid.celltypes[0], 3)   # VTK_LINE

    # ------------------------------------------------------------------
    # Matrix variable support
    # ------------------------------------------------------------------

    def test_matrix_variable_on_element(self):
        """Matrix-type variables (e.g. KM.LOCAL_AXES_MATRIX) are flattened row-major."""
        mp = self.model.CreateModelPart("MatrixPart")
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        prop = mp.GetProperties()[0]
        elem = mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        # Fill a 3x3 identity matrix into LOCAL_AXES_MATRIX (defined in core variables.h)
        mat = KM.Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                mat[i, j] = 1.0 if i == j else 0.0
        elem.SetValue(KM.LOCAL_AXES_MATRIX, mat)

        grid = pv_utils.ModelPartToPyVista(
            mp, elementVariables=[KM.LOCAL_AXES_MATRIX]
        )

        self.assertIn("LOCAL_AXES_MATRIX", grid.cell_data)
        stored = grid.cell_data["LOCAL_AXES_MATRIX"]
        self.assertEqual(stored.shape, (1, 9))       # 1 cell × 3×3 flat

        # Diagonal elements should be 1 (row-major: [0]=0,0; [4]=1,1; [8]=2,2)
        self.assertAlmostEqual(stored[0, 0], 1.0)
        self.assertAlmostEqual(stored[0, 4], 1.0)
        self.assertAlmostEqual(stored[0, 8], 1.0)
        self.assertAlmostEqual(stored[0, 1], 0.0)

    def test_matrix_variable_on_node(self):
        """Non-historical matrix variables on nodes are stored as flat row-major arrays."""
        mp = self.model.CreateModelPart("MatrixNodePart")
        n1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        n3 = mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        mat = KM.Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                mat[i, j] = float(i * 3 + j)
        for node in mp.Nodes:
            node.SetValue(KM.LOCAL_AXES_MATRIX, mat)

        grid = pv_utils.ModelPartToPyVista(
            mp, nodalVariables=[KM.LOCAL_AXES_MATRIX]
        )

        self.assertIn("LOCAL_AXES_MATRIX", grid.point_data)
        stored = grid.point_data["LOCAL_AXES_MATRIX"]
        self.assertEqual(stored.shape, (3, 9))
        self.assertAlmostEqual(stored[0, 0], 0.0)
        self.assertAlmostEqual(stored[0, 4], 4.0)
        self.assertAlmostEqual(stored[0, 8], 8.0)

    # ------------------------------------------------------------------
    # File I/O
    # ------------------------------------------------------------------

    def test_SaveModelPart(self):
        mp = self.model.CreateModelPart("SavePart")
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = os.path.join(tmp_dir, "mesh.vtu")
            pv_utils.SaveModelPart(mp, file_path)
            self.assertTrue(os.path.exists(file_path))

    # ------------------------------------------------------------------
    # PlotModelPart
    # ------------------------------------------------------------------

    def _make_2d_mp(self, name):
        mp = self.model.CreateModelPart(name)
        mp.AddNodalSolutionStepVariable(KM.PRESSURE)
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        for i, (x, y) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)], 1):
            n = mp.CreateNewNode(i, float(x), float(y), 0.0)
            n.SetSolutionStepValue(KM.PRESSURE, 0, float(i) * 10.0)
            n.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.1 * x, 0.1 * y, 0.0])
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewElement("Element2D3N", 2, [2, 4, 3], prop)
        return mp

    def test_PlotModelPart_scalar(self):
        """PlotModelPart returns a Plotter and colours by a scalar variable."""
        mp = self._make_2d_mp("Plot2DScalar")
        plotter = pv_utils.PlotModelPart(
            mp, name=KM.PRESSURE, offScreen=True
        )
        self.assertIsInstance(plotter, pv.Plotter)
        plotter.close()

    def test_PlotModelPart_vector_component(self):
        """Component selection extracts a single column from a vector variable."""
        mp = self._make_2d_mp("Plot2DVecComp")
        for comp in [0, 1, 2, None]:
            plotter = pv_utils.PlotModelPart(
                mp, name=KM.DISPLACEMENT, component=comp, offScreen=True
            )
            self.assertIsInstance(plotter, pv.Plotter)
            plotter.close()

    def test_PlotModelPart_warp_by_vector(self):
        """warpByVector warps the mesh and optionally shows the undeformed ghost."""
        mp = self._make_2d_mp("PlotWarp")

        # With undeformed ghost
        plotter = pv_utils.PlotModelPart(
            mp, name=KM.PRESSURE, warpByVector=KM.DISPLACEMENT,
            factor=2.0, showUndeformed=True, offScreen=True
        )
        self.assertIsInstance(plotter, pv.Plotter)
        plotter.close()

        # Without ghost
        plotter = pv_utils.PlotModelPart(
            mp, name=KM.PRESSURE, warpByVector=KM.DISPLACEMENT,
            factor=2.0, showUndeformed=False, offScreen=True
        )
        self.assertIsInstance(plotter, pv.Plotter)
        plotter.close()

    def test_PlotModelPart_custom_label(self):
        """Custom label is accepted without error."""
        mp = self._make_2d_mp("PlotCustomLabel")
        plotter = pv_utils.PlotModelPart(
            mp, name=KM.PRESSURE, label="My Pressure [Pa]", offScreen=True
        )
        self.assertIsInstance(plotter, pv.Plotter)
        plotter.close()

    def test_PlotModelPart_2d_auto_camera(self):
        """Flat 2-D meshes trigger xy camera with parallel projection."""
        mp = self._make_2d_mp("PlotCam2D")
        plotter = pv_utils.PlotModelPart(mp, offScreen=True, view="default")
        self.assertIsInstance(plotter, pv.Plotter)
        plotter.close()

    def test_PlotModelPart_explicit_view(self):
        """An explicit view string is accepted without error."""
        mp = self._make_2d_mp("PlotExplicitView")
        plotter = pv_utils.PlotModelPart(mp, offScreen=True, view="xy")
        self.assertIsInstance(plotter, pv.Plotter)
        plotter.close()

    def test_PlotModelPart_composed_plotter(self):
        """Passing an existing plotter accumulates actors on it."""
        mp = self._make_2d_mp("PlotCompose")
        base_plotter = pv.Plotter(off_screen=True)
        returned = pv_utils.PlotModelPart(
            mp, name=KM.PRESSURE, plotter=base_plotter, offScreen=True
        )
        self.assertIs(returned, base_plotter)
        base_plotter.close()

    def test_PlotModelPart_no_variable(self):
        """PlotModelPart works with name=None (topology only)."""
        mp = self._make_2d_mp("PlotNoVar")
        plotter = pv_utils.PlotModelPart(mp, offScreen=True)
        self.assertIsInstance(plotter, pv.Plotter)
        plotter.close()

    def test_PlotModelPart_theme(self):
        """A named theme is applied without error."""
        mp = self._make_2d_mp("PlotTheme")
        plotter = pv_utils.PlotModelPart(
            mp, name=KM.PRESSURE, theme="document", offScreen=True
        )
        self.assertIsInstance(plotter, pv.Plotter)
        plotter.close()

    # ------------------------------------------------------------------
    # ScreenshotModelPart
    # ------------------------------------------------------------------

    def test_ScreenshotModelPart_to_file(self):
        """ScreenshotModelPart writes a PNG file."""
        mp = self._make_2d_mp("Screenshot2D")
        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = os.path.join(tmp_dir, "screenshot.png")
            pv_utils.ScreenshotModelPart(mp, file_path, name=KM.PRESSURE)
            self.assertTrue(os.path.exists(file_path))
            self.assertGreater(os.path.getsize(file_path), 0)

    def test_ScreenshotModelPart_return_array(self):
        """Passing filename=None returns a numpy image array."""
        mp = self._make_2d_mp("ScreenshotArr")
        img = pv_utils.ScreenshotModelPart(mp, filename=None, name=KM.PRESSURE)
        self.assertIsInstance(img, np.ndarray)
        self.assertEqual(img.ndim, 3)   # height × width × channels

    # ------------------------------------------------------------------
    # CreateExtractedSurface
    # ------------------------------------------------------------------

    def test_CreateExtractedSurface(self):
        """CreateExtractedSurface returns a PolyData surface from a 3-D mesh."""
        mp = self.model.CreateModelPart("SurfacePart")
        mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        for i, coords in enumerate([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]
        ], 1):
            n = mp.CreateNewNode(i, *coords)
            n.SetSolutionStepValue(KM.TEMPERATURE, 0, float(i))
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element3D8N", 1, list(range(1, 9)), prop)

        surface = pv_utils.CreateExtractedSurface(
            mp, nodalVariables=[KM.TEMPERATURE]
        )
        self.assertIsInstance(surface, pv.PolyData)
        self.assertGreater(surface.n_points, 0)
        self.assertGreater(surface.n_cells,  0)
        self.assertIn("TEMPERATURE", surface.point_data)

    def test_CreateExtractedSurface_empty(self):
        """CreateExtractedSurface on an empty ModelPart returns empty PolyData."""
        mp = self.model.CreateModelPart("EmptySurfacePart")
        surface = pv_utils.CreateExtractedSurface(mp)
        self.assertIsInstance(surface, pv.PolyData)

    # ------------------------------------------------------------------
    # Helper label/component utilities
    # ------------------------------------------------------------------

    def test_BuildScalarLabel_scalar(self):
        label = pv_utils._BuildScalarLabel("PRESSURE", 0, 1)
        self.assertEqual(label, "PRESSURE")

    def test_BuildScalarLabel_vector(self):
        self.assertEqual(pv_utils._BuildScalarLabel("VELOCITY", 0, 3), "VELOCITY X")
        self.assertEqual(pv_utils._BuildScalarLabel("VELOCITY", 1, 3), "VELOCITY Y")
        self.assertEqual(pv_utils._BuildScalarLabel("VELOCITY", 2, 3), "VELOCITY Z")
        self.assertEqual(pv_utils._BuildScalarLabel("VELOCITY", None, 3), "VELOCITY Magnitude")

    def test_BuildScalarLabel_voigt(self):
        self.assertEqual(pv_utils._BuildScalarLabel("STRESS", 0, 6), "STRESS XX")
        self.assertEqual(pv_utils._BuildScalarLabel("STRESS", 3, 6), "STRESS XY")

    def test_BuildScalarLabel_principal(self):
        label = pv_utils._BuildScalarLabel("Principal Values of STRAIN", 0, 3)
        self.assertIn("Max.", label)
        label = pv_utils._BuildScalarLabel("Principal Values of STRAIN", 2, 3)
        self.assertIn("Min.", label)

    def test_ExtractComponent_scalar(self):
        arr = np.array([1.0, 2.0, 3.0])
        out, ndim = pv_utils._ExtractComponent(arr, 0, "PRESSURE")
        self.assertEqual(ndim, 1)
        self.assertTrue(np.allclose(out, arr))

    def test_ExtractComponent_vector_component(self):
        arr = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        out, ndim = pv_utils._ExtractComponent(arr, 1, "VELOCITY")
        self.assertEqual(ndim, 3)
        self.assertTrue(np.allclose(out, [2.0, 5.0]))

    def test_ExtractComponent_magnitude(self):
        arr = np.array([[3.0, 4.0, 0.0]])
        out, ndim = pv_utils._ExtractComponent(arr, None, "VELOCITY")
        self.assertAlmostEqual(out[0], 5.0)

    def test_ExtractComponent_principal_sorting(self):
        """Principal-value arrays must be sorted in descending order."""
        arr = np.array([[1.0, 3.0, 2.0], [5.0, 4.0, 6.0]])
        out, ndim = pv_utils._ExtractComponent(arr, 0, "Principal Values of STRAIN")
        # Max principal should be the largest value per row
        self.assertAlmostEqual(out[0], 3.0)
        self.assertAlmostEqual(out[1], 6.0)

    # ------------------------------------------------------------------
    # Existing post-processing utilities
    # ------------------------------------------------------------------

    def test_ComputeWarpedMesh(self):
        mp = self.model.CreateModelPart("WarpPart")
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        n1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        n3 = mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        for n in [n1, n2, n3]:
            n.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.1, 0.2, 0.0])

        warped = pv_utils.ComputeWarpedMesh(mp, KM.DISPLACEMENT, factor=2.0)
        self.assertIsInstance(warped, pv.UnstructuredGrid)
        self.assertEqual(warped.n_points, 3)
        expected = np.array([[0.2, 0.4, 0.0], [1.2, 0.4, 0.0], [0.2, 1.4, 0.0]])
        self.assertTrue(np.allclose(warped.points, expected))

    def test_CreateOrthogonalSlices(self):
        mp = self.model.CreateModelPart("SlicePart")
        for i, coords in enumerate([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]
        ], 1):
            mp.CreateNewNode(i, *coords)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element3D8N", 1, list(range(1, 9)), prop)

        slices = pv_utils.CreateOrthogonalSlices(mp, x=0.5, y=0.5, z=0.5)
        self.assertIsInstance(slices, pv.MultiBlock)
        self.assertGreater(len(slices), 0)

    def test_CreateIsosurfaces(self):
        mp = self.model.CreateModelPart("ContourPart")
        mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        for i, coords in enumerate([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]
        ], 1):
            n = mp.CreateNewNode(i, *coords)
            n.SetSolutionStepValue(KM.TEMPERATURE, 0, float(i))
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element3D8N", 1, list(range(1, 9)), prop)

        contours = pv_utils.CreateIsosurfaces(mp, KM.TEMPERATURE, valuesOrNumber=3)
        self.assertIsInstance(contours, pv.PolyData)

    def test_CreateStreamlines(self):
        mp = self.model.CreateModelPart("StreamlinePart")
        mp.AddNodalSolutionStepVariable(KM.VELOCITY)
        for i, coords in enumerate([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]
        ], 1):
            n = mp.CreateNewNode(i, *coords)
            n.SetSolutionStepValue(KM.VELOCITY, 0, [1.0, 0.0, 0.0])
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element3D8N", 1, list(range(1, 9)), prop)

        streamlines = pv_utils.CreateStreamlines(
            mp, velocityVariable=KM.VELOCITY,
            sourceCenter=[0.1, 0.5, 0.5], sourceRadius=0.1, nPoints=5
        )
        self.assertIsInstance(streamlines, pv.PolyData)

    def test_CreateThresholdedMesh(self):
        mp = self.model.CreateModelPart("ThresholdPart")
        mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        for i, coords in enumerate([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]
        ], 1):
            n = mp.CreateNewNode(i, *coords)
            n.SetSolutionStepValue(KM.TEMPERATURE, 0, float(i))
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element3D8N", 1, list(range(1, 9)), prop)

        self.assertIsInstance(
            pv_utils.CreateThresholdedMesh(mp, KM.TEMPERATURE, 4.5, "above"),
            pv.UnstructuredGrid
        )
        self.assertIsInstance(
            pv_utils.CreateThresholdedMesh(mp, KM.TEMPERATURE, 4.5, "below"),
            pv.UnstructuredGrid
        )
        self.assertIsInstance(
            pv_utils.CreateThresholdedMesh(mp, KM.TEMPERATURE, [2.5, 6.5], "between"),
            pv.UnstructuredGrid
        )

    def test_CreateThresholdedMesh_invalid_type(self):
        mp = self.model.CreateModelPart("ThresholdInvalid")
        mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        for i, coords in enumerate([[0, 0, 0], [1, 0, 0], [0, 1, 0]], 1):
            n = mp.CreateNewNode(i, *coords)
            n.SetSolutionStepValue(KM.TEMPERATURE, 0, float(i))
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        with self.assertRaises(ValueError):
            pv_utils.CreateThresholdedMesh(mp, KM.TEMPERATURE, 1.0, "unknown")

    def test_CreateVectorGlyphs(self):
        mp = self.model.CreateModelPart("GlyphPart")
        mp.AddNodalSolutionStepVariable(KM.VELOCITY)
        for i, coords in enumerate([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], 1):
            n = mp.CreateNewNode(i, *coords)
            n.SetSolutionStepValue(KM.VELOCITY, 0, [1.0, 2.0, 0.0])
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D4N", 1, [1, 2, 3, 4], prop)

        glyphs = pv_utils.CreateVectorGlyphs(mp, KM.VELOCITY, scaleFactor=0.5, glyphType="cone")
        self.assertIsInstance(glyphs, pv.PolyData)
        self.assertIn("VELOCITY", glyphs.point_data)
        self.assertIn("VELOCITY_magnitude", glyphs.point_data)

    def test_CreateClippedMesh(self):
        mp = self.model.CreateModelPart("ClipPart")
        for i, coords in enumerate([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]
        ], 1):
            mp.CreateNewNode(i, *coords)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element3D8N", 1, list(range(1, 9)), prop)

        clipped = pv_utils.CreateClippedMesh(
            mp, normal=[1.0, 0.0, 0.0], origin=[0.5, 0.5, 0.5]
        )
        self.assertIsInstance(clipped, pv.UnstructuredGrid)

    # ------------------------------------------------------------------
    # Transient helpers
    # ------------------------------------------------------------------

    def _make_transient_mp(self, name):
        """Build the 2-D fixture with a solution-step buffer, ready to be advanced."""
        mp = self.model.CreateModelPart(name)
        mp.AddNodalSolutionStepVariable(KM.PRESSURE)
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.SetBufferSize(2)
        for i, (x, y) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)], 1):
            mp.CreateNewNode(i, float(x), float(y), 0.0)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewElement("Element2D3N", 2, [2, 4, 3], prop)
        mp.ProcessInfo[KM.STEP] = 0
        mp.ProcessInfo[KM.TIME] = 0.0
        return mp

    def _advance_transient_mp(self, mp, step, dt=0.1):
        """Advance *mp* one time step and refresh its nodal results."""
        mp.CloneTimeStep(dt * step)
        mp.ProcessInfo[KM.STEP] = step
        for node in mp.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, float(node.Id) * step)
            node.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.01 * step, 0.02 * step, 0.0])
        return mp

    # ------------------------------------------------------------------
    # Transient helpers - pure unit tests
    # ------------------------------------------------------------------

    def test_NaturalSortKey(self):
        """Digit runs sort numerically, so step 2 precedes step 10."""
        names = ["res_10.vtu", "res_2.vtu", "res_1.vtu", "res_21.vtu"]
        self.assertEqual(
            sorted(names, key=pv_utils._NaturalSortKey),
            ["res_1.vtu", "res_2.vtu", "res_10.vtu", "res_21.vtu"]
        )

    def test_PrettyTime(self):
        """Float noise from repeated CloneTimeStep is rounded away."""
        self.assertEqual(pv_utils._PrettyTime(0.30000000000000004), 0.3)
        self.assertEqual(pv_utils._PrettyTime(1), 1.0)

    def test_IsAnimationAvailable(self):
        """PNG needs nothing extra; the encoders are probed by import."""
        self.assertTrue(pv_utils.IsAnimationAvailable("png"))
        self.assertTrue(pv_utils.IsAnimationAvailable(".png"))
        self.assertTrue(pv_utils.IsAnimationAvailable("pvd"))
        self.assertIsInstance(pv_utils.IsAnimationAvailable("gif"), bool)
        self.assertIsInstance(pv_utils.IsAnimationAvailable(".mp4"), bool)

    # ------------------------------------------------------------------
    # TransientPlotter
    # ------------------------------------------------------------------

    def test_TransientPlotter_png_sequence(self):
        """A .png target writes one zero-padded frame per AddFrame call."""
        mp = self._make_transient_mp("TransientPng")
        with tempfile.TemporaryDirectory() as tmp_dir:
            target = os.path.join(tmp_dir, "anim.png")
            with pv_utils.TransientPlotter(target, name=KM.PRESSURE) as recorder:
                for step in range(1, 5):
                    self._advance_transient_mp(mp, step)
                    recorder.AddFrame(mp)
                self.assertEqual(recorder.numberOfFrames, 4)

            written = sorted(os.path.basename(f)
                             for f in glob.glob(os.path.join(tmp_dir, "anim_*.png")))
            self.assertEqual(
                written,
                ["anim_0000.png", "anim_0001.png", "anim_0002.png", "anim_0003.png"]
            )

    def test_TransientPlotter_directory_target(self):
        """A name without a media suffix becomes a directory of frame_XXXX.png."""
        mp = self._make_transient_mp("TransientDir")
        with tempfile.TemporaryDirectory() as tmp_dir:
            target = os.path.join(tmp_dir, "frames")
            with pv_utils.TransientPlotter(target, name=KM.PRESSURE) as recorder:
                for step in range(1, 3):
                    self._advance_transient_mp(mp, step)
                    recorder.AddFrame(mp)
            self.assertTrue(os.path.isdir(target))
            self.assertEqual(len(glob.glob(os.path.join(target, "frame_*.png"))), 2)

    def test_TransientPlotter_no_filename_returns_frames(self):
        """Without a filename the frames are kept in memory as image arrays."""
        mp = self._make_transient_mp("TransientMem")
        recorder = pv_utils.TransientPlotter(None, name=KM.PRESSURE)
        for step in range(1, 4):
            self._advance_transient_mp(mp, step)
            recorder.AddFrame(mp)
        self.assertEqual(recorder.Close(), 3)
        self.assertEqual(len(recorder.frames), 3)
        self.assertEqual(recorder.frames[0].ndim, 3)

    def test_TransientPlotter_reads_time_from_process_info(self):
        """The time annotation defaults to ProcessInfo[TIME]."""
        mp = self._make_transient_mp("TransientTime")
        self._advance_transient_mp(mp, 3)
        self.assertAlmostEqual(pv_utils._PrettyTime(mp.ProcessInfo[KM.TIME]), 0.3)

        recorder = pv_utils.TransientPlotter(None, name=KM.PRESSURE, timeAnnotation=True)
        recorder.AddFrame(mp)
        recorder.AddFrame(mp, time=1.25)
        self.assertEqual(recorder.Close(), 2)

    def test_TransientPlotter_topology_change(self):
        """Frames survive a topology change between steps (adaptive remeshing)."""
        mp = self.model.CreateModelPart("TransientRemesh")
        mp.AddNodalSolutionStepVariable(KM.PRESSURE)
        mp.SetBufferSize(2)
        for i, (x, y) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)], 1):
            mp.CreateNewNode(i, float(x), float(y), 0.0)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        recorder = pv_utils.TransientPlotter(None, name=KM.PRESSURE)
        for step in range(1, 5):
            mp.CloneTimeStep(0.1 * step)
            for node in mp.Nodes:
                node.SetSolutionStepValue(KM.PRESSURE, 0, float(node.Id) * step)
            if step == 2:
                # Cell count changes mid-animation
                mp.CreateNewElement("Element2D3N", 2, [2, 4, 3], prop)
            recorder.AddFrame(mp)

        self.assertEqual(recorder.Close(), 4)
        self.assertEqual(len(recorder.frames), 4)

    def test_TransientPlotter_explicit_clim(self):
        """An explicit color range is accepted and kept for every frame."""
        mp = self._make_transient_mp("TransientClim")
        recorder = pv_utils.TransientPlotter(None, name=KM.PRESSURE, clim=(0.0, 100.0))
        for step in range(1, 3):
            self._advance_transient_mp(mp, step)
            recorder.AddFrame(mp)
        self.assertEqual(recorder.clim, (0.0, 100.0))
        recorder.Close()

    def test_TransientPlotter_lock_clim(self):
        """lockClim freezes the range at the first frame instead of rescaling."""
        mp = self._make_transient_mp("TransientLockClim")
        recorder = pv_utils.TransientPlotter(None, name=KM.PRESSURE, lockClim=True)
        self._advance_transient_mp(mp, 1)
        recorder.AddFrame(mp)
        locked = recorder.clim
        self.assertIsNotNone(locked)

        self._advance_transient_mp(mp, 5)   # much larger values
        recorder.AddFrame(mp)
        self.assertEqual(recorder.clim, locked)
        recorder.Close()

    def test_TransientPlotter_warp_and_component(self):
        """The full PlotModelPart vocabulary is accepted by the recorder."""
        mp = self._make_transient_mp("TransientWarp")
        recorder = pv_utils.TransientPlotter(
            None, name=KM.DISPLACEMENT, component=1, warpByVector=KM.DISPLACEMENT,
            factor=2.0, showUndeformed=True, cmap="viridis", view="xy"
        )
        for step in range(1, 3):
            self._advance_transient_mp(mp, step)
            recorder.AddFrame(mp)
        self.assertEqual(recorder.Close(), 2)

    def test_TransientPlotter_grid_source(self):
        """A frame may be a PyVista grid instead of a ModelPart."""
        mp = self._make_transient_mp("TransientGridSrc")
        self._advance_transient_mp(mp, 1)
        grid = pv_utils.ModelPartToPyVista(mp, nodalVariables=[KM.PRESSURE])

        recorder = pv_utils.TransientPlotter(None, name=KM.PRESSURE)
        recorder.AddFrame(grid)
        recorder.AddFrame(grid, time=0.5)
        self.assertEqual(recorder.Close(), 2)

    def test_TransientPlotter_closed_rejects_frames(self):
        """Adding a frame after Close is an error rather than a silent no-op."""
        mp = self._make_transient_mp("TransientClosed")
        self._advance_transient_mp(mp, 1)
        recorder = pv_utils.TransientPlotter(None, name=KM.PRESSURE)
        recorder.AddFrame(mp)
        recorder.Close()
        with self.assertRaises(RuntimeError):
            recorder.AddFrame(mp)

    def test_TransientPlotter_plotter_property(self):
        """The underlying plotter is exposed for custom compositing."""
        mp = self._make_transient_mp("TransientPlotterProp")
        self._advance_transient_mp(mp, 1)
        recorder = pv_utils.TransientPlotter(None, name=KM.PRESSURE)
        self.assertIsInstance(recorder.plotter, pv.Plotter)
        recorder.AddFrame(mp)
        recorder.Close()

    @KratosUnittest.skipIf(missing_imageio, "Missing python libraries (imageio)")
    def test_TransientPlotter_gif(self):
        """A .gif target produces a readable animated GIF."""
        mp = self._make_transient_mp("TransientGif")
        with tempfile.TemporaryDirectory() as tmp_dir:
            target = os.path.join(tmp_dir, "anim.gif")
            with pv_utils.TransientPlotter(target, name=KM.PRESSURE, fps=8) as recorder:
                for step in range(1, 4):
                    self._advance_transient_mp(mp, step)
                    recorder.AddFrame(mp)
            self.assertTrue(os.path.exists(target))
            self.assertGreater(os.path.getsize(target), 0)

    @KratosUnittest.skipIf(missing_imageio_ffmpeg, "Missing python libraries (imageio-ffmpeg)")
    def test_TransientPlotter_movie(self):
        """An .mp4 target produces a movie file."""
        mp = self._make_transient_mp("TransientMovie")
        with tempfile.TemporaryDirectory() as tmp_dir:
            target = os.path.join(tmp_dir, "anim.mp4")
            with pv_utils.TransientPlotter(target, name=KM.PRESSURE, fps=8) as recorder:
                for step in range(1, 4):
                    self._advance_transient_mp(mp, step)
                    recorder.AddFrame(mp)
            self.assertTrue(os.path.exists(target))
            self.assertGreater(os.path.getsize(target), 0)

    @KratosUnittest.skipIf(not missing_imageio, "imageio is installed")
    def test_TransientPlotter_gif_without_imageio_raises(self):
        """Without imageio the GIF path fails with an actionable message."""
        mp = self._make_transient_mp("TransientNoImageio")
        self._advance_transient_mp(mp, 1)
        with tempfile.TemporaryDirectory() as tmp_dir:
            recorder = pv_utils.TransientPlotter(
                os.path.join(tmp_dir, "anim.gif"), name=KM.PRESSURE
            )
            with self.assertRaises(ImportError) as context:
                recorder.AddFrame(mp)
            self.assertIn("imageio", str(context.exception))
            self.assertIn(".png", str(context.exception))
            recorder.Close()

    # ------------------------------------------------------------------
    # Animate
    # ------------------------------------------------------------------

    def test_Animate_from_grids(self):
        """An in-memory list of grids animates without touching the disk first."""
        mp = self._make_transient_mp("AnimateGrids")
        grids = []
        for step in range(1, 4):
            self._advance_transient_mp(mp, step)
            grids.append(pv_utils.ModelPartToPyVista(mp, nodalVariables=[KM.PRESSURE]))

        with tempfile.TemporaryDirectory() as tmp_dir:
            target = os.path.join(tmp_dir, "grids.png")
            self.assertEqual(pv_utils.Animate(grids, target, name=KM.PRESSURE), 3)
            self.assertEqual(len(glob.glob(os.path.join(tmp_dir, "grids_*.png"))), 3)

    def test_Animate_from_model_parts(self):
        """A list of ModelPart snapshots is also a valid source."""
        mp = self._make_transient_mp("AnimateModelParts")
        self._advance_transient_mp(mp, 1)
        with tempfile.TemporaryDirectory() as tmp_dir:
            target = os.path.join(tmp_dir, "mps.png")
            self.assertEqual(pv_utils.Animate([mp, mp], target, name=KM.PRESSURE), 2)

    def test_Animate_from_files(self):
        """A glob source is read back in natural (not lexicographic) order."""
        mp = self._make_transient_mp("AnimateFiles")
        self._advance_transient_mp(mp, 1)
        grid = pv_utils.ModelPartToPyVista(mp, nodalVariables=[KM.PRESSURE])

        with tempfile.TemporaryDirectory() as tmp_dir:
            for index in range(12):
                grid.save(os.path.join(tmp_dir, f"res_{index}.vtu"))

            paths = sorted(glob.glob(os.path.join(tmp_dir, "res_*.vtu")),
                           key=pv_utils._NaturalSortKey)
            self.assertEqual(os.path.basename(paths[2]), "res_2.vtu")
            self.assertEqual(os.path.basename(paths[10]), "res_10.vtu")

            target = os.path.join(tmp_dir, "out", "series.png")
            self.assertEqual(
                pv_utils.Animate(os.path.join(tmp_dir, "res_*.vtu"), target, name=KM.PRESSURE),
                12
            )

    def test_Animate_no_match_raises(self):
        """An empty source is reported rather than producing an empty animation."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            with self.assertRaises(ValueError):
                pv_utils.Animate(os.path.join(tmp_dir, "nothing_*.vtu"),
                                 os.path.join(tmp_dir, "out.png"))
            with self.assertRaises(ValueError):
                pv_utils.Animate([], os.path.join(tmp_dir, "out.png"))

    # ------------------------------------------------------------------
    # PvdWriter
    # ------------------------------------------------------------------

    def test_PvdWriter(self):
        """A .vtu per step plus a .pvd collection indexing them by time."""
        mp = self._make_transient_mp("PvdPart")
        with tempfile.TemporaryDirectory() as tmp_dir:
            target = os.path.join(tmp_dir, "run.pvd")
            with pv_utils.PvdWriter(target, nodalVariables=[KM.PRESSURE]) as writer:
                for step in range(1, 5):
                    self._advance_transient_mp(mp, step)
                    writer.AddStep(mp)
                self.assertEqual(writer.numberOfSteps, 4)

            self.assertTrue(os.path.exists(target))
            self.assertEqual(len(glob.glob(os.path.join(tmp_dir, "run", "run_*.vtu"))), 4)

            tree = ET.parse(target)
            datasets = tree.getroot().find("Collection").findall("DataSet")
            self.assertEqual(len(datasets), 4)
            self.assertEqual(
                [float(d.get("timestep")) for d in datasets], [0.1, 0.2, 0.3, 0.4]
            )
            self.assertEqual(datasets[0].get("file"), "run/run_0000.vtu")

    def test_PvdWriter_appends_suffix(self):
        """A filename without the .pvd suffix gets one."""
        mp = self._make_transient_mp("PvdSuffix")
        self._advance_transient_mp(mp, 1)
        with tempfile.TemporaryDirectory() as tmp_dir:
            writer = pv_utils.PvdWriter(os.path.join(tmp_dir, "run"))
            writer.AddStep(mp)
            self.assertTrue(writer.Close().endswith("run.pvd"))

    def test_PvdWriter_closed_rejects_steps(self):
        """Adding a step after Close is an error."""
        mp = self._make_transient_mp("PvdClosed")
        self._advance_transient_mp(mp, 1)
        with tempfile.TemporaryDirectory() as tmp_dir:
            writer = pv_utils.PvdWriter(os.path.join(tmp_dir, "run.pvd"))
            writer.AddStep(mp)
            writer.Close()
            with self.assertRaises(RuntimeError):
                writer.AddStep(mp)

    def test_Animate_from_pvd(self):
        """A .pvd round trip: write the series, then animate it with its own times."""
        mp = self._make_transient_mp("PvdRoundTrip")
        with tempfile.TemporaryDirectory() as tmp_dir:
            collection = os.path.join(tmp_dir, "run.pvd")
            with pv_utils.PvdWriter(collection, nodalVariables=[KM.PRESSURE]) as writer:
                for step in range(1, 4):
                    self._advance_transient_mp(mp, step)
                    writer.AddStep(mp)

            reader = pv.PVDReader(collection)
            self.assertEqual(reader.number_time_points, 3)
            self.assertEqual(list(reader.time_values), [0.1, 0.2, 0.3])

            target = os.path.join(tmp_dir, "out", "run.png")
            self.assertEqual(pv_utils.Animate(collection, target, name=KM.PRESSURE), 3)
            self.assertEqual(len(glob.glob(os.path.join(tmp_dir, "out", "run_*.png"))), 3)


@KratosUnittest.skipIf(missing_pyvista, "Missing python libraries (pyvista)")
@KratosUnittest.skipIf(missing_animation_process,
                       "KratosMultiphysics.pyvista_animation_output_process is not installed")
class TestPyVistaAnimationOutputProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()

    def _make_model_part(self, name="Main"):
        mp = self.model.CreateModelPart(name)
        mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.SetBufferSize(2)
        for i, (x, y) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)], 1):
            mp.CreateNewNode(i, float(x), float(y), 0.0)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewElement("Element2D3N", 2, [2, 4, 3], prop)
        mp.ProcessInfo[KM.IS_RESTARTED] = False
        mp.ProcessInfo[KM.STEP] = 0
        mp.ProcessInfo[KM.TIME] = 0.0
        return mp

    def _run_steps(self, process, mp, number_of_steps):
        """Drive the process the way AnalysisStage.OutputSolutionStep does."""
        printed = 0
        for step in range(1, number_of_steps + 1):
            mp.CloneTimeStep(0.1 * step)
            mp.ProcessInfo[KM.STEP] = step
            for node in mp.Nodes:
                node.SetSolutionStepValue(KM.TEMPERATURE, 0, float(node.Id) * step)
                node.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.01 * step, 0.0, 0.0])
            if process.IsOutputStep():
                process.PrintOutput()
                printed += 1
        return printed

    def _build_process(self, tmp_dir, extra_parameters=""):
        settings = KM.Parameters("""{
            "python_module" : "pyvista_animation_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "PyVistaAnimationOutputProcess",
            "Parameters"    : {
                "model_part_name"    : "Main",
                "output_path"        : "PLACEHOLDER",
                "file_name"          : "anim.png",
                "variable_name"      : "TEMPERATURE",
                "output_control_type": "step",
                "output_interval"    : 1
                """ + extra_parameters + """
            }
        }""")
        settings["Parameters"]["output_path"].SetString(tmp_dir)
        return pv_animation_process.Factory(settings, self.model)

    def test_Factory_returns_output_process(self):
        self._make_model_part()
        with tempfile.TemporaryDirectory() as tmp_dir:
            process = self._build_process(tmp_dir)
            self.assertIsInstance(process, KM.OutputProcess)
            self.assertEqual(process.Check(), 0)

    def test_writes_a_frame_per_output_step(self):
        mp = self._make_model_part()
        with tempfile.TemporaryDirectory() as tmp_dir:
            process = self._build_process(tmp_dir)
            process.ExecuteInitialize()
            printed = self._run_steps(process, mp, 4)
            process.ExecuteFinalize()

            self.assertEqual(printed, 4)
            self.assertEqual(len(glob.glob(os.path.join(tmp_dir, "anim_*.png"))), 4)

    def test_output_interval_is_respected(self):
        """output_interval=2 over 8 steps yields 4 frames."""
        mp = self._make_model_part()
        with tempfile.TemporaryDirectory() as tmp_dir:
            settings = KM.Parameters("""{
                "python_module" : "pyvista_animation_output_process",
                "kratos_module" : "KratosMultiphysics",
                "Parameters"    : {
                    "model_part_name"    : "Main",
                    "output_path"        : "PLACEHOLDER",
                    "file_name"          : "anim.png",
                    "variable_name"      : "TEMPERATURE",
                    "output_control_type": "step",
                    "output_interval"    : 2
                }
            }""")
            settings["Parameters"]["output_path"].SetString(tmp_dir)
            process = pv_animation_process.Factory(settings, self.model)

            process.ExecuteInitialize()
            printed = self._run_steps(process, mp, 8)
            process.ExecuteFinalize()

            self.assertEqual(printed, 4)
            self.assertEqual(len(glob.glob(os.path.join(tmp_dir, "anim_*.png"))), 4)

    def test_warp_and_clim_parameters(self):
        """Warping and an explicit color range come through from the JSON."""
        mp = self._make_model_part()
        with tempfile.TemporaryDirectory() as tmp_dir:
            process = self._build_process(
                tmp_dir,
                extra_parameters=""",
                "warp_by_vector" : "DISPLACEMENT",
                "warp_factor"    : 5.0,
                "clim"           : [0.0, 40.0],
                "component"      : 0,
                "cmap"           : "viridis\""""
            )
            process.ExecuteInitialize()
            self._run_steps(process, mp, 3)
            process.ExecuteFinalize()
            self.assertEqual(len(glob.glob(os.path.join(tmp_dir, "anim_*.png"))), 3)

    def test_default_parameters_are_complete(self):
        """Every documented key is present in the defaults."""
        self._make_model_part()
        with tempfile.TemporaryDirectory() as tmp_dir:
            process = self._build_process(tmp_dir)
            defaults = process.GetDefaultParameters()
            for key in ["model_part_name", "output_control_type", "output_interval",
                        "output_path", "file_name", "fps", "variable_name", "cmap",
                        "clim", "warp_by_vector", "time_annotation"]:
                self.assertTrue(defaults.Has(key), msg=f"missing default '{key}'")


if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
