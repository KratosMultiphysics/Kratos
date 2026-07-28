# Import Kratos Multiphysics and the unittest framework
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import NumPy and standard modules
import numpy as np
import os
import tempfile

# Attempt to import PyVista and the utilities module, set a flag if missing
try:
    import pyvista as pv
    import KratosMultiphysics.pyvista_utilities as pv_utils
    missing_pyvista = False
except ImportError:
    missing_pyvista = True


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


if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
