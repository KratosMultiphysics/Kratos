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
        # We need a clean Kratos Model and ModelPart for each test
        self.model = KM.Model()

    def test_convert_2d_mesh(self):
        mp = self.model.CreateModelPart("Main2D")
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KM.PRESSURE)

        # Create nodes
        n1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        n3 = mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        n4 = mp.CreateNewNode(4, 1.0, 1.0, 0.0)

        # Assign solution step data values
        for node in mp.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 0, 10.0 * node.Id)
            node.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.1 * node.Id, 0.2 * node.Id, 0.0])

        # Assign initial vs current coordinates to test deformed/undeformed
        for node in mp.Nodes:
            node.X = node.X0 + 0.5
            node.Y = node.Y0 + 0.5

        # Create elements
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewElement("Element2D3N", 2, [2, 4, 3], prop)

        # 1. Test undeformed configuration
        grid = pv_utils.ModelPartToPyVista(
            mp,
            useDeformedConfiguration=False,
            nodalVariables=[KM.PRESSURE, KM.DISPLACEMENT]
        )

        self.assertIsInstance(grid, pv.UnstructuredGrid)
        self.assertEqual(grid.n_points, 4)
        self.assertEqual(grid.n_cells, 2)

        # Verify initial coordinates
        expected_points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0]
        ])
        self.assertTrue(np.allclose(grid.points, expected_points))

        # Verify point data (nodal variables)
        self.assertIn("PRESSURE", grid.point_data)
        self.assertIn("DISPLACEMENT", grid.point_data)

        expected_pressure = np.array([10.0, 20.0, 30.0, 40.0])
        expected_displacement = np.array([
            [0.1, 0.2, 0.0],
            [0.2, 0.4, 0.0],
            [0.3, 0.6, 0.0],
            [0.4, 0.8, 0.0]
        ])
        self.assertTrue(np.allclose(grid.point_data["PRESSURE"], expected_pressure))
        self.assertTrue(np.allclose(grid.point_data["DISPLACEMENT"], expected_displacement))

        # 2. Test deformed configuration
        grid_deformed = pv_utils.ModelPartToPyVista(
            mp,
            useDeformedConfiguration=True,
            nodalVariables=[KM.PRESSURE]
        )
        expected_points_deformed = np.array([
            [0.5, 0.5, 0.0],
            [1.5, 0.5, 0.0],
            [0.5, 1.5, 0.0],
            [1.5, 1.5, 0.0]
        ])
        self.assertTrue(np.allclose(grid_deformed.points, expected_points_deformed))

    def test_convert_3d_mesh(self):
        mp = self.model.CreateModelPart("Main3D")
        mp.AddNodalSolutionStepVariable(KM.TEMPERATURE)

        # Create nodes
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        mp.CreateNewNode(4, 0.0, 0.0, 1.0)

        # Create element
        prop = mp.GetProperties()[0]
        elem = mp.CreateNewElement("Element3D4N", 1, [1, 2, 3, 4], prop)

        # Assign non-historical elemental variable
        elem.SetValue(KM.DENSITY, 7800.0)

        grid = pv_utils.ModelPartToPyVista(
            mp,
            elementVariables=[KM.DENSITY]
        )

        self.assertIsInstance(grid, pv.UnstructuredGrid)
        self.assertEqual(grid.n_points, 4)
        self.assertEqual(grid.n_cells, 1)

        # Verify density in cell_data
        self.assertIn("DENSITY", grid.cell_data)
        self.assertAlmostEqual(grid.cell_data["DENSITY"][0], 7800.0)

    def test_multiblock_conversion(self):
        mp = self.model.CreateModelPart("MultiBlockPart")

        # Nodes
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)

        # Elements & Conditions
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewCondition("LineCondition2D2N", 1, [1, 2], prop)

        # Convert both
        blocks = pv_utils.ModelPartToPyVista(
            mp,
            exportElements=True,
            exportConditions=True
        )

        self.assertIsInstance(blocks, pv.MultiBlock)
        self.assertIn("elements", blocks.keys())
        self.assertIn("conditions", blocks.keys())

        elem_grid = blocks["elements"]
        cond_grid = blocks["conditions"]

        self.assertEqual(elem_grid.n_cells, 1)
        self.assertEqual(cond_grid.n_cells, 1)

        # Element cell type VTK_TRIANGLE (5)
        self.assertEqual(elem_grid.celltypes[0], 5)
        # Condition cell type VTK_LINE (3)
        self.assertEqual(cond_grid.celltypes[0], 3)

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

    def test_ComputeWarpedMesh(self):
        mp = self.model.CreateModelPart("WarpPart")
        mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        n1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        n3 = mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        prop = mp.GetProperties()[0]
        mp.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)

        n1.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.1, 0.2, 0.0])
        n2.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.1, 0.2, 0.0])
        n3.SetSolutionStepValue(KM.DISPLACEMENT, 0, [0.1, 0.2, 0.0])

        warped_grid = pv_utils.ComputeWarpedMesh(mp, KM.DISPLACEMENT, factor=2.0)
        self.assertIsInstance(warped_grid, pv.UnstructuredGrid)
        self.assertEqual(warped_grid.n_points, 3)
        expected_warped_points = np.array([
            [0.2, 0.4, 0.0],
            [1.2, 0.4, 0.0],
            [0.2, 1.4, 0.0]
        ])
        self.assertTrue(np.allclose(warped_grid.points, expected_warped_points))

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
        self.assertTrue(len(slices) > 0)

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
            mp,
            velocityVariable=KM.VELOCITY,
            sourceCenter=[0.1, 0.5, 0.5],
            sourceRadius=0.1,
            nPoints=5
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

        # Test above
        threshed_above = pv_utils.CreateThresholdedMesh(mp, KM.TEMPERATURE, thresholdValue=4.5, thresholdType="above")
        self.assertIsInstance(threshed_above, pv.UnstructuredGrid)
        # Test below
        threshed_below = pv_utils.CreateThresholdedMesh(mp, KM.TEMPERATURE, thresholdValue=4.5, thresholdType="below")
        self.assertIsInstance(threshed_below, pv.UnstructuredGrid)
        # Test between
        threshed_between = pv_utils.CreateThresholdedMesh(mp, KM.TEMPERATURE, thresholdValue=[2.5, 6.5], thresholdType="between")
        self.assertIsInstance(threshed_between, pv.UnstructuredGrid)

    def test_CreateVectorGlyphs(self):
        mp = self.model.CreateModelPart("GlyphPart")
        mp.AddNodalSolutionStepVariable(KM.VELOCITY)
        for i, coords in enumerate([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]
        ], 1):
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

        clipped = pv_utils.CreateClippedMesh(mp, normal=[1.0, 0.0, 0.0], origin=[0.5, 0.5, 0.5])
        self.assertIsInstance(clipped, pv.UnstructuredGrid)

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
