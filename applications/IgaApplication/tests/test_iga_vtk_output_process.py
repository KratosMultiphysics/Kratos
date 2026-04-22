import os
import h5py
import numpy as np
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication
from KratosMultiphysics import KratosUnittest

# Returns path to files stored in the dedicated test data folder
def GetFilePath(file_name):
    return os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "iga_vtk_output_process_test",
        file_name
    )

class TestIgaVTKHDFOutputProcess(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Base folder of this test
        cls.work_folder = os.path.dirname(os.path.realpath(__file__))

        # Reference file (checked against) and output file (generated during test)
        cls.ref_file = os.path.join(
            cls.work_folder,
            "iga_vtk_output_process_test",
            "reference_output.vtkhdf"
        )
        cls.out_file = os.path.join(cls.work_folder, "test_output.vtkhdf")

        # Create model and model part
        cls.model = KM.Model()
        cls.model_part = cls.model.CreateModelPart("ModelPart")

        # Minimal buffer is enough since we only access current step values
        cls.model_part.SetBufferSize(1)

        # Register variables to be written by the output process
        cls.model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        # Read Brep geometry from JSON
        input_file = GetFilePath("2_patch_geometry.cad.json")
        input_file = os.path.abspath(input_file)
        KM.CadJsonInput(input_file).ReadModelPart(cls.model_part)

        # Use first geometry (single-surface test case)
        cls.brep_id = next(iter(cls.model_part.Geometries)).Id

        # Assign a deterministic displacement field (varies per node)
        # This avoids constant fields and makes the test more sensitive
        i = 0
        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT, [1.0*i, 2.0*i, 3.0*i])
            i += 1

    @classmethod
    def tearDownClass(cls):
        # Clean up generated output file
        if os.path.exists(cls.out_file):
            os.remove(cls.out_file)

    def test_output_against_reference(self):
        # Ensure reference exists before comparing
        self.assertTrue(os.path.exists(self.ref_file))

        # Initialize time
        self.model_part.ProcessInfo[KM.TIME] = 0.0

        # Settings for the VTKHDF output process
        settings = KM.Parameters(f"""
        {{
            "model_part_name" : "ModelPart",
            "output_file_name" : "{self.out_file}",
            "brep_surface_ids" : [{self.brep_id}],
            "nodal_solution_step_data_variables" : ["DISPLACEMENT"],
            "output_refinement" : [3,3],
            "output_control_type" : "step",
            "output_frequency" : 1
        }}
        """)

        # Create process
        from KratosMultiphysics.IgaApplication.iga_vtk_output_process import IgaVTKOutputProcess
        process = IgaVTKOutputProcess(self.model, settings)

        # Standard Kratos process lifecycle
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        # First output (time = 0.0)
        if process.IsOutputStep():
            process.PrintOutput()

        # Advance time and trigger second output
        self.model_part.ProcessInfo[KM.TIME] = 1.0
        process.ExecuteFinalizeSolutionStep()

        if process.IsOutputStep():
            process.PrintOutput()

        # Compare generated file with reference
        with h5py.File(self.ref_file, "r") as f_ref, \
             h5py.File(self.out_file, "r") as f_out:

            self._compare_groups(f_ref, f_out)

    def _compare_groups(self, g_ref, g_out):
        # Compare group structure
        self.assertEqual(set(g_ref.keys()), set(g_out.keys()))

        # Compare attributes (handle numpy arrays explicitly)
        for key in g_ref.attrs.keys():
            self.assertIn(key, g_out.attrs)

            val_ref = g_ref.attrs[key]
            val_out = g_out.attrs[key]

            if isinstance(val_ref, np.ndarray):
                np.testing.assert_array_equal(val_ref, val_out)
            else:
                self.assertEqual(val_ref, val_out)

        # Recursively compare datasets and subgroups
        for key in g_ref.keys():

            obj_ref = g_ref[key]
            obj_out = g_out[key]

            if isinstance(obj_ref, h5py.Dataset):

                data_ref = obj_ref[()]
                data_out = obj_out[()]

                # Shape must match exactly
                self.assertEqual(data_ref.shape, data_out.shape)

                # Floating-point comparison with tolerance
                if np.issubdtype(data_ref.dtype, np.floating):
                    np.testing.assert_allclose(data_ref, data_out, rtol=1e-9, atol=1e-12)
                else:
                    np.testing.assert_array_equal(data_ref, data_out)

            elif isinstance(obj_ref, h5py.Group):
                self._compare_groups(obj_ref, obj_out)

if __name__ == "__main__":
    KratosUnittest.main()