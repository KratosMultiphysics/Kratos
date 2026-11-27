import os
import numpy as np
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

try:
    import h5py
    import KratosMultiphysics.gid_hdf5_utilities as gid_hdf5
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False
    print("Skipping tests that require h5py.")


@KratosUnittest.skipUnless(HAS_H5PY, "h5py not available")
class TestGidHdf5Utilities(KratosUnittest.TestCase):

    def setUp(self):
        self.file_path = "gid_hdf5_utilities_test.h5"
        KratosUtilities.DeleteFileIfExisting(self.file_path)

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(self.file_path)

    def _write_fixture(self, str_attributes=True):
        """Create a fixture mimicking the GID-HDF5 output layout:
        Results/<group> with 'Step' and 'Name' attributes and numbered datasets.
        Dataset '1' holds ids, and a group with a '4' dataset is treated as vector data.
        """
        def set_attr(group, name, value):
            if str_attributes:
                group.attrs[name] = value
            else:
                # Force a bytes-typed attribute (older h5py versions read these back as bytes)
                group.attrs.create(name, data=np.bytes_(value))

        with h5py.File(self.file_path, "w") as f:
            results = f.create_group("Results")

            # Scalar PRESSURE at step 0.0
            g_scalar = results.create_group("group_1")
            set_attr(g_scalar, "Step", "0.0")
            set_attr(g_scalar, "Name", "PRESSURE")
            g_scalar.create_dataset("1", data=np.array([1, 2, 3], dtype=np.int32))
            g_scalar.create_dataset("2", data=np.array([10.0, 20.0, 30.0]))

            # Vector VELOCITY at step 0.0
            g_vec = results.create_group("group_2")
            set_attr(g_vec, "Step", "0.0")
            set_attr(g_vec, "Name", "VELOCITY")
            g_vec.create_dataset("1", data=np.array([1, 2], dtype=np.int32))
            g_vec.create_dataset("2", data=np.array([1.0, 4.0]))
            g_vec.create_dataset("3", data=np.array([2.0, 5.0]))
            g_vec.create_dataset("4", data=np.array([3.0, 6.0]))

    def test_decode_attribute_handles_bytes_and_str(self):
        self.assertEqual(gid_hdf5._decode_attribute(b"0.0"), "0.0")
        self.assertEqual(gid_hdf5._decode_attribute("0.0"), "0.0")

    def test_find_scalar_dataset(self):
        self._write_fixture(str_attributes=True)
        ids, data = gid_hdf5.find_dataset_by_step_and_variable(self.file_path, 0.0, "PRESSURE")
        np.testing.assert_array_equal(ids, np.array([1, 2, 3], dtype=np.int32))
        self.assertEqual(ids.dtype, np.int32)
        self.assertEqual(data.shape, (3, 1))
        np.testing.assert_allclose(data[:, 0], [10.0, 20.0, 30.0])

    def test_find_vector_dataset(self):
        self._write_fixture(str_attributes=True)
        ids, data = gid_hdf5.find_dataset_by_step_and_variable(self.file_path, 0.0, "VELOCITY")
        np.testing.assert_array_equal(ids, np.array([1, 2], dtype=np.int32))
        self.assertEqual(data.shape, (2, 3))
        np.testing.assert_allclose(data, np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))

    def test_find_with_bytes_attributes(self):
        # Regression: the previous implementation called .decode() unconditionally, which
        # raises AttributeError on str-typed attributes (h5py >= 3) and silently returned None
        # through the broad exception handler. Both attribute storage types must work.
        self._write_fixture(str_attributes=False)
        ids, data = gid_hdf5.find_dataset_by_step_and_variable(self.file_path, 0.0, "PRESSURE")
        np.testing.assert_array_equal(ids, np.array([1, 2, 3], dtype=np.int32))
        np.testing.assert_allclose(data[:, 0], [10.0, 20.0, 30.0])

    def test_dtype_argument_is_applied(self):
        self._write_fixture(str_attributes=True)
        _, data = gid_hdf5.find_dataset_by_step_and_variable(
            self.file_path, 0.0, "VELOCITY", dtype=np.float64)
        np.testing.assert_allclose(data, np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))

    def test_no_match_raises_eof_error(self):
        # Regression: EOFError used to be swallowed by the broad 'except Exception' and the
        # function returned None instead of reporting the well-defined 'not found' condition.
        self._write_fixture(str_attributes=True)
        with self.assertRaises(EOFError):
            gid_hdf5.find_dataset_by_step_and_variable(self.file_path, 99.0, "PRESSURE")
        with self.assertRaises(EOFError):
            gid_hdf5.find_dataset_by_step_and_variable(self.file_path, 0.0, "NOPE")

    def test_missing_file_returns_none(self):
        self.assertIsNone(gid_hdf5.find_dataset_by_step_and_variable(
            "this_file_does_not_exist.h5", 0.0, "PRESSURE"))

    def test_missing_results_group_returns_none(self):
        with h5py.File(self.file_path, "w") as f:
            f.create_group("Meshes")
        self.assertIsNone(gid_hdf5.find_dataset_by_step_and_variable(
            self.file_path, 0.0, "PRESSURE"))

    def test_list_by_variable(self):
        self._write_fixture(str_attributes=True)
        listing = gid_hdf5.list_by_variable(self.file_path)
        self.assertEqual(sorted(listing.keys()), ["PRESSURE", "VELOCITY"])
        self.assertEqual(listing["PRESSURE"][0], ["0.0"])
        self.assertEqual(listing["PRESSURE"][1], ["group_1"])
        self.assertEqual(listing["VELOCITY"][0], ["0.0"])
        self.assertEqual(listing["VELOCITY"][1], ["group_2"])

    def _write_two_step_fixture(self):
        """Fixture with two PRESSURE steps and one VELOCITY step, plus a 'Meshes' group
        (required by ExtractFileFromBack), to check the from-back selection."""
        with h5py.File(self.file_path, "w") as f:
            f.create_group("Meshes").attrs["info"] = "dummy"
            results = f.create_group("Results")

            for group_name, step, values in [("group_1", "0.0", [10.0, 20.0]),
                                             ("group_3", "1.0", [30.0, 40.0])]:
                g = results.create_group(group_name)
                g.attrs["Step"] = step
                g.attrs["Name"] = "PRESSURE"
                g.create_dataset("1", data=np.array([1, 2], dtype=np.int32))
                g.create_dataset("2", data=np.array(values))

            g_vec = results.create_group("group_2")
            g_vec.attrs["Step"] = "0.0"
            g_vec.attrs["Name"] = "VELOCITY"
            g_vec.create_dataset("1", data=np.array([1, 2], dtype=np.int32))
            g_vec.create_dataset("2", data=np.array([1.0, 4.0]))
            g_vec.create_dataset("3", data=np.array([2.0, 5.0]))
            g_vec.create_dataset("4", data=np.array([3.0, 6.0]))

    def test_list_by_variable_sorts_by_step(self):
        self._write_two_step_fixture()
        listing = gid_hdf5.list_by_variable(self.file_path)
        self.assertEqual(listing["PRESSURE"][0], ["0.0", "1.0"])
        self.assertEqual(listing["PRESSURE"][1], ["group_1", "group_3"])

    def test_extract_from_back_default_takes_most_recent(self):
        # Regression: steps_from_back=0 (default) used to index var_data[1][-0], i.e. the FIRST
        # (oldest) group instead of the most recent one.
        self._write_two_step_fixture()
        dst_path = "gid_hdf5_utilities_test_copy.h5"
        KratosUtilities.DeleteFileIfExisting(dst_path)
        gid_hdf5.ExtractFileFromBack(self.file_path, dst_path)
        try:
            ids, data = gid_hdf5.find_dataset_by_step_and_variable(dst_path, 1.0, "PRESSURE")
            np.testing.assert_allclose(data[:, 0], [30.0, 40.0])
        finally:
            KratosUtilities.DeleteFileIfExisting(dst_path)

    def test_extract_from_back_second_to_last(self):
        self._write_two_step_fixture()
        dst_path = "gid_hdf5_utilities_test_copy.h5"
        KratosUtilities.DeleteFileIfExisting(dst_path)
        gid_hdf5.ExtractFileFromBack(self.file_path, dst_path, steps_from_back=2)
        try:
            ids, data = gid_hdf5.find_dataset_by_step_and_variable(dst_path, 0.0, "PRESSURE")
            np.testing.assert_allclose(data[:, 0], [10.0, 20.0])
        finally:
            KratosUtilities.DeleteFileIfExisting(dst_path)

    def test_extract_from_back_copies_meshes(self):
        self._write_two_step_fixture()
        dst_path = "gid_hdf5_utilities_test_copy.h5"
        KratosUtilities.DeleteFileIfExisting(dst_path)
        gid_hdf5.ExtractFileFromBack(self.file_path, dst_path)
        try:
            with h5py.File(dst_path, "r") as f:
                self.assertIn("Meshes", f)
                self.assertIn("Results", f)
                # One group per variable: PRESSURE (most recent) + VELOCITY
                self.assertEqual(sorted(f["Results"].keys()), ["group_2", "group_3"])
        finally:
            KratosUtilities.DeleteFileIfExisting(dst_path)


if __name__ == '__main__':
    KratosUnittest.main()
