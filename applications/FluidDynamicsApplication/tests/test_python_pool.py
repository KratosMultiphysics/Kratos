import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.FluidDynamicsApplication
import numpy as np
import KratosMultiphysics.FluidDynamicsApplication.python_pool

try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

class TestBufferPoolNumPy(UnitTest.TestCase):

    def setUp(self):
        self.pool = KratosMultiphysics.FluidDynamicsApplication.python_pool.BufferPool(
            sizes=[100, 200],
            xp=np,
            dtype=np.float64
        )

    def test_Get_basic(self):
        buf = self.pool.Get(0, (10, 10))
        self.assertEqual(buf.shape, (10, 10))
        self.assertEqual(buf.dtype, np.float64)
        self.assertTrue(self.pool.in_use[0])

    def test_Get_reshape_is_view(self):
        buf = self.pool.Get(0, (20, 5))
        self.assertIs(buf.base, self.pool.buffers[0])
        self.pool.Release(0)

    def test_Get_exceeds_size(self):
        with self.assertRaises(ValueError):
            self.pool.Get(0, (50, 50))  # requires 2500 > 100

    def test_Get_smaller_than_required(self):
        # Request shape that fits but check that the returned buffer
        # is exactly the required size, not smaller.
        buf = self.pool.Get(0, (5, 5))  # requires 25
        self.assertEqual(buf.size, 25)
        self.assertEqual(buf.shape, (5, 5))
        self.pool.Release(0)

    def test_double_Get_raises(self):
        self.pool.Get(0, (5, 5))
        with self.assertRaises(RuntimeError):
            self.pool.Get(0, (2, 2))

    def test_Release_basic(self):
        self.pool.Get(1, (20, 10))
        self.pool.Release(1)
        self.assertFalse(self.pool.in_use[1])

    def test_Release_not_in_use(self):
        with self.assertRaises(RuntimeError):
            self.pool.Release(0)

    def test_multiple_buffers_independent(self):
        a = self.pool.Get(0, (10, 10))
        b = self.pool.Get(1, (20, 10))
        self.assertTrue(self.pool.in_use[0])
        self.assertTrue(self.pool.in_use[1])
        self.pool.Release(0)
        self.pool.Release(1)

    def test_shape_from_numpy_array(self):
        shape = np.array([10, 10])
        buf = self.pool.Get(0, shape)
        self.assertEqual(buf.shape, (10, 10))
        self.pool.Release(0)

    def test_release_by_array_numpy(self):
        """Test releasing a buffer by passing the array itself (auto-detect index)."""
        buf = self.pool.Get(0, (10, 10))
        self.assertTrue(self.pool.in_use[0])
        # Release by passing the array instead of index
        self.pool.Release(buf)
        self.assertFalse(self.pool.in_use[0])

    def test_release_by_array_different_buffers(self):
        """Test releasing different buffers by passing their arrays."""
        buf0 = self.pool.Get(0, (5, 5))
        buf1 = self.pool.Get(1, (10, 5))
        self.assertTrue(self.pool.in_use[0])
        self.assertTrue(self.pool.in_use[1])
        # Release buffer 1 by passing its array
        self.pool.Release(buf1)
        self.assertFalse(self.pool.in_use[1])
        self.assertTrue(self.pool.in_use[0])
        # Release buffer 0 by passing its array
        self.pool.Release(buf0)
        self.assertFalse(self.pool.in_use[0])

    def test_release_by_array_raises_if_not_in_use(self):
        """Test that releasing an array that is not in use raises RuntimeError."""
        buf = self.pool.Get(0, (5, 5))
        self.pool.Release(buf)  # First release succeeds
        with self.assertRaises(RuntimeError):
            self.pool.Release(buf)  # Second release should raise


# --- CuPy Test Suite (Kratos style, auto-skip) ---

@UnitTest.skipUnless(HAS_CUPY, "CuPy not available")
class TestBufferPoolCuPy(UnitTest.TestCase):

    def setUp(self):
        self.pool = KratosMultiphysics.FluidDynamicsApplication.python_pool.BufferPool(
            sizes=[500, 1000],
            xp=cp,
            dtype=cp.float32
        )

    def test_Get_basic(self):
        buf = self.pool.Get(0, (10, 10))
        self.assertEqual(buf.shape, (10, 10))
        self.assertEqual(buf.dtype, cp.float32)
        self.assertTrue(self.pool.in_use[0])
        self.pool.Release(0)

    def test_Get_exceeds_size(self):
        with self.assertRaises(ValueError):
            self.pool.Get(0, (30, 30))  # requires 900 > 500

    def test_shape_from_cupy_array(self):
        shape = cp.array([10, 10])
        buf = self.pool.Get(0, shape)
        self.assertEqual(buf.shape, (10, 10))
        self.pool.Release(0)

    def test_Get_smaller_than_required(self):
        buf = self.pool.Get(1, (5, 5))
        self.assertEqual(buf.size, 25)
        self.assertEqual(buf.shape, (5, 5))
        self.pool.Release(1)

    def test_double_Get_raises(self):
        self.pool.Get(0, (5, 5))
        with self.assertRaises(RuntimeError):
            self.pool.Get(0, (2, 2))
        self.pool.Release(0)

    def test_release_by_array_cupy(self):
        """Test releasing a buffer by passing the CuPy array itself (auto-detect index)."""
        buf = self.pool.Get(0, (10, 10))
        self.assertTrue(self.pool.in_use[0])
        # Release by passing the array instead of index
        self.pool.Release(buf)
        self.assertFalse(self.pool.in_use[0])

    def test_release_by_array_different_buffers_cupy(self):
        """Test releasing different CuPy buffers by passing their arrays."""
        buf0 = self.pool.Get(0, (5, 5))
        buf1 = self.pool.Get(1, (10, 5))
        self.assertTrue(self.pool.in_use[0])
        self.assertTrue(self.pool.in_use[1])
        # Release buffer 1 by passing its array
        self.pool.Release(buf1)
        self.assertFalse(self.pool.in_use[1])
        self.assertTrue(self.pool.in_use[0])
        # Release buffer 0 by passing its array
        self.pool.Release(buf0)
        self.assertFalse(self.pool.in_use[0])

    def test_release_by_array_raises_if_not_in_use_cupy(self):
        """Test that releasing a CuPy array that is not in use raises RuntimeError."""
        buf = self.pool.Get(0, (5, 5))
        self.pool.Release(buf)  # First release succeeds
        with self.assertRaises(RuntimeError):
            self.pool.Release(buf)  # Second release should raise


if __name__ == "__main__":
    UnitTest.main()