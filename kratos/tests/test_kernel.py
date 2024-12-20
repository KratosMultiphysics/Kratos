import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

class TestKernel(KratosUnittest.TestCase):
    def setUp(self):
      pass

    def testIsLibraryAvailable(self):
      # Check if the libraries are available
      self.assertFalse(KM.KratosGlobals.Kernel.IsLibraryAvailable("pikachu_lib"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("boost"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("amgcl"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("benchmark"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("clipper"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("concurrentqueue"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("ghc"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("gidpost"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("intrusive_ptr"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("json"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("pybind11"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("span"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("tinyexpr"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("vexcl"))
      self.assertTrue(KM.KratosGlobals.Kernel.IsLibraryAvailable("zlib"))
      self.assertNotEqual(KM.KratosGlobals.Kernel.IsLibraryAvailable("triangle"), KM.KratosGlobals.Kernel.IsLibraryAvailable("delaunator-cpp"))

if __name__ == "__main__":
    KratosUnittest.main()
