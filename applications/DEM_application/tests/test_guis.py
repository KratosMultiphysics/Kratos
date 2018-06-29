import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestGUIs(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_GUIs_1(self):
        import filecmp
        f1 = "../custom_problemtype/DEMpack/G-DEMPack/kratos.gid/python/KratosDEM.py"
        f2 = "../python_scripts/KratosDEM.py"
        self.assertTrue(filecmp.cmp(f1, f2))

    def tearDown(self):
        pass


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()