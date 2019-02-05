import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestGUIs(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_GUIs_1(self):
        import filecmp
        f1 = "../custom_problemtype/DEMpack/G-DEMPack/kratos.gid/python/KratosDEM.py"
        f2 = "../python_scripts/KratosDEM.py"
        f1_absolute_path = GetFilePath(f1)
        f2_absolute_path = GetFilePath(f2)
        self.assertTrue(filecmp.cmp(f1_absolute_path, f2_absolute_path))

    def tearDown(self):
        pass


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()