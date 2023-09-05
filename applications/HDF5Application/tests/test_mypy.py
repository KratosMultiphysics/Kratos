import subprocess
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTypeHinting(KratosUnittest.TestCase):
    def test_TypeHinting(self):
        install_path = Kratos.KratosPaths.kratos_install_path
        self.assertEqual(subprocess.run(["mypy", f"{install_path}/KratosMultiphysics/HDF5Application", "--check-untyped-defs"]).returncode, 0)

if __name__ == "__main__":
    KratosUnittest.main()
