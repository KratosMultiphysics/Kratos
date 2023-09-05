import subprocess
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestTypeHinting(KratosUnittest.TestCase):
    def test_TypeHinting(self):
        try:
            import mypy
        except:
            self.skipTest("Requires h5py.")

        install_path = Kratos.KratosPaths.kratos_install_path
        self.assertEqual(subprocess.run(["mypy", f"{install_path}/KratosMultiphysics/HDF5Application"]).returncode, 0)

if __name__ == "__main__":
    KratosUnittest.main()
