import os

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils

dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
if dependencies_are_available:
    import KratosMultiphysics.StructuralMechanicsApplication


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMaterialsInputWithAccessors(KratosUnittest.TestCase):
    def _prepare_test(self, input_file = "material_with_accessor.json"):
        # Define a Model
        self.current_model = KratosMultiphysics.Model()

        self.model_part = self.current_model.CreateModelPart("Main")

        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/cube_few_elements")) #reusing the file that is already in the directory
        self.model_part_io.ReadModelPart(self.model_part)

        self.test_settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                    "materials_filename": "material_with_accessor.json"
            }
        }
        """)
        
        # Assign the real path
        self.test_settings["Parameters"]["materials_filename"].SetString(GetFilePath("auxiliar_files_for_python_unittest/materials_files/" + input_file))

    @KratosUnittest.skipUnless(dependencies_are_available,"StructuralMechanicsApplication not available")
    def test_input_cpp(self):
        self._prepare_test()

        KratosMultiphysics.ReadMaterialsUtility(self.test_settings, self.current_model)

        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, node.Y**2)
            node.SetSolutionStepValue(KratosMultiphysics.VISCOSITY, node.X + node.Y + node.Z)

        prop = self.model_part.GetProperties()[1]
        N = KratosMultiphysics.Vector([0.25, 0.25, 0.25, 0.25])
        for elem in self.model_part.Elements:
            geom = elem.GetGeometry()
            reference_temperature = 0.0
            reference_viscosity = 0.0
            for i in range(geom.PointsNumber()):
                reference_temperature += geom[i].Y**2 * N[i]
                reference_viscosity += (geom[i].X + geom[i].Y + geom[i].Z )* N[i]

            #check database accessor
            visc = prop.GetValue(KratosMultiphysics.VISCOSITY, geom, N, self.model_part.ProcessInfo)
            self.assertAlmostEqual(reference_viscosity, visc)

            #check table accessor
            Eref = prop.GetTable(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.YOUNG_MODULUS).GetValue(reference_temperature)
            E = prop.GetValue(KratosMultiphysics.YOUNG_MODULUS, geom, N, self.model_part.ProcessInfo)
            self.assertAlmostEqual(reference_viscosity, visc)



if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()