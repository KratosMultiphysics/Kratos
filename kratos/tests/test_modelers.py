import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.modeler_factory import KratosModelerFactory


def run_modelers(current_model, modelers_list):
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(
        current_model,
        modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()


class TestModelers(KratosUnittest.TestCase):
    def test_mdpa_io_modeler(self):

        modelers_list = KM.Parameters(
        """ [{
            "modeler_name": "MdpaIoModeler",
            "Parameters": {
                "echo_level"                                 : 0,
                "input_file_name"                            : "auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read",
                "model_part_name"                            : "model_part",
                "skip_timer"                                 : true,
                "ignore_variables_not_in_solution_step_data" : false
            }
        }] """)

        current_model = KM.Model()

        model_part = current_model.CreateModelPart("model_part")
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)

        run_modelers(current_model, modelers_list)

        self.assertEqual(model_part.NumberOfSubModelParts(), 2)
        self.assertEqual(model_part.NumberOfTables(), 1)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)

        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KM.DISPLACEMENT_Y), 0.000974)

        # The other features of the model part import are tested in test_model_part_io.py

if __name__ == '__main__':
    KratosUnittest.main()
