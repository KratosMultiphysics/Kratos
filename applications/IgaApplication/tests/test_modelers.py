import KratosMultiphysics
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

def run_modelers(current_model, modelers_list):
    from KratosMultiphysics.modeler_factory import KratosModelerFactory
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()

class TestModelers(KratosUnittest.TestCase):
    def _strong_support_surface_test(self, current_model):

        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "CadIoModeler",
            "Parameters": {
                "echo_level": 0,
                "cad_model_part_name": "IgaModelPart",
                "geometry_file_name": "modeler_tests/surface_geometry.cad.json"
            } }, {
            "modeler_name": "IgaModeler",
            "Parameters": {
                "echo_level":  0,
                "cad_model_part_name": "IgaModelPart",
                "analysis_model_part_name": "IgaModelPart",
                "physics_file_name": "modeler_tests/strong_support_physics.iga.json"
            } }] """ )

        run_modelers(current_model, modelers_list)

        support_1_model_part = current_model.GetModelPart("IgaModelPart.Support_1")
        support_1_Variation_model_part = current_model.GetModelPart("IgaModelPart.Support_1_Variation")
        support_2_model_part = current_model.GetModelPart("IgaModelPart.Support_2")
        support_2_Variation_model_part = current_model.GetModelPart("IgaModelPart.Support_2_Variation")

        # Check if all needed node are within the model parts
        self.assertEqual(support_1_model_part.NumberOfNodes(), 2)
        self.assertEqual(support_1_model_part.GetNodes()[6].Id, 6)
        self.assertEqual(support_1_model_part.GetNodes()[12].Id, 12)

        self.assertEqual(support_1_Variation_model_part.NumberOfNodes(), 2)
        self.assertEqual(support_1_Variation_model_part.GetNodes()[5].Id, 5)
        self.assertEqual(support_1_Variation_model_part.GetNodes()[11].Id, 11)

        self.assertEqual(support_2_model_part.NumberOfNodes(), 6)
        self.assertEqual(support_2_model_part.GetNodes()[7].Id, 7)
        self.assertEqual(support_2_model_part.GetNodes()[8].Id, 8)
        self.assertEqual(support_2_model_part.GetNodes()[9].Id, 9)
        self.assertEqual(support_2_model_part.GetNodes()[10].Id, 10)
        self.assertEqual(support_2_model_part.GetNodes()[11].Id, 11)
        self.assertEqual(support_2_model_part.GetNodes()[12].Id, 12)

        self.assertEqual(support_2_Variation_model_part.NumberOfNodes(), 6)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[1].Id, 1)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[2].Id, 2)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[3].Id, 3)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[4].Id, 4)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[5].Id, 5)
        self.assertEqual(support_2_Variation_model_part.GetNodes()[6].Id, 6)

    def test_StrongSupportSurface(self):
        current_model = KratosMultiphysics.Model()
        self._strong_support_surface_test(current_model)

if __name__ == '__main__':
    KratosUnittest.main()
