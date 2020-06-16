from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.IgaApplication as IgaApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestModelers(KratosUnittest.TestCase):
    def _run_modelers(current_model, modelers_list):
        list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

        for modeler in list_of_modelers:
            modeler.SetupGeometryModel()

        for modeler in list_of_modelers:
            modeler.PrepareGeometryModel()

        for modeler in list_of_modelers:
            modeler.SetupModelPart()


    def _strong_support_surface_test(self, current_model):

        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "CadIoModeler",
            "Parameters": {
                "echo_level": 4,
                "cad_model_part_name": "IgaModelPart",
                "physics_file_name": "modeler_tests/surface_geometry.cad.json"
            } }, {
            "modeler_name": "IgaModeler",
            "Parameters": {
                "echo_level":  4,
                "cad_model_part_name": "IgaModelPart",
                "analysis_model_part_name": "IgaModelPart",
                "physics_file_name": "modeler_tests/strong_support_physics.iga.json"
            } }] """ )

        _run_modelers(current_model, modelers_list)

        support_1_model_part = current_model.GetModelPart("IgaModelPart.Support_1")
        support_2_model_part = current_model.GetModelPart("IgaModelPart.Support_2")

        self.assertEqual(support_1_model_part.NumberOfNodes(), 2)
        self.assertEqual(support_1_model_part.Nodes()[0], 6)
        self.assertEqual(support_1_model_part.Nodes()[0], 12)

        self.assertEqual(support_2_model_part.NumberOfNodes(), 6)

    def test_StrongSupportSurface(self):
        current_model = KratosMultiphysics.Model()
        self._strong_support_surface_test(current_model)

if __name__ == '__main__':
    KratosUnittest.main()
