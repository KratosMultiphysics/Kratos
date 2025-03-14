import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

def run_modelers(current_model, modelers_list):
    from KratosMultiphysics.modeler_factory import KratosModelerFactory
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

class TestCadTessellationModeler(KratosUnittest.TestCase):
    @KratosUnittest.skipUnless(KratosMultiphysics.KratosGlobals.Kernel.IsLibraryAvailable("triangle"), "Kratos compiled without triangle")
    def test_cad_tessellation_modeler(self):

        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "CadIoModeler",
            "Parameters": {
                "echo_level": 0,
                "cad_model_part_name": "CadModelPart",
                "geometry_file_name":
                    "auxiliar_files_for_python_unittest/cad_json_files/cylinder.cad.json"
            } }, {
            "modeler_name": "CadTessellationModeler",
            "Parameters": {
<<<<<<< HEAD
                "echo_level": 0,
                "cad_model_part_name": "CadModelPart",
                "skin_model_part_name": "SkinModelPart",
                "absolute_chordal_error"   : 1e-1,
=======
                "echo_level"                    : 1,
                "cad_model_part_name"           : "CadModelPart",
                "skin_model_part_name"          : "SkinModelPart",
                "absolute_chordal_error"        : 1e-1,
>>>>>>> master
                "absolute_triangulation_error"  : 1e-1,
                "initial_triangle_area"         : 5,
                "max_triangulation_iteration"   : 1
            } }] """)

        current_model = KratosMultiphysics.Model()
        run_modelers(current_model, modelers_list)

        skin_model_part = current_model.GetModelPart("SkinModelPart")

        # Check if all needed nodes are within the SkinModelPart
        self.assertEqual(skin_model_part.NumberOfNodes(), 141)
        self.assertEqual(skin_model_part.NumberOfElements(), 47)

        self.assertEqual(skin_model_part.GetNodes()[32].Id, 32)
        self.assertEqual(skin_model_part.GetElements()[11].Id, 11)


        self.assertAlmostEqual(skin_model_part.GetNodes()[0].X, -3.53553390593274)
        self.assertAlmostEqual(skin_model_part.GetNodes()[1].X, -3.8272589577240357)
        self.assertAlmostEqual(skin_model_part.GetNodes()[2].X, -1.6682976074902058)
        self.assertAlmostEqual(skin_model_part.GetNodes()[3].X, 3.589564999704995)

        self.assertAlmostEqual(skin_model_part.GetNodes()[0].Y, 3.535533905932736)
        self.assertAlmostEqual(skin_model_part.GetNodes()[1].Y, 3.217466219017866)
        self.assertAlmostEqual(skin_model_part.GetNodes()[2].Y, 4.713468265814724)
        self.assertAlmostEqual(skin_model_part.GetNodes()[3].Y, 3.480664176977272)

        self.assertAlmostEqual(skin_model_part.GetNodes()[0].Z, 0.0)
        self.assertAlmostEqual(skin_model_part.GetNodes()[1].Z, 2.5227380849575387)
        self.assertAlmostEqual(skin_model_part.GetNodes()[2].Z, 2.554143073766855)
        self.assertAlmostEqual(skin_model_part.GetNodes()[3].Z, 5.0)

        self.assertAlmostEqual(skin_model_part.GetNodes()[50].X, -4.929871219389899)
        self.assertAlmostEqual(skin_model_part.GetNodes()[51].X, -4.662394769221466)
        self.assertAlmostEqual(skin_model_part.GetNodes()[52].X, -4.929871219389899)
        self.assertAlmostEqual(skin_model_part.GetNodes()[53].X, -5.0)

        self.assertAlmostEqual(skin_model_part.GetNodes()[50].Y, 0.8344877232357282)
        self.assertAlmostEqual(skin_model_part.GetNodes()[51].Y, -1.8061215396357764)
        self.assertAlmostEqual(skin_model_part.GetNodes()[52].Y, 0.8344877232357282)
        self.assertAlmostEqual(skin_model_part.GetNodes()[53].Y, -6.123233995736766e-16)

        self.assertAlmostEqual(skin_model_part.GetNodes()[50].Z, 2.455843399391155)
        self.assertAlmostEqual(skin_model_part.GetNodes()[51].Z, 1.0104686519865287)
        self.assertAlmostEqual(skin_model_part.GetNodes()[52].Z, 2.455843399391155)
        self.assertAlmostEqual(skin_model_part.GetNodes()[53].Z, 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
