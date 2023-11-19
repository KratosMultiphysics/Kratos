import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities


class TestModelPartUtiliites(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()

        self.data_settings = KM.Parameters("""{
            "disp" : {
                "model_part_name" : "Structure",
                "variable_name"   : "DISPLACEMENT"
            },
            "load" : {
                "model_part_name" : "Structure.Interface",
                "variable_name"   : "FORCE"
            },
            "reaction" : {
                "model_part_name" : "Structure.GENERIC_FSI.subinterface",
                "variable_name"   : "REACTION"
            },
            "vel" : {
                "model_part_name" : "aux_model_part",
                "variable_name"   : "VELOCITY_X"
            }
        }""")

    def test_AllocateHistoricalVariablesFromCouplingDataSettings(self):
        mp_structure = self.model.CreateModelPart("Structure")
        mp_aux = self.model.CreateModelPart("aux_model_part")

        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.data_settings, self.model, "dummy_solver_name")

        # make sure the original settings are not modified
        self.assertFalse(self.data_settings["disp"].Has("location"))

        self.assertTrue(mp_structure.HasNodalSolutionStepVariable(KM.DISPLACEMENT))
        self.assertTrue(mp_structure.HasNodalSolutionStepVariable(KM.FORCE))
        self.assertTrue(mp_structure.HasNodalSolutionStepVariable(KM.REACTION))
        self.assertFalse(mp_structure.HasNodalSolutionStepVariable(KM.VELOCITY_X))

        self.assertTrue(mp_aux.HasNodalSolutionStepVariable(KM.VELOCITY_X))
        self.assertFalse(mp_aux.HasNodalSolutionStepVariable(KM.DISPLACEMENT))
        self.assertFalse(mp_aux.HasNodalSolutionStepVariable(KM.FORCE))
        self.assertFalse(mp_aux.HasNodalSolutionStepVariable(KM.REACTION))

    def test_CreateMainModelPartsFromCouplingDataSettings(self):
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.data_settings, self.model, "dummy_solver_name")

        # make sure the original settings are not modified
        self.assertFalse(self.data_settings["disp"].Has("location"))

        self.assertTrue(self.model.HasModelPart("Structure"))
        self.assertTrue(self.model.HasModelPart("aux_model_part"))

        self.assertFalse(self.model.HasModelPart("Structure.Interface"))
        self.assertFalse(self.model.HasModelPart("Structure.GENERIC_FSI"))

        self.assertFalse(self.model.HasModelPart("Interface"))
        self.assertFalse(self.model.HasModelPart("GENERIC_FSI"))
        self.assertFalse(self.model.HasModelPart("subinterface"))

    def test_RecursiveCreateModelParts(self):
        model_part = self.model.CreateModelPart("for_test")

        self.assertFalse(model_part.HasSubModelPart("sub_model_part"))
        model_part_utilities.RecursiveCreateModelParts(model_part, "sub_model_part")
        self.assertTrue(model_part.HasSubModelPart("sub_model_part"))

        self.assertFalse(model_part.GetSubModelPart("sub_model_part").HasSubModelPart("subsub"))
        model_part_utilities.RecursiveCreateModelParts(model_part, "sub_model_part.subsub")
        self.assertTrue(model_part.GetSubModelPart("sub_model_part").HasSubModelPart("subsub"))

    def test_CreateModelPartsFromCouplingDataSettings(self):
        model_part_utilities.CreateModelPartsFromCouplingDataSettings(self.data_settings, self.model, "dummy_solver_name")

        # make sure the original settings are not modified
        self.assertFalse(self.data_settings["disp"].Has("location"))

        self.assertTrue(self.model.HasModelPart("Structure"))
        self.assertTrue(self.model.HasModelPart("Structure.Interface"))
        self.assertTrue(self.model.HasModelPart("Structure.GENERIC_FSI.subinterface"))
        self.assertTrue(self.model.HasModelPart("aux_model_part"))


if __name__ == '__main__':
    KratosUnittest.main()
