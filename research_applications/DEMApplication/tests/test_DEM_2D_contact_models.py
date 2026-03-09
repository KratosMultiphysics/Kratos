import os
import shutil
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication as DEMApplication
import KratosMultiphysics.DEMApplication.DEM_analysis_stage
import auxiliary_functions_for_tests

# To add a new contact model test, create a "MaterialsDEM_Contact_Model_Name.json" file that uses the target contact model to be tested,
# and update the reference results in the GetReferenceResults function accordingly.

class ContactModels2DTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):
    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_contact_models_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        self.CompareResults()
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), "")
        super().Finalize()

    def CompareResults(self):
        node = 16
        tol = 1e-2
        ref = self.GetReferenceResults()
        f_ref = ref["force"]
        t_ref = ref["torque"]
        f = self.spheres_model_part.GetNode(node).GetSolutionStepValue(DEMApplication.CONTACT_FORCES)
        t = self.spheres_model_part.GetNode(node).GetSolutionStepValue(DEMApplication.PARTICLE_MOMENT)
        #print(self.DEM_parameters["solver_settings"]["material_import_settings"]["materials_filename"].GetString())
        #print("Force:  [{:.15f}, {:.15f}, {:.15f}]".format(f[0],f[1],f[2]))
        #print("Moment: [{:.15f}, {:.15f}, {:.15f}]".format(t[0],t[1],t[2]))
        self.assertAlmostEqual(f[0], f_ref[0], delta=tol)
        self.assertAlmostEqual(f[1], f_ref[1], delta=tol)
        self.assertAlmostEqual(f[2], f_ref[2], delta=tol)
        self.assertAlmostEqual(t[0], t_ref[0], delta=tol)
        self.assertAlmostEqual(t[1], t_ref[1], delta=tol)
        self.assertAlmostEqual(t[2], t_ref[2], delta=tol)

    def GetReferenceResults(self):
        materials_filename = self.DEM_parameters["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
        if materials_filename == "MaterialsDEM_Linear_Simple_Coulomb.json":
            return {
                "force": [-1789.794623844300077, -2435.638073447114039, 0.000000000000000],
                "torque": [0.000000000000000, 0.000000000000000, 96.918194075810305]
                }
        elif materials_filename == "MaterialsDEM_Linear_Viscous_Coulomb2D.json":
            return {
                "force": [3522.406660549475419, -254.294001706025028, 0.000000000000000],
                "torque": [0.000000000000000, 0.000000000000000, 489.529890369310749]
                }
        elif materials_filename == "MaterialsDEM_Hertz_Viscous_Coulomb2D.json":
            return {
                "force": [3027.754322894346842, -270.827205739553392, 0.000000000000000],
                "torque": [0.000000000000000, 0.000000000000000, 368.864261863613478]
                }
        else:
            raise ValueError(f"Unknown materials file name: {materials_filename}")

class Test2DContactModels(KratosUnittest.TestCase):
    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_ContactModels(self):
        materials_files_name = [f for f in os.listdir(os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_contact_models_tests_files")) if f.startswith("MaterialsDEM_") and f.endswith(".json")]
        parameters_file_name = os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_contact_models_tests_files"), "ProjectParametersDEM.json")

        for materials_file in materials_files_name:
            with open(parameters_file_name, "r") as file:
                parameters = Kratos.Parameters(file.read())
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString(materials_file)
            with open(parameters_file_name, "w") as file:
                file.write(parameters.PrettyPrintJsonString())
            model = Kratos.Model()
            auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(ContactModels2DTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()