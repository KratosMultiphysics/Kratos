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

class ContactModels3DTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):
    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_contact_models_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def Finalize(self):
        self.CompareResults()
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), "")
        super().Finalize()

    def CompareResults(self):
        node = 335
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
                "force": [-756.509124672903454, -1168.369818390392993, -145.804749769301509],
                "torque": [-119.790077944505612, -26.819218395491443, -7.595296848374547]
                }
        elif materials_filename == "MaterialsDEM_Linear_Classic.json":
            return {
                "force": [-47.679531905699207, -41.377312890971552, 2500.104846833379270],
                "torque": [-85.559759618387559, 18.107378192547465, 5.168848884958166]
                }
        elif materials_filename == "MaterialsDEM_Linear_Viscous_Coulomb.json":
            return {
                "force": [-56.794788109531510, -53.047528652410165, 1408.214206348576681],
                "torque": [-29.206861740144873, 11.506795739843547, -22.952174514486565]
                }
        elif materials_filename == "MaterialsDEM_Hertz_Viscous_Coulomb.json":
            return {
                "force": [-29.997387927285470, -59.838478020708408, 799.969409659972371],
                "torque": [-14.527599771313696, 7.441196776833555, -4.917194429095725]
                }
        else:
            raise ValueError(f"Unknown materials file name: {materials_filename}")

class Test3DContactModels(KratosUnittest.TestCase):
    def setUp(self):
        pass

    @classmethod
    def test_DEM3D_ContactModels(self):
        materials_files_name = [f for f in os.listdir(os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_contact_models_tests_files")) if f.startswith("MaterialsDEM_") and f.endswith(".json")]
        parameters_file_name = os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM3D_contact_models_tests_files"), "ProjectParametersDEM.json")

        for materials_file in materials_files_name:
            with open(parameters_file_name, "r") as file:
                parameters = Kratos.Parameters(file.read())
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString(materials_file)
            with open(parameters_file_name, "w") as file:
                file.write(parameters.PrettyPrintJsonString())
            model = Kratos.Model()
            auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(ContactModels3DTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()