import os
import math
import json
import KratosMultiphysics
from KratosMultiphysics import ParallelEnvironment, IsDistributedRun
import KratosMultiphysics.MeshingApplication
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestMPIParMmg(KratosUnittest.TestCase):

    @KratosUnittest.skipUnless(IsDistributedRun() and ParallelEnvironment.GetDefaultSize() <= 4,  "Test designed to be run with max. 4 ranks.")
    def test_mpi_sphere(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # We create the model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)

        # We add the variables needed
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)


        # We import the model main_model_part
        file_path = GetFilePath("/parmmg_eulerian_test/background_mesh_sphere")
        ReadModelPart(file_path, main_model_part)

        communicator = main_model_part.GetCommunicator().GetDataCommunicator()

        for node in main_model_part.Nodes:
            distance = math.sqrt(node.X**2+node.Y**2+node.Z**2) - 1.0/2.0
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance)

        # Setting some flags on submodelparts entities arbitrarily
        for cond in main_model_part.GetSubModelPart("SurfaceLoad3D_Load_on_surfaces_Auto3").Conditions:
            cond.Set(KratosMultiphysics.STRUCTURE)
        for elem in main_model_part.GetSubModelPart("Parts_Solid_Solid_Auto1").Elements:
            elem.Set(KratosMultiphysics.VISITED)

        ##COMPUTE DISTANCE GRADIENT AND NODAL_H
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)

        local_gradient.Execute()
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()

        ##COMPUTE LEVEL SET METRIC
        metric_parameters = KratosMultiphysics.Parameters("""
        {
            "minimal_size"                             : 0.5,
            "sizing_parameters": {
                "reference_variable_name"               : "DISTANCE",
                "boundary_layer_max_distance"           : 2.0,
                "interpolation"                         : "constant"
            },
            "enforce_current"                      : false,
            "anisotropy_remeshing"                 : false
        }
        """)
        metric_process = KratosMultiphysics.MeshingApplication.ComputeLevelSetSolMetricProcess3D(
            main_model_part,
            KratosMultiphysics.DISTANCE_GRADIENT,
            metric_parameters)
        metric_process.Execute()

        ##PERFORM REMESHING
        pmmg_parameters = KratosMultiphysics.Parameters("""
        {
            "filename"                         : "output",
            "save_external_files"              : true,
            "save_colors_files"                : true,
            "initialize_entities"              : false,
            "preserve_flags"                   : true,
            "echo_level"                       : 0,
            "advanced_parameters" : {
                "number_of_iterations"         : 4
            }
        }
        """)
        pmmg_parameters["filename"].SetString(GetFilePath(pmmg_parameters["filename"].GetString()))
        pmmg_process = KratosMultiphysics.MeshingApplication.ParMmgProcess3D(main_model_part.GetRootModelPart(), pmmg_parameters)
        pmmg_process.Execute()

        reference_file_name = GetFilePath("parmmg_eulerian_test/cond_ref_map.json")
        result_file_name = GetFilePath("output_step=0_"+str(communicator.Rank())+".cond.ref.json")
        self._CompareColorFiles(reference_file_name, result_file_name)

        reference_file_name = GetFilePath("parmmg_eulerian_test/elem_ref_map.json")
        result_file_name = GetFilePath("output_step=0_"+str(communicator.Rank())+".elem.ref.json")
        self._CompareColorFiles(reference_file_name, result_file_name)

        result_dict_file_name=GetFilePath("parmmg_eulerian_test/reference_parmmg_spehere_mdpa_hierarchy.json")
        with open(result_dict_file_name, 'r') as f:
            reference_hierarchy = json.load(f)
        self.CheckModelPartHierarchie(main_model_part, reference_hierarchy[str(communicator.Size())])

        # Check flags are correctly set on the corresponding submodelparts only
        for cond in main_model_part.Conditions:
            if cond in main_model_part.GetSubModelPart("SurfaceLoad3D_Load_on_surfaces_Auto3").Conditions:
                self.assertTrue(cond.Is(KratosMultiphysics.STRUCTURE))
            else:
                self.assertTrue(cond.IsNot(KratosMultiphysics.STRUCTURE))

        for elem in main_model_part.Elements:
            if elem in main_model_part.GetSubModelPart("Parts_Solid_Solid_Auto1").Elements:
                self.assertTrue(elem.Is(KratosMultiphysics.VISITED))
            else:
                self.assertTrue(elem.IsNot(KratosMultiphysics.VISITED))

        for file_name in os.listdir(GetFilePath("")):
            if file_name.endswith(".json") or file_name.endswith(".mdpa") or file_name.endswith(".mesh") or  file_name.endswith(".sol"):
                kratos_utilities.DeleteFileIfExisting(GetFilePath(file_name))
        kratos_utilities.DeleteTimeFiles(os.getcwd())

    def _CompareColorFiles(self, ref_dict_filename, result_dict_file_name):

        with open(ref_dict_filename, 'r') as f:
            reference_values = json.load(f)

        with open(result_dict_file_name, 'r') as f:
            result_values = json.load(f)

        self.assertEqual(len(reference_values.keys()), len(result_values.keys()))
        for key_ref, key_result in zip(reference_values.keys(), result_values.keys()):
            self.assertEqual(reference_values[key_ref], result_values[key_result])

    def _CheckModelPart(self, ref_model_part, result_model_part):
        self.assertEqual(ref_model_part.NumberOfNodes(), result_model_part.NumberOfNodes())
        self.assertEqual(ref_model_part.NumberOfElements(), result_model_part.NumberOfElements())
        self.assertEqual(ref_model_part.NumberOfConditions(), result_model_part.NumberOfConditions())

    def CheckModelPartHierarchie(self, model_part, hierarchie):
        """Checking if the hierarchie of a ModelPart matches the expected one
        This is intended to check larger models, where it is not feasible
        save large mdpa-files as references
        the hierarchie is a dict with the structure of the ModelPart. E.g.:
        {
            "name_model_part" : {
                "nodes": 15
                "elements": 11
                "conditions": 5
                "properties": 2,
                "sub_model_parts" : {
                    "domain" : {
                        "nodes": 15,
                        "elements" : 11,
                        "properties" :1
                        "sub_model_parts" : {
                            "sub_domain" : {
                                "nodes" : 3,
                                "elements" : 2
                            }
                        }
                    },
                    "boundary" : {
                        "nodes": 6
                        "conditions" : 5,
                        "properties" : 1
                    }
                    }
                }
            }
        }
        """
        def CheckModelPartHierarchieNumbers(smp, smp_hierarchie):
            comm = smp.GetCommunicator().GetDataCommunicator()
            local_number_of_nodes = smp.GetCommunicator().LocalMesh().NumberOfNodes()
            local_number_of_elem = smp.GetCommunicator().LocalMesh().NumberOfElements()
            local_number_of_cond = smp.GetCommunicator().LocalMesh().NumberOfConditions()
            local_number_of_prop = smp.GetCommunicator().LocalMesh().NumberOfProperties()

            exp_num = smp_hierarchie.get("nodes", 0)
            self.assertEqual(comm.SumAll(local_number_of_nodes), exp_num, msg='ModelPart "{}" is expected to have {} nodes but has {}'.format(smp.FullName(), exp_num, smp.NumberOfNodes()))

            exp_num = smp_hierarchie.get("elements", 0)
            self.assertEqual(comm.SumAll(local_number_of_elem), exp_num, msg='ModelPart "{}" is expected to have {} elements but has {}'.format(smp.FullName(), exp_num, smp.NumberOfElements()))

            exp_num = smp_hierarchie.get("conditions", 0)
            self.assertEqual(comm.SumAll(local_number_of_cond), exp_num, msg='ModelPart "{}" is expected to have {} conditions but has {}'.format(smp.FullName(), exp_num, smp.NumberOfConditions()))

            exp_num = smp_hierarchie.get("properties", 0)
            self.assertEqual(comm.SumAll(local_number_of_prop), exp_num, msg='ModelPart "{}" is expected to have {} properties but has {}'.format(smp.FullName(), exp_num, smp.NumberOfProperties()))

            if "sub_model_parts" in smp_hierarchie:
                smp_hierarchie = smp_hierarchie["sub_model_parts"]
                for name_smp in smp_hierarchie:
                    self.assertTrue(smp.HasSubModelPart(name_smp), msg='ModelPart "{}" does not have SubModelPart with name "{}"'.format(smp.FullName(), name_smp))
                    CheckModelPartHierarchieNumbers(smp.GetSubModelPart(name_smp), smp_hierarchie[name_smp])

        # check name of MainModelPart
        self.assertEqual(len(hierarchie), 1)
        name_main_model_part = hierarchie.__iter__().__next__()
        self.assertEqual(model_part.Name, name_main_model_part)

        CheckModelPartHierarchieNumbers(model_part, hierarchie[name_main_model_part])

if __name__ == '__main__':
    KratosUnittest.main()
