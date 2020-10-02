import os
import sys
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.mpi import distributed_import_model_part_utility
import json

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestMPIProcesses(KratosUnittest.TestCase):

    def _ReadSerialModelPart(self,model_part, mdpa_file_name):
        import_flags = KratosMultiphysics.ModelPartIO.READ | KratosMultiphysics.ModelPartIO.SKIP_TIMER
        KratosMultiphysics.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

    def _ReadDistributedModelPart(self,model_part, mdpa_file_name):

        importer_settings = KratosMultiphysics.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + mdpa_file_name + """\",
                "partition_in_memory" : true
            },
            "echo_level" : 0
        }""")

        model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, importer_settings)
        model_part_import_util.ImportModelPart()
        model_part_import_util.CreateCommunicators()

    def _ReadModelPart(self, input_filename, mp):
        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()

        if kratos_comm.IsDistributed():
            self._ReadDistributedModelPart(mp, input_filename)
        else:
            self._ReadSerialModelPart(mp, input_filename)

    def testComputeNodalGradientProcess(self):

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()

        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("main_model_part")
        main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]=2

        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

        ## Serial partition of the original .mdpa file
        input_filename = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_mpi_serializer")
        self._ReadModelPart(input_filename, main_model_part)

        ParallelFillCommunicator = KratosMultiphysics.mpi.ParallelFillCommunicator(main_model_part)
        ParallelFillCommunicator.Execute()

        for node in main_model_part.Nodes:
            distance = node.X**2+node.Y**2+node.Z**2 - 1
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,distance)

        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        ## Read reference
        file_name = "auxiliar_files_for_python_unittest/reference_files/test_compute_nodal_gradient_process_results.json"
        reference_file_name = GetFilePath(file_name)
        reference_values = json.load(open(reference_file_name, 'r'))

        for node in main_model_part.Nodes:
            distance_gradient = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for gradient_i, gradient_i_ref in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(gradient_i, gradient_i_ref)

    def testComputeNonHistoricalNodalGradientProcess(self):

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()

        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("main_model_part")
        main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]=2

        ## Add variables to the model part
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        ## Serial partition of the original .mdpa file
        input_filename = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_mpi_serializer")
        self._ReadModelPart(input_filename, main_model_part)

        ParallelFillCommunicator = KratosMultiphysics.mpi.ParallelFillCommunicator(main_model_part)
        ParallelFillCommunicator.Execute()

        for node in main_model_part.Nodes:
            distance = node.X**2+node.Y**2+node.Z**2 - 1
            node.SetValue(KratosMultiphysics.DISTANCE,distance)

        non_historical_gradient_variable = True
        local_gradient = KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(main_model_part,
        KratosMultiphysics.DISTANCE,
        KratosMultiphysics.DISTANCE_GRADIENT,
        KratosMultiphysics.NODAL_AREA,
        non_historical_gradient_variable)
        local_gradient.Execute()

        ## Read reference
        file_name = "auxiliar_files_for_python_unittest/reference_files/test_compute_nodal_gradient_process_results.json"
        reference_file_name = GetFilePath(file_name)
        reference_values = json.load(open(reference_file_name, 'r'))

        for node in main_model_part.Nodes:
            distance_gradient = node.GetValue(KratosMultiphysics.DISTANCE_GRADIENT)
            ref_gradient = reference_values[str(node.Id)]
            for gradient_i, gradient_i_ref in zip(distance_gradient, ref_gradient):
                self.assertAlmostEqual(gradient_i, gradient_i_ref)

if __name__ == '__main__':
    KratosUnittest.main()
