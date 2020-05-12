from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
from KratosMultiphysics.FluidDynamicsApplication.apply_mass_conservation_check_process import ApplyMassConservationCheckProcess
import re

if KratosMultiphysics.DataCommunicator.GetDefault().IsDistributed():
    import KratosMultiphysics.mpi.distributed_import_model_part_utility as distributed_import_model_part_utility

def GetFilePath(fileName):
    return KratosMultiphysics.os.path.dirname(KratosMultiphysics.os.path.realpath(__file__)) + "/" + fileName

class MassConservationUtility(KratosUnittest.TestCase):
    def setUp(self):
        self.comm = KratosMultiphysics.DataCommunicator.GetDefault()
        self.size = self.comm.Size()
        self.rank = self.comm.Rank()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.comm.Barrier()
            if self.rank == 0:
                kratos_utils.DeleteFileIfExisting(self.file_name + ".time")
            kratos_utils.DeleteFileIfExisting(self.file_name + "_" + str(self.rank) + ".mdpa")
            self.comm.Barrier()

    @staticmethod
    def _SetParameters(name, use_memory):
        parameters = """{
            "echo_level" : 0,
            "model_import_settings" : {
                "input_type" : "mdpa",
                "input_filename" :\""""+ name +"""\",
                "partition_in_memory" : """ + use_memory + """
            }
        }"""

        return parameters

    def _CreateModelPart(self):
        self.parameters = self._SetParameters(self.file_name, use_memory="false")
        self.model = KratosMultiphysics.Model()
        model_part = self.model.CreateModelPart("ModelPart")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        if self.comm.IsDistributed():
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

        model_part.SetBufferSize(3)

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            import_settings = KratosMultiphysics.Parameters(self.parameters)
            self._ImportModelPart(model_part, import_settings)

        return model_part

    def _ImportModelPart(self, model_part, import_settings):
        if self.comm.IsDistributed():
            distributed_model_part_importer = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, import_settings)
            distributed_model_part_importer.ImportModelPart()
            distributed_model_part_importer.CreateCommunicators()
        else:
            import_flags = KratosMultiphysics.ModelPartIO.READ
            KratosMultiphysics.ModelPartIO(import_settings['model_import_settings']['input_filename'].GetString(), import_flags).ReadModelPart(model_part)

    def _SetInletAndOutlet(self, model_part):
        inlet_conds = [5,6,17,18,29,30,41,42]
        outlet_conds = [11,12,23,24,35,36,47,48]
        for cond in model_part.Conditions:
            if cond.Id in inlet_conds:
                cond.Set(KratosMultiphysics.INLET)
            if cond.Id in outlet_conds:
                cond.Set(KratosMultiphysics.OUTLET)



    def test_Initialize(self):
        self.work_folder = "auxiliary_files"
        self.file_name = "cube_tetrahedra_elements_coarse"
        model_part = self._CreateModelPart()
        for node in model_part.GetCommunicator().LocalMesh().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.Y - 0.5)
        settings = KratosMultiphysics.Parameters()
        mass_conservation_process = KratosFluid.MassConservationUtility(model_part, settings)
        message = mass_conservation_process.Initialize()
        filtered_results = re.findall(r"[-+]?\d*\.\d+|\d+", message)
        self.assertAlmostEqual(float(filtered_results[0]), 0.5)
        self.assertAlmostEqual(float(filtered_results[1]), 0.5)
        self.assertAlmostEqual(float(filtered_results[2]), 1.0)

    def test_ComputeBalancedVolume(self):
        self.work_folder = "auxiliary_files"
        self.file_name = "cube_tetrahedra_elements_coarse"
        model_part = self._CreateModelPart()
        for node in model_part.GetCommunicator().LocalMesh().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.Y - 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 1.0 - 0.5*node.X)
        self._SetInletAndOutlet(model_part)
        settings = KratosMultiphysics.Parameters()
        mass_conservation_process = KratosFluid.MassConservationUtility(model_part, settings)
        model_part.CloneTimeStep(0.1)
        model_part.ProcessInfo[KratosMultiphysics.TIME_CRT] = 0.1
        mass_conservation_process.Initialize()
        balance_message = mass_conservation_process.ComputeBalancedVolume()
        filtered_results = re.findall(r"[-+]?\d*\.\d+|\d+", balance_message)
        current_time = float(filtered_results[0])
        neg_vol = float(filtered_results[1])
        pos_vol = float(filtered_results[2])
        volume_error = float(filtered_results[3])
        net_inflow_inlet = float(filtered_results[4])
        net_inflow_outlet = float(filtered_results[5])
        inlet_area = float(filtered_results[6])
        self.assertAlmostEqual(current_time, 0.1)
        self.assertAlmostEqual(neg_vol, 0.5)
        self.assertAlmostEqual(pos_vol, 0.5)
        self.assertAlmostEqual(volume_error, 0.025)
        self.assertAlmostEqual(net_inflow_inlet, 0.5)
        self.assertAlmostEqual(net_inflow_outlet, -0.25)
        self.assertAlmostEqual(inlet_area, 1.0)

    def test_OrthogonalFlowIntoAir(self):
        self.work_folder = "auxiliary_files"
        self.file_name = "cube_tetrahedra_elements_coarse"
        model_part = self._CreateModelPart()
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.Y - 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 2.0) #This term should not affect
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 3.0 - node.Z)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_Y, 1.0)

        settings = KratosMultiphysics.Parameters()
        mass_conservation_process = KratosFluid.MassConservationUtility(model_part, settings)
        orthogonal_flow = mass_conservation_process.OrthogonalFlowIntoAir(1.0)
        self.assertAlmostEqual(orthogonal_flow,1.5)

    def test_NoOrthogonalFlowIntoAir(self):
        self.work_folder = "auxiliary_files"
        self.file_name = "cube_tetrahedra_elements_coarse"
        model_part = self._CreateModelPart()
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.Y - 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 2.0) #This term should not affect
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, -3.0 + node.Z)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_Y, -1.0)

        settings = KratosMultiphysics.Parameters()
        mass_conservation_process = KratosFluid.MassConservationUtility(model_part, settings)
        orthogonal_flow = mass_conservation_process.OrthogonalFlowIntoAir(1.0)
        self.assertAlmostEqual(orthogonal_flow,0.0)

    def test_ReverseOrthogonalFlowIntoAir(self):
        self.work_folder = "auxiliary_files"
        self.file_name = "cube_tetrahedra_elements_coarse"
        model_part = self._CreateModelPart()
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.Y - 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 2.0) #This term should not affect
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, -3.0 + node.Z)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY_Y, -1.0)

        settings = KratosMultiphysics.Parameters()
        mass_conservation_process = KratosFluid.MassConservationUtility(model_part, settings)
        orthogonal_flow = mass_conservation_process.OrthogonalFlowIntoAir(-1.0)
        self.assertAlmostEqual(orthogonal_flow,1.5)

    def test_ComputeDtForConvection(self):
        self.work_folder = "auxiliary_files"
        self.file_name = "cube_tetrahedra_elements_coarse"
        model_part = self._CreateModelPart()
        for node in model_part.GetCommunicator().LocalMesh().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.X - 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 1.0 - node.X)

        self._SetInletAndOutlet(model_part)
        settings = KratosMultiphysics.Parameters()
        mass_conservation_process = KratosFluid.MassConservationUtility(model_part, settings)
        model_part.CloneTimeStep(0.1)
        mass_conservation_process.Initialize()
        mass_conservation_process.ComputeBalancedVolume()
        dt = mass_conservation_process.ComputeDtForConvection()
        self.assertAlmostEqual(dt, 0.2)

    def test_CorrectMass(self):
        self.work_folder = "auxiliary_files"
        self.file_name = "cube_tetrahedra_elements_coarse"
        model_part = self._CreateModelPart()
        for node in model_part.GetCommunicator().LocalMesh().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.X - 0.5)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 1.0 - node.Y)

        self._SetInletAndOutlet(model_part)
        settings = KratosMultiphysics.Parameters("""{
            "model_part_name" : "ModelPart"
        }""")
        mass_conservation_process = ApplyMassConservationCheckProcess(model_part.GetModel(), settings)
        mass_conservation_process.ExecuteInitialize()
        model_part.CloneTimeStep(0.1)
        mass_conservation_process.ExecuteFinalizeSolutionStep()
        results = {1: -0.5938516359222957,2: -0.09699581511190579,3: 0.4014031390113931,4: -0.5517656205702922,
            5: -0.05119945641905187, 6: 0.4492730587186531,7: -0.5033021519665317,8: -1e-07,9: 0.4991147948829323,
            10: -0.5917444339533421,11: -0.09536505385930948,12: 0.4020325332379011,13: -0.5508310379528296,
            14: -0.05344725574988601,15: 0.4487503962770498,16: -0.5042584883219816,17: -0.00027154466803484535,
            18: 0.5,19: -0.5867198935872885,20: -0.09537629393576269,21: 0.40337851653918727,22: -0.5466292414006754,
            23: -0.05774088594966558,24: 0.4494118475239448,25: -0.5256186306939375,26: -1e-07,27: 0.49984147049824706}
        for node in model_part.GetCommunicator().LocalMesh().Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE), results.get(node.Id))
        self.assertAlmostEqual(model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME], 0.1)


if __name__ == '__main__':
    KratosUnittest.main()
