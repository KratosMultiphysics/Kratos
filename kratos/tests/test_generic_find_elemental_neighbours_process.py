import os
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.testing.utilities import ReadModelPart
import KratosMultiphysics.kratos_utilities as kratos_utils

class TestGenericFindElementalNeighboursProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        self.size = self.comm.Size()
        self.rank = self.comm.Rank()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.comm.Barrier()
            if self.rank == 0:
                kratos_utils.DeleteDirectoryIfExisting(self.file_name + "_partitioned")
            self.comm.Barrier()

    def testGenericFindElementalNeighboursProcess(self):
        self.work_folder = "auxiliar_files_for_python_unittest/mdpa_files/"
        self.file_name = "small_cube"
        model = self._ImportModelPart(GetFilePath(self.work_folder+ "/"+self.file_name)) #MdpaTriangleDummyModel()
        model_part = model["main"]
        proc = KratosMultiphysics.GenericFindElementalNeighboursProcess(model_part)
        proc.ExecuteInitialize()

        #expected results. faces which have a null pointer (no neigh elem in face)
        null_pointers_faces_dict = {
            1 : [False,False,True,True],
            2 : [True,False,False,True],
            3 : [False,False,True,True],
            4 : [True,False,False,True],
            5 : [True,False,False,True],
            6 : [False,False,True,True]
        }

        for elem in model_part.Elements:
            elem_id = elem.Id
            neigh_elems_pointers_str = proc.HelperTestFunctionReturnPointersInStr(elem_id)
            this_elem_null_pointers_expected = null_pointers_faces_dict.get(elem_id)
            self.assertEqual(len(neigh_elems_pointers_str),4)
            for face in range(len(neigh_elems_pointers_str)):
                is_null_pointer = neigh_elems_pointers_str[face]=="0000000000000000"
                self.assertEqual(is_null_pointer,this_elem_null_pointers_expected[face])
                
        

    def _ImportModelPart(self,model_path):
        Model = KratosMultiphysics.Model()
        model_part = Model.CreateModelPart("main")
        if self.comm.IsDistributed():
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        ReadModelPart(model_path,model_part)
        return Model

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(__file__), fileName)

if __name__ == '__main__':
    KratosUnittest.main()
