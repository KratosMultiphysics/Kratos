from __future__ import print_function, absolute_import, division

from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as UnitTest

import KratosMultiphysics.kratos_utilities as kratos_utils

try:
    from KratosMultiphysics.FluidDynamicsApplication import *
    have_fluid_dynamics = True
except ImportError:
    have_fluid_dynamics = False

import filecmp
import os

class WorkFolderScope(object):
    '''Auxiliary class to define a work folder for the tests.'''
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(have_fluid_dynamics,"Missing required application: FluidDynamicsApplication")
class TestGiDIOGaussPoints(UnitTest.TestCase):
    '''Tests related to GiD I/O Gauss point results printing.'''

    def setUp(self):
        self.setModelPart()
        self.workFolder = "auxiliar_files_for_python_unittest/gid_io"

    def tearDown(self):
        with WorkFolderScope(self.workFolder):
            for suffix in ['_0.post.res', '_0.post.msh']:
                kratos_utils.DeleteFileIfExisting(self.output_file_name+suffix)


    def setModelPart(self):
        self.model = Model()
        modelPart = self.model.CreateModelPart("Test ModelPart")

        modelPart.AddNodalSolutionStepVariable(DISTANCE)
        modelPart.AddNodalSolutionStepVariable(VELOCITY)

        nodes = list()
        nodes.append( modelPart.CreateNewNode(1, 0.0, 0.0, 0.0) )
        nodes.append( modelPart.CreateNewNode(2, 1.0, 0.0, 0.0) )
        nodes.append( modelPart.CreateNewNode(3, 1.0, 1.0, 0.0) )
        nodes.append( modelPart.CreateNewNode(4, 0.0, 1.0, 0.0) )
        nodes.append( modelPart.CreateNewNode(5, 0.5, 0.5, 1.0) )

        for node in nodes:
            rx = node.X - 0.5
            rz = node.Z - 0.5
            r = (rx**2 + rz**2)**0.5
            vel = Array3()
            vel[0] = - rz/r
            vel[1] = 0.0
            vel[2] = rx/r
            node.SetSolutionStepValue(VELOCITY,0,vel)
            node.SetSolutionStepValue(DISTANCE,0,r)

        properties = modelPart.GetProperties()[0]

        modelPart.CreateNewElement("VMS3D",1,[1,2,4,5],properties)
        modelPart.CreateNewElement("VMS3D",2,[2,3,4,5],properties)

        modelPart.CreateNewCondition("MonolithicWallCondition3D",1,[1,5,4],properties)
        modelPart.CreateNewCondition("MonolithicWallCondition3D",2,[1,2,5],properties)
        modelPart.CreateNewCondition("MonolithicWallCondition3D",3,[2,3,5],properties)
        modelPart.CreateNewCondition("MonolithicWallCondition3D",4,[3,4,5],properties)

        modelPart.SetBufferSize(2)

        self.modelPart = modelPart

    def deactivateSome(self):
        for elem in self.modelPart.Elements:
            if elem.Id % 2 == 0:
                elem.Set(ACTIVE,False)

        for cond in self.modelPart.Conditions:
            if cond.Id % 2 == 0:
                cond.Set(ACTIVE,False)

    def initializeOutputFile(self):
        self.gid_io = GidIO(
            self.output_file_name,
            self.post_mode,
            MultiFileFlag.SingleFile,
            WriteDeformedMeshFlag.WriteUndeformed,
            WriteConditionsFlag.WriteConditions)

        self.gid_io.InitializeMesh(0)
        self.gid_io.WriteMesh(self.modelPart.GetMesh())
        self.gid_io.FinalizeMesh()

        self.gid_io.InitializeResults(0.0, self.modelPart.GetMesh())

    def writeResults(self,label):
        self.gid_io.WriteNodalResults(VELOCITY, self.modelPart.Nodes, label, 0)
        self.gid_io.PrintOnGaussPoints(VORTICITY, self.modelPart, label)
        self.gid_io.PrintOnGaussPoints(NORMAL, self.modelPart, label)
        self.gid_io.PrintFlagsOnGaussPoints(ACTIVE, "ACTIVE", self.modelPart, label)


    def finalizeOutputFile(self):
        self.gid_io.FinalizeResults()

    def outputMatchesReferenceSolution(self):
        msh_file_matches = filecmp.cmp(self.reference_file_name+'_0.post.msh',self.output_file_name+'_0.post.msh')
        res_file_matches = filecmp.cmp(self.reference_file_name+'_0.post.res',self.output_file_name+'_0.post.res')
        return msh_file_matches and res_file_matches


    def test_write_active_only(self):

        self.post_mode = GiDPostMode.GiD_PostAscii
        self.output_file_name = "test_gid_io_gp_active_only"
        self.reference_file_name = "ref_gid_io_gp_active_only"

        self.deactivateSome()

        with WorkFolderScope(self.workFolder):
            self.initializeOutputFile()
            self.writeResults(0.0)
            self.finalizeOutputFile()

            self.assertTrue(self.outputMatchesReferenceSolution())

    def test_write_dynamic_deactivation(self):

        self.post_mode = GiDPostMode.GiD_PostAscii
        self.output_file_name = "test_gid_io_gp_dynamic_deactivation"
        self.reference_file_name = "ref_gid_io_gp_dynamic_deactivation"

        with WorkFolderScope(self.workFolder):
            self.initializeOutputFile()
            self.writeResults(0.0)

            self.deactivateSome()

            self.writeResults(1.0)
            self.finalizeOutputFile()

            self.assertTrue(self.outputMatchesReferenceSolution())


if __name__ == '__main__':
    test = TestGiDIOGaussPoints()
    test.setUp()
    test.test_write_active_only()
    test.tearDown()
    test.setUp()
    test.test_write_dynamic_deactivation()
    test.tearDown()
