from __future__ import print_function, absolute_import, division

from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
import KratosMultiphysics.KratosUnittest as UnitTest

class TestGiDIOGaussPoints(UnitTest.TestCase):

    def setUp(self):
        self.modelPart = self.setModelPart()

    def tearDown(self):
        pass

    def setModelPart(self):
        modelPart = ModelPart("Test ModelPart")

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

        elements = list()
        elements.append( modelPart.CreateNewElement("VMS3D",1,[1,2,4,5],properties) )
        elements.append( modelPart.CreateNewElement("VMS3D",2,[2,3,4,5],properties) )
        
        modelPart.SetBufferSize(2)

        return modelPart

    def disableSome(self):
        for elem in self.modelPart.Elements:
            if elem.Id % 2 == 0:
                elem.Set(ACTIVE,False)

    def writeTestFile(self):

        gid_io = GidIO(
            "test_git_io_gauss_points",
            GiDPostMode.GiD_PostBinary,
            MultiFileFlag.SingleFile,
            WriteDeformedMeshFlag.WriteUndeformed,
            WriteConditionsFlag.WriteConditions)

        gid_io.InitializeMesh(0)
        gid_io.WriteMesh(self.modelPart.GetMesh())
        gid_io.FinalizeMesh()

        gid_io.InitializeResults(0.0, self.modelPart.GetMesh())
        gid_io.WriteNodalResults(VELOCITY, self.modelPart.Nodes, 0.0, 0)
        gid_io.PrintOnGaussPoints(VORTICITY, self.modelPart, 0.0)
        gid_io.FinalizeResults()

if __name__ == '__main__':
    test = TestGiDIOGaussPoints()
    test.setUp()
    test.disableSome()
    test.writeTestFile()
    test.tearDown()