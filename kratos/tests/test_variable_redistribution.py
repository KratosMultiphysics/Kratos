from KratosMultiphysics import *
try:
    from KratosMultiphysics.FSIApplication import *
    have_fsi = True
except ImportError:
    have_fsi = False

import KratosMultiphysics.KratosUnittest as UnitTest

from os import remove

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(have_fsi,"Missing required application: FSIApplication")
class VariableRedistributionTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.input_file = "redistribution_test"
        self.work_folder = "RedistributionTest"

        self.redistribution_iterations = 100
        self.redistribution_tolerance = 1e-7

        self.check_tolerance = 1e-6
        self.print_output = False

        self.mapped_function = None

    def tearDown(self):
        self.deleteOutFile(self.input_file+'.time')

    def deleteOutFile(self,filename):
        with WorkFolderScope(self.work_folder):
            try:
                remove(filename)
            except FileNotFoundError as e:
                pass

    def testLinearFunction(self):
        def Flag1Check(node):
            return node.GetSolutionStepValue(FLAG_VARIABLE) == 1.0

        def ReferenceSolution(node,variable):
            node.SetSolutionStepValue(variable,node.X)

        self.RunTestCase(
            Flag1Check,
            ReferenceSolution,
            PRESSURE,
            NODAL_PAUX,
            TEMPERATURE
        )

    def testSharpCorners(self):

        def FlagDefinedCheck(node):
            return node.GetSolutionStepValue(FLAG_VARIABLE) > 0.0

        def ReferenceSolution(node,variable):
            node.SetSolutionStepValue(variable, 10.*node.X + node.Y)


        self.RunTestCase(
            FlagDefinedCheck,
            ReferenceSolution,
            PRESSURE,
            NODAL_PAUX,
            TEMPERATURE
        )

    def testQuadratic(self):

        def FlagDefinedCheck(node):
            return node.GetSolutionStepValue(FLAG_VARIABLE) > 0.0

        def ReferenceSolution(node,variable):
            node.SetSolutionStepValue(variable, node.X*node.X + node.Y*node.Z)


        self.RunTestCase(
            FlagDefinedCheck,
            ReferenceSolution,
            PRESSURE,
            NODAL_PAUX,
            TEMPERATURE
        )

    def testVector(self):

        def FlagDefinedCheck(node):
            return node.GetSolutionStepValue(FLAG_VARIABLE) > 0.0

        def ReferenceSolution(node,variable):
            value = Array3()
            value[0] = node.Y
            value[1] = 10.*node.X + node.Z
            value[2] = 500
            node.SetSolutionStepValue(variable, value)


        self.RunTestCase(
            FlagDefinedCheck,
            ReferenceSolution,
            VELOCITY,
            VORTICITY,
            ACCELERATION
        )


    def RunTestCase(self,inteface_check,set_reference,reference_variable,intermediate_variable,result_variable):
        with WorkFolderScope(self.work_folder):

            self.SetUpProblem()

            for node in self.model_part.Nodes:
                if inteface_check(node):
                    node.Set(INTERFACE,True)
                    set_reference(node,reference_variable)

            self.GenerateInterface()

            VariableRedistributionUtility.ConvertDistributedValuesToPoint(
                self.interface_model_part,
                reference_variable,
                intermediate_variable)

            VariableRedistributionUtility.DistributePointValues(
                self.interface_model_part,
                intermediate_variable,
                result_variable,
                self.redistribution_tolerance,
                self.redistribution_iterations)

            if self.print_output:
                self.InitializeOutput()
                self.PrintOutput()
                self.FinalizeOutput()

            if KratosGlobals.Kernel.HasDoubleVariable(reference_variable.Name()):
                self.CheckDoubleResults(reference_variable,result_variable)
            elif KratosGlobals.Kernel.HasArrayVariable(reference_variable.Name()):
                self.CheckArrayResults(reference_variable,result_variable)
            else:
                self.fail("Failing due to incorrect test definition: Wrong variable type")

    def testNodalArea(self):
        self.input_file = "square10"

        with WorkFolderScope(self.work_folder):

            self.SetUpProblem()

            for node in self.model_part.Nodes:
                node.Set(INTERFACE,True)
                node.SetSolutionStepValue(PRESSURE,1.0)

            self.GenerateInterface()

            VariableRedistributionUtility.ConvertDistributedValuesToPoint(
                self.interface_model_part,
                PRESSURE,
                TEMPERATURE)

            for cond in self.model_part.Conditions:
                area = cond.GetArea()
                for node in cond.GetNodes():
                    nodal_area = node.GetSolutionStepValue(NODAL_PAUX)
                    node.SetSolutionStepValue(NODAL_PAUX,nodal_area+area/3.0)

            if self.print_output:
                self.InitializeOutput()
                self.PrintOutput()
                self.FinalizeOutput()

            self.CheckDoubleResults(TEMPERATURE,NODAL_PAUX)

    def SetUpProblem(self):
        current_model = Model()
        self.model_part = current_model.CreateModelPart("Model")
        self.model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
        self.model_part.AddNodalSolutionStepVariable(PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(NODAL_PAUX)
        self.model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(MAPPER_SCALAR_PROJECTION_RHS)
        self.model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(VORTICITY)
        self.model_part.AddNodalSolutionStepVariable(ACCELERATION)
        self.model_part.AddNodalSolutionStepVariable(MAPPER_VECTOR_PROJECTION_RHS)

        model_part_io = ModelPartIO(self.input_file)
        model_part_io.ReadModelPart(self.model_part)

        self.model_part.SetBufferSize(1)

    def GenerateInterface(self):
        current_model = Model()
        self.interface_model_part = current_model.CreateModelPart("Interface")

        interface_model_part_generator = InterfacePreprocess()
        interface_model_part_generator.GenerateTriangleInterfacePart(self.model_part,self.interface_model_part)


    def CheckDoubleResults(self,reference_variable,result_variable):

        for node in self.interface_model_part.Nodes:
            reference = node.GetSolutionStepValue(reference_variable)
            result = node.GetSolutionStepValue(result_variable)
            self.assertAlmostEqual(reference, result, delta=self.check_tolerance)

    def CheckArrayResults(self,reference_variable,result_variable):
        for node in self.interface_model_part.Nodes:
            reference = node.GetSolutionStepValue(reference_variable)
            result = node.GetSolutionStepValue(result_variable)
            self.assertAlmostEqual(reference[0], result[0], delta=self.check_tolerance)
            self.assertAlmostEqual(reference[1], result[1], delta=self.check_tolerance)
            self.assertAlmostEqual(reference[2], result[2], delta=self.check_tolerance)

    def InitializeOutput(self):
        gid_mode = GiDPostMode.GiD_PostBinary
        multifile = MultiFileFlag.SingleFile
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteConditions
        self.gid_io = GidIO(self.input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)

        mesh_name = 0.0
        self.gid_io.InitializeMesh( mesh_name)
        self.gid_io.WriteMesh( self.model_part.GetMesh() )
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults(mesh_name,(self.model_part).GetMesh())

    def FinalizeOutput(self):
        self.gid_io.FinalizeResults()

    def PrintOutput(self):
        label = self.model_part.ProcessInfo[TIME]
        self.gid_io.WriteNodalResults(PRESSURE,self.model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(NODAL_PAUX,self.model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(TEMPERATURE,self.model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(MAPPER_SCALAR_PROJECTION_RHS,self.model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(VELOCITY,self.model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(VORTICITY,self.model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(ACCELERATION,self.model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(MAPPER_VECTOR_PROJECTION_RHS,self.model_part.Nodes,label,0)


if __name__ == '__main__':
    test = VariableRedistributionTest()
    test.setUp()
    test.print_output = True
    #test.testLinearFunction()
    #test.testSharpCorners()
    #test.testVector()
    #test.testQuadratic()
    test.testNodalArea()
    test.tearDown()