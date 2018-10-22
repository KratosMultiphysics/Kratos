from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

import KratosMultiphysics.KratosUnittest as UnitTest
        
import vms_monolithic_solver

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

class FluidElementTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.input_file = "cavity10"
        self.reference_file = "reference10_qasgs"
        self.work_folder = "FluidElementTest"
        self.element = "QSVMS2D3N"

        self.dt = 0.1
        self.nsteps = 10

        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

        self.oss_switch = 0

    def tearDown(self):
        import os
        with WorkFolderScope(self.work_folder):
            try:
                os.remove(self.input_file+'.time')
            except FileNotFoundError as e:
                pass

    def testCavity(self):

        with WorkFolderScope(self.work_folder):
            self.setUpModel()
            self.setUpProblem()
            self.setUpSolvers()

            self.runTest()

            self.checkResults()

    
    def testCavityQSASGS(self):
        self.reference_file = "reference10_qsasgs"
        self.element = "QSVMS2D3N"
        self.oss_switch = 0
        self.testCavity()

    def testCavityQSOSS(self):
        self.reference_file = "reference10_qsoss"
        self.element = "QSVMS2D3N"
        self.oss_switch = 1
        self.testCavity()
    
    def testCavityDASGS(self):
        self.reference_file = "reference10_dasgs"
        self.element = "DVMS2D3N"
        self.oss_switch = 0
        self.testCavity()

    def testCavityDOSS(self):
        self.reference_file = "reference10_doss"
        self.element = "DVMS2D3N"
        self.oss_switch = 1
        self.testCavity()

    def setUpModel(self):
        self.model = Model()
        self.fluid_model_part = self.model.CreateModelPart("Fluid")

        vms_monolithic_solver.AddVariables(self.fluid_model_part)

        model_part_io = ModelPartIO(self.input_file)
        model_part_io.ReadModelPart(self.fluid_model_part)

        self.fluid_model_part.SetBufferSize(2)
        vms_monolithic_solver.AddDofs(self.fluid_model_part)

        ## Replace element and conditions
        replace_settings = Parameters("""{
            "condition_name": "MonolithicWallCondition2D"
        }""")
        replace_settings.AddEmptyValue("element_name")
        replace_settings["element_name"].SetString(self.element)

        ReplaceElementsAndConditionsProcess(self.fluid_model_part, replace_settings).Execute()

    def setUpSolvers(self):

        # Building custom fluid solver
        self.fluid_solver = vms_monolithic_solver.MonolithicSolver(self.fluid_model_part,self.domain_size)
        rel_vel_tol = 1e-5
        abs_vel_tol = 1e-7
        rel_pres_tol = 1e-5
        abs_pres_tol = 1e-7
        self.fluid_solver.conv_criteria = VelPrCriteria(rel_vel_tol,abs_vel_tol,rel_pres_tol,abs_pres_tol)
        self.fluid_solver.conv_criteria.SetEchoLevel(0)

        alpha = -0.3
        move_mesh = 0
        self.fluid_solver.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(alpha,move_mesh,self.domain_size)
        precond = DiagonalPreconditioner()
        self.fluid_solver.linear_solver = BICGSTABSolver(1e-6, 5000, precond)
        builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.fluid_solver.linear_solver)
        self.fluid_solver.max_iter = 50
        self.fluid_solver.compute_reactions = False
        self.fluid_solver.ReformDofSetAtEachStep = False
        self.fluid_solver.MoveMeshFlag = False

        self.fluid_solver.solver = ResidualBasedNewtonRaphsonStrategy(\
                self.fluid_model_part,
                self.fluid_solver.time_scheme,
                self.fluid_solver.linear_solver,
                self.fluid_solver.conv_criteria,
                builder_and_solver,
                self.fluid_solver.max_iter,
                self.fluid_solver.compute_reactions,
                self.fluid_solver.ReformDofSetAtEachStep,
                self.fluid_solver.MoveMeshFlag)

        self.fluid_solver.solver.SetEchoLevel(0)
        self.fluid_solver.solver.Check()

        self.fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH,self.oss_switch)

        self.fluid_solver.divergence_clearance_steps = 0
        self.fluid_solver.use_slip_conditions = 0


    def setUpProblem(self):
        xmin = 0.0
        xmax = 1.0
        ymin = 0.0
        ymax = 1.0

        ux = 1.0

        for element in self.fluid_model_part.Elements:
            element.Initialize()

        ## Set initial and boundary conditions
        for node in self.fluid_model_part.Nodes:

            if node.X == xmin or node.X == xmax or node.Y == ymin or node.Y == ymax:
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                if node.X == xmin and node.Y == ymin:
                    node.Fix(PRESSURE)

            if node.Y == ymax:
                node.SetSolutionStepValue(VELOCITY_X,ux)

    def runTest(self):
        time = 0.0

        if self.print_output:
            self.__InitializeOutput()

        for step in range(self.nsteps):
            time = time+self.dt
            self.fluid_model_part.CloneTimeStep(time)
            self.fluid_solver.Solve()

            if self.print_output:
                self.__PrintResults()

        if self.print_output:
            self.__FinializeOutput()


    def checkResults(self):

        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, VELOCITY_X, VELOCITY_Y, PRESSURE\n")
                for node in self.fluid_model_part.Nodes:
                    vel = node.GetSolutionStepValue(VELOCITY,0)
                    temp = node.GetSolutionStepValue(PRESSURE,0)
                    ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], temp))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()

                for node in self.fluid_model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    #node_id = values[0]
                    reference_vel_x = values[1]
                    reference_vel_y = values[2]
                    reference_press = values[3]

                    velocity = node.GetSolutionStepValue(VELOCITY)
                    self.assertAlmostEqual(reference_vel_x, velocity[0], delta=self.check_tolerance)
                    self.assertAlmostEqual(reference_vel_y, velocity[1], delta=self.check_tolerance)
                    pressure = node.GetSolutionStepValue(PRESSURE)
                    self.assertAlmostEqual(reference_press, pressure, delta=self.check_tolerance)

                    line = reference_file.readline()
                if line != '': # If we did not reach the end of the reference file
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

    def __InitializeOutput(self):
        gid_mode = GiDPostMode.GiD_PostBinary
        multifile = MultiFileFlag.SingleFile
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteElementsOnly
        self.gid_io = GidIO(self.input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)

        mesh_name = 0.0
        self.gid_io.InitializeMesh( mesh_name)
        self.gid_io.WriteMesh( self.fluid_model_part.GetMesh() )
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults(mesh_name,(self.fluid_model_part).GetMesh())

    def __PrintResults(self):
        label = self.fluid_model_part.ProcessInfo[TIME]
        self.gid_io.WriteNodalResults(VELOCITY,self.fluid_model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(PRESSURE,self.fluid_model_part.Nodes,label,0)
        self.gid_io.PrintOnGaussPoints(SUBSCALE_VELOCITY,self.fluid_model_part,label)
        self.gid_io.PrintOnGaussPoints(SUBSCALE_PRESSURE,self.fluid_model_part,label)

    def __FinializeOutput(self):
        self.gid_io.FinalizeResults()

if __name__ == '__main__':
    test = FluidElementTest()
    test.setUp()
    test.print_reference_values = False
    test.print_output = True
    #test.testCavityQSASGS()
    #test.testCavityQSOSS()
    #test.testCavityDASGS()
    test.testCavityDOSS()
    test.tearDown()
