from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

import KratosMultiphysics.KratosUnittest as UnitTest

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class DarcyChannelTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 3
        self.input_file = "darcy_channel_test"
        self.reference_file = "darcy_channel_reference"
        self.work_folder = "TwoFluidDarcyTest"

        self.xmin = 0.0
        self.xmax = 5.0

        self.dt = 1.0
        self.nsteps = 10
        self.linear_darcy_coefficient = 0.0
        self.nonlinear_darcy_coefficient = 0.0

        self.rho = 1.0
        self.nu = 1.0
        self.u0 = 2.0

        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

    def tearDown(self):
        import os
        with WorkFolderScope(self.work_folder):
            try:
                os.remove(self.input_file+'.time')
            except FileNotFoundError as e:
                pass

    def testDarcyChannel(self):

        with WorkFolderScope(self.work_folder):
            self.setUpModel()
            self.setUpSolver()
            self.setUpProblem()

            if self.print_output:
                self.InitializeOutput()

            self.runTest()

            if self.print_output:
                self.FinalizeOutput()

            self.checkResults()

    def testDarcyLinear(self):
        self.u0 = 2.0
        self.linear_darcy_coefficient = 1.0
        self.nonlinear_darcy_coefficient = 0.0
        self.reference_file = "reference_darcy_linear"
        self.testDarcyChannel()

    def testDarcyNonLinear(self):
        self.u0 = 2.0
        self.linear_darcy_coefficient = 0.0
        self.nonlinear_darcy_coefficient = 1.0
        self.reference_file = "reference_darcy_nonlinear"
        self.testDarcyChannel()

    def testDarcyDensity(self):
        self.u0 = 2.0
        self.linear_darcy_coefficient = 1.0
        self.nonlinear_darcy_coefficient = 1.0
        self.rho = 1000.0
        self.reference_file = "reference_darcy_density"
        self.testDarcyChannel()

    

    def testReferenceValues(self):
        '''
        Table III Dynamic viscosity of water, pure aluminium and A356 alloy.
                                    Water       Pure Aluminium  A356
        Temperature (C)             7           710             710
        Dynamic Viscosity (Pa S)    1.38E-03    1.25E-03        1.03E-03
        Density (kg/m3)             1000        2386            2340
        Kinematic Viscosity (m2/s)  1.38E-06    5.25E-07        4.41E-07
        '''

        class Fluid(object):
            def __init__(self,name,temperature,density,kinematic_viscosity):
                self.name = str(name)
                self.density = float(density)
                self.temperature = float(temperature)
                self.kinematic_viscosity = float(kinematic_viscosity)
 
        fluids = [
            Fluid("Water",7,1000,1.38E-06),
            #Fluid("Pure Aluminium",710,2386,5.25E-07),
            #Fluid("A356",710,2340,4.41E-07),
        ]

        '''
        Table XII Average value of k1 and k2 for different types of filters
        Filter Type Forchheimer k1 (m)  Forchheimer k2 (m2)
        30 ppi      4.339E-08           5.086E-04
        40 ppi      3.099E-08           3.379E-04
        50 ppi      1.748E-08           1.960E-04
        80 ppi      6.352E-09           1.094E-04
        '''       

        class Filter(object):
            def __init__(self,name,k1,k2):
                self.name = str(name)
                self.k1 = float(k1)
                self.k2 = float(k2)

        filters = [
            Filter("30 ppi",4.339E-08,5.086E-04),
            #Filter("40 ppi",3.099E-08,3.379E-04),
            #Filter("50 ppi",1.748E-08,1.960E-04),
            #Filter("80 ppi",6.352E-09,1.094E-04)
        ]

        for fluid in fluids:
            for filt in filters:
                self.runTestCase(fluid,filt)

    def runTestCase(self,fluid,filt):
        self.rho = fluid.density
        self.nu = fluid.kinematic_viscosity
        self.dt = 0.01
        self.nsteps = 100

        self.linear_darcy_coefficient = self.rho * self.nu / filt.k1
        self.nonlinear_darcy_coefficient = self.rho / filt.k2

        # dP/dX = (mu/k1) * u0 + (rho/k2) * u0**2
        self.testDarcyChannel()

        for node in self.fluid_model_part.Nodes:
            if node.Id == 13:
                p_in = node.GetSolutionStepValue(PRESSURE)
            elif node.Id == 435:
                p_out = node.GetSolutionStepValue(PRESSURE)

        expected_pressure_drop = (self.xmax - self.xmin) * (self.linear_darcy_coefficient*self.u0 + self.nonlinear_darcy_coefficient*self.u0**2)
        measured_pressure_drop = p_in - p_out
        print("Expected drop: {0} Measured drop: {1}".format(expected_pressure_drop,measured_pressure_drop))



    def setUpModel(self):

        self.fluid_model_part = ModelPart("Fluid")

        self.fluid_model_part.Properties[0].SetValue(LIN_DARCY_COEF,self.linear_darcy_coefficient)
        self.fluid_model_part.Properties[0].SetValue(NONLIN_DARCY_COEF,self.nonlinear_darcy_coefficient)

    def setUpSolver(self):
        oss_switch = 0

        import vms_monolithic_solver
        vms_monolithic_solver.AddVariables(self.fluid_model_part)
        self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

        model_part_io = ModelPartIO(self.input_file)
        model_part_io.ReadModelPart(self.fluid_model_part)

        self.fluid_model_part.SetBufferSize(3)
        vms_monolithic_solver.AddDofs(self.fluid_model_part)

        # Building custom fluid solver
        self.fluid_solver = vms_monolithic_solver.MonolithicSolver(self.fluid_model_part,self.domain_size)
        rel_vel_tol = 1e-5
        abs_vel_tol = 1e-7
        rel_pres_tol = 1e-5
        abs_pres_tol = 1e-7
        self.fluid_solver.conv_criteria = VelPrCriteria(rel_vel_tol,abs_vel_tol,rel_pres_tol,abs_pres_tol)

        alpha = -0.3
        move_mesh = 0
        self.fluid_solver.time_scheme = ResidualBasedPredictorCorrectorBDFSchemeTurbulentNoReaction(self.domain_size)
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

        self.fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH,oss_switch)

        self.fluid_solver.divergence_clearance_steps = 0
        self.fluid_solver.use_slip_conditions = 0


    def setUpProblem(self):

        ymin = 0.0
        ymax = 1.0
        zmin = 0.0
        zmax = 1.0

        ## Set initial and boundary conditions
        for node in self.fluid_model_part.Nodes:
            node.SetSolutionStepValue(DENSITY,self.rho)
            node.SetSolutionStepValue(VISCOSITY,self.nu)
            #node.SetSolutionStepValue(DISTANCE,1000.0)
            #node.SetSolutionStepValue(BODY_FORCE_X,self.rho*self.u0)

            if node.X == self.xmin:
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)
                node.SetSolutionStepValue(VELOCITY_X,0,self.u0)
            elif node.X == self.xmax:
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)


            else:
                if node.Y == ymin or node.Y == ymax:
                    node.Fix(VELOCITY_Y)
                if node.Z == zmin or node.Z == zmax:
                    node.Fix(VELOCITY_Z)

    def runTest(self):
        time = 0.0

        for step in range(self.nsteps):
            time = time+self.dt
            self.fluid_model_part.CloneTimeStep(time)
            if step > 0:
                self.fluid_solver.Solve()

            if self.print_output:
                self.printOutput()

    def checkResults(self):

        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, PRESSURE\n")
                for node in self.fluid_model_part.Nodes:
                    value = node.GetSolutionStepValue(PRESSURE,0)
                    ref_file.write("{0}, {1}\n".format(node.Id, value))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()
                node_iter = self.fluid_model_part.Nodes

                for node in self.fluid_model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    node_id = values[0]
                    reference_pressure = values[1]

                    solution_pressure = node.GetSolutionStepValue(PRESSURE)
                    self.assertAlmostEqual(reference_pressure, solution_pressure, delta=self.check_tolerance)

                    line = reference_file.readline()
                if line != '': # If we did not reach the end of the reference file
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

    def InitializeOutput(self):
        gid_mode = GiDPostMode.GiD_PostBinary
        multifile = MultiFileFlag.SingleFile
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteElementsOnly
        self.gid_io = GidIO(self.input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)

        mesh_name = 0.0
        self.gid_io.InitializeMesh( mesh_name)
        self.gid_io.WriteMesh( self.fluid_model_part.GetMesh() )
        self.gid_io.FinalizeMesh()

    def printOutput(self):
        label = self.fluid_model_part.ProcessInfo[TIME]
        self.gid_io.WriteNodalResults(VELOCITY,self.fluid_model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(PRESSURE,self.fluid_model_part.Nodes,label,0)
    
    def FinalizeOutput(self):
        self.gid_io.FinalizeResults()

if __name__ == '__main__':
    test = DarcyChannelTest()
    test.setUp()
    test.print_reference_values = True
    test.print_output = True
    #test.testDarcyLinear()
    #test.testDarcyNonLinear()
    #test.testDarcyDensity()
    test.testReferenceValues()
    test.tearDown()
