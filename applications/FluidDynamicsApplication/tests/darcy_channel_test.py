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
        self.ymin = 0.0
        self.ymax = 1.0
        self.zmin = 0.0
        self.zmax = 1.0

        self.dt = 1.0
        self.nsteps = 10
        self.linear_darcy_coefficient = 0.0
        self.nonlinear_darcy_coefficient = 0.0
        self.dynamic_tau = 1.0

        self.rho = 1.0
        self.nu = 1.0
        self.u0 = 2.0

        self.do_check = True
        self.check_tolerance = 1e-6
        self.print_output = False

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

            if self.do_check:
                self.checkResults()

    def testDarcyLinear(self):
        self.u0 = 2.0
        self.linear_darcy_coefficient = 1.0
        self.nonlinear_darcy_coefficient = 0.0
        self.testDarcyChannel()

    def testDarcyNonLinear(self):
        self.u0 = 2.0
        self.linear_darcy_coefficient = 0.0
        self.nonlinear_darcy_coefficient = 1.0
        self.testDarcyChannel()

    def testDarcyDensity(self):
        self.u0 = 2.0
        self.linear_darcy_coefficient = 1.0
        self.nonlinear_darcy_coefficient = 1.0
        self.rho = 1000.0
        self.testDarcyChannel()

    def testReferenceValues(self):
        '''
        Reference values from K. Zhang (2012) "Liquid Permeability of Ceramic Foam Filters"
        Master thesis at the Norwegian University of Science and Technology
        Retrieved from
        https://daim.idi.ntnu.no/masteroppgaver/008/8358/masteroppgave.pdf
        '''

        '''
        Table III Dynamic viscosity of water, pure aluminium and A356 alloy.
                                    Water       Pure Aluminium  A356
        Temperature (C)             7           710             710
        Dynamic Viscosity (Pa s)    1.38E-03    1.25E-03        1.03E-03
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
            Fluid("Pure Aluminium",710,2386,5.25E-07),
            Fluid("A356",710,2340,4.41E-07),
        ]

        '''
        These are ceramic filters:
        Table XII Average value of k1 and k2 for different types of filters
        Filter Type Forchheimer k1 (m2) Forchheimer k2 (m)
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
            Filter("Ceramic 30 ppi",4.339E-08,5.086E-04),
            Filter("Ceramic 40 ppi",3.099E-08,3.379E-04),
            Filter("Ceramic 50 ppi",1.748E-08,1.960E-04),
            Filter("Ceramic 80 ppi",6.352E-09,1.094E-04),
            Filter("Cellular 100 csi",3.52E-08,1.41E-02),
            Filter("Cellular 200 csi",1.87E-08,3.00E-02),
            Filter("Cellular 300 csi",1.71E-08,8.84E-03),
            Filter("Foam 10 ppi",2.14E-08,2.66E-03),
            Filter("Foam 20 ppi",2.62E-08,1.69E-03),
            Filter("Foam 30 ppi",1.79E-08,1.48E-03),
        ]

        self.dt = 1e-2 #1.0
        self.nsteps = 100
        self.u0 = 0.5

        self.xmin = 0.0
        self.xmax = 0.5
        self.ymin = 0.0
        self.ymax = 0.05
        self.zmin = 0.0
        self.zmax = 0.05

        self.dynamic_tau = 0.0

        self.work_folder = "TwoFluidDarcyValidation"
        self.input_file = "darcy_pressure_drop" #_fine"
        self.reference_file = "pressure_drop.csv"

        with open(self.reference_file,'a') as outfile:
            outfile.write("#Fluid; Filter; dt (s); Dynamic Tau; Expected Pressure drop (Pa); Measured Pressure drop (Pa); Relative error (%)\n")
            for fluid in fluids:
                for filt in filters:
                    self.runTestCase(fluid,filt,outfile)
            outfile.write("\n")

    def runTestCase(self,fluid,filt,outfile):
        self.rho = fluid.density
        self.nu = fluid.kinematic_viscosity

        self.linear_darcy_coefficient = self.rho * self.nu / filt.k1
        self.nonlinear_darcy_coefficient = self.rho / filt.k2

        print("A: {0} B: {1}".format(self.linear_darcy_coefficient,self.nonlinear_darcy_coefficient))

        self.do_check = False # override default verification function
        self.testDarcyChannel()

        # for the mesh in 'TwoFluidDarcyTest', node 13 is in the inlet, node 453 in the outlet
        for node in self.fluid_model_part.Nodes:
            if node.GetSolutionStepValue(FLAG_VARIABLE) == 1.0:
                p_in = node.GetSolutionStepValue(PRESSURE)
            elif node.GetSolutionStepValue(FLAG_VARIABLE) == 2.0:
                p_out = node.GetSolutionStepValue(PRESSURE)


        # dP/dX = (mu/k1) * u0 + (rho/k2) * u0**2
        expected_pressure_drop = (self.xmax - self.xmin) * (self.linear_darcy_coefficient*self.u0 + self.nonlinear_darcy_coefficient*self.u0**2)
        measured_pressure_drop = p_in - p_out
        rel_error = 100. * (measured_pressure_drop-expected_pressure_drop)/expected_pressure_drop
        outfile.write("{0}; {1}; {2}; {3}; {4}; {5}; {6}\n".format(fluid.name,filt.name,self.dt,self.dynamic_tau,expected_pressure_drop,measured_pressure_drop,rel_error))

    def setUpModel(self):
        self.model = Model()
        self.fluid_model_part = self.model.CreateModelPart("Fluid")

        self.fluid_model_part.Properties[0].SetValue(LIN_DARCY_COEF,self.linear_darcy_coefficient)
        self.fluid_model_part.Properties[0].SetValue(NONLIN_DARCY_COEF,self.nonlinear_darcy_coefficient)

    def setUpSolver(self):
        oss_switch = 0

        import vms_monolithic_solver
        vms_monolithic_solver.AddVariables(self.fluid_model_part)
        self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
        self.fluid_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)

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
        self.fluid_solver.conv_criteria.SetEchoLevel(0)

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
        self.fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU,self.dynamic_tau)

        self.fluid_solver.divergence_clearance_steps = 0
        self.fluid_solver.use_slip_conditions = 0


    def setUpProblem(self):
        ## Set initial and boundary conditions
        for node in self.fluid_model_part.Nodes:
            node.SetSolutionStepValue(DENSITY,self.rho)
            node.SetSolutionStepValue(VISCOSITY,self.nu)

            if node.X == self.xmin:
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)
                node.SetSolutionStepValue(VELOCITY_X,0,self.u0)
            elif node.X == self.xmax:
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)
            else:
                if node.Y == self.ymin or node.Y == self.ymax:
                    node.Fix(VELOCITY_Y)
                if node.Z == self.zmin or node.Z == self.zmax:
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

        node_in = None
        node_out = None
        for node in self.fluid_model_part.Nodes:
            if node.X == self.xmin:
                node_in = node
            elif node.X == self.xmax:
                node_out = node

        if node_in is not None and node_out is not None:
            p_in = node_in.GetSolutionStepValue(PRESSURE)
            p_out = node_out.GetSolutionStepValue(PRESSURE)

            # dP/dX = (mu/k1) * u0 + (rho/k2) * u0**2
            expected_pressure_drop = (self.xmax-self.xmin) * (self.linear_darcy_coefficient*self.u0 + self.nonlinear_darcy_coefficient*self.u0**2)
            measured_pressure_drop = p_in - p_out
            self.assertAlmostEqual(expected_pressure_drop,measured_pressure_drop,6)
        else:
            self.fail("Could not find inlet and/or outlet nodes: something is wrong with the input mesh.")

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
    test.print_output = True
    #test.testDarcyLinear()
    #test.testDarcyNonLinear()
    test.testDarcyDensity()
    #test.testReferenceValues()
    test.tearDown()
