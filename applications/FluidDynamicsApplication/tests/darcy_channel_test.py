from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities

import KratosMultiphysics.FluidDynamicsApplication.navier_stokes_two_fluids_solver as two_fluids_solver

class TwoFluidNoRedistanceSolver(two_fluids_solver.NavierStokesTwoFluidsSolver):
    """
    Ad-hoc solver skipping distance levelset re-calculation.
    Since this test does not involve a free surface, it is more practical to just skip its convection.
    """
    def __init__(self, model, settings):
        super(TwoFluidNoRedistanceSolver,self).__init__(model,settings)

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # Recompute the BDF2 coefficients
            (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

            # Initialize the solver current step
            self._GetSolutionStrategy().InitializeSolutionStep()


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
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting(self.input_file+'.time')

    def runDarcyChannelTest(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
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
        self.runDarcyChannelTest()

    def testDarcyNonLinear(self):
        self.u0 = 2.0
        self.linear_darcy_coefficient = 0.0
        self.nonlinear_darcy_coefficient = 1.0
        self.runDarcyChannelTest()

    def testDarcyDensity(self):
        self.u0 = 2.0
        self.rho = 1000.0
        self.linear_darcy_coefficient = 1.0/self.rho
        self.nonlinear_darcy_coefficient = 1.0/self.rho
        self.runDarcyChannelTest()

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

        self.linear_darcy_coefficient = 1.0 / filt.k1
        self.nonlinear_darcy_coefficient = 1.0 / filt.k2

        Logger.PrintInfo("Darcy Test","A: {0} B: {1}".format(self.linear_darcy_coefficient,self.nonlinear_darcy_coefficient))

        self.do_check = False # override default verification function
        self.runDarcyChannelTest()

        # for the mesh in 'TwoFluidDarcyTest', node 13 is in the inlet, node 453 in the outlet
        for node in self.fluid_model_part.Nodes:
            if node.GetSolutionStepValue(FLAG_VARIABLE) == 1.0:
                p_in = node.GetSolutionStepValue(PRESSURE)
            elif node.GetSolutionStepValue(FLAG_VARIABLE) == 2.0:
                p_out = node.GetSolutionStepValue(PRESSURE)


        # dP/dX = (mu/k1) * u0 + (rho/k2) * u0**2
        expected_pressure_drop = (self.xmax - self.xmin) * (self.nu*self.rho*self.linear_darcy_coefficient*self.u0 + self.rho*self.nonlinear_darcy_coefficient*self.u0**2)
        measured_pressure_drop = p_in - p_out
        rel_error = 100. * (measured_pressure_drop-expected_pressure_drop)/expected_pressure_drop
        outfile.write("{0}; {1}; {2}; {3}; {4}; {5}; {6}\n".format(fluid.name,filt.name,self.dt,self.dynamic_tau,expected_pressure_drop,measured_pressure_drop,rel_error))

    def setUpModel(self):
        self.model = Model()

        with open("solver_settings.json","r") as settings_file:
            settings = Parameters(settings_file.read())

        settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(self.input_file)
        settings["solver_settings"]["formulation"]["dynamic_tau"].SetDouble(self.dynamic_tau)
        settings["solver_settings"]["time_stepping"]["time_step"].SetDouble(self.dt)

        self.fluid_solver = TwoFluidNoRedistanceSolver(self.model, settings["solver_settings"])
        self.fluid_model_part = self.model.GetModelPart("Fluid")

    def setUpSolver(self):

        self.fluid_solver.AddVariables()
        self.fluid_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)

        self.fluid_solver.ImportModelPart()

        mu = self.rho * self.nu
        for prop in self.fluid_model_part.Properties:
            prop.SetValue(DENSITY, self.rho)
            prop.SetValue(DYNAMIC_VISCOSITY, mu)
            prop.SetValue(LIN_DARCY_COEF, self.linear_darcy_coefficient)
            prop.SetValue(NONLIN_DARCY_COEF, self.nonlinear_darcy_coefficient)
            prop.SetValue(CONSTITUTIVE_LAW,NewtonianTwoFluid3DLaw())

        self.fluid_solver.PrepareModelPart()
        self.fluid_solver.AddDofs()

        self.fluid_solver.Initialize()

    def setUpProblem(self):
        ## Set initial and boundary conditions
        for node in self.fluid_model_part.Nodes:
            node.SetSolutionStepValue(DISTANCE, 0, -100.0)

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
            time = self.fluid_solver.AdvanceInTime(time)

            self.fluid_solver.InitializeSolutionStep()
            self.fluid_solver.Predict()
            self.fluid_solver.SolveSolutionStep()
            self.fluid_solver.FinalizeSolutionStep()

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
            expected_pressure_drop = (self.xmax-self.xmin) * (self.nu*self.rho*self.linear_darcy_coefficient*self.u0 + self.rho*self.nonlinear_darcy_coefficient*self.u0**2)
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
        self.gid_io.WriteNodalResults(DISTANCE,self.fluid_model_part.Nodes,label,0)

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
