from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities

have_convection_diffusion = KratosUtilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication")
if have_convection_diffusion:
    import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff

import KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic as navier_stokes_solver

@UnitTest.skipUnless(have_convection_diffusion,"Missing required application: ConvectionDiffusionApplication")
class BuoyancyTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.reference_file = "reference10"
        self.input_file = "cavity10"

        self.convection_diffusion_solver = "eulerian"
        self.nsteps = 3
        self.thermal_expansion_coefficient = None # If set, it will be used instead of 1./AmbientTemperature

        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

    def tearDown(self):
        with UnitTest.WorkFolderScope("BuoyancyTest", __file__):
            KratosUtilities.DeleteFileIfExisting(self.input_file+'.time')

    def testEulerian(self):
        self.convection_diffusion_solver = "eulerian"
        self.reference_file = "reference10_eulerian"
        self.testBuoyancy()


    def testThermalExpansionCoefficient(self):
        self.convection_diffusion_solver = "eulerian"
        self.reference_file = "reference10_eulerian"
        self.thermal_expansion_coefficient = 1./293.15
        self.testBuoyancy()

    def testBFECC(self):
        self.convection_diffusion_solver = "bfecc"
        self.reference_file = "reference10_bfecc"
        self.check_tolerance = 1e-1 # The bfecc solver shows some variation between runs, we cannot be too strict here
        self.testBuoyancy()

    def validationEulerian(self):
        self.input_file = "cavity80"
        self.reference_file = "reference80_eulerian"

        self.convection_diffusion_solver = "eulerian"
        self.nsteps = 200

        self.testBuoyancy()

    def testBuoyancy(self):
        with UnitTest.WorkFolderScope("BuoyancyTest", __file__):
            self.setUpModel()
            self.setUpSolvers()
            self.setUpProblem()

            self.runTest()

            self.checkResults()
            if self.print_output:
                self.printOutput()

    def setUpModel(self):
        self.model = Model()

        with open("solver_settings.json","r") as settings_file:
            settings = Parameters(settings_file.read())

        settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(self.input_file)
        thermal_settings = ConvectionDiffusionSettings()
        thermal_settings.SetUnknownVariable(TEMPERATURE)
        thermal_settings.SetDensityVariable(DENSITY)
        thermal_settings.SetSpecificHeatVariable(SPECIFIC_HEAT)
        thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
        #thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
        #thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
        thermal_settings.SetVelocityVariable(VELOCITY)
        thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
        if self.convection_diffusion_solver == 'bfecc':
            thermal_settings.SetProjectionVariable(ConvDiff.PROJECTED_SCALAR1)

        self.fluid_solver = navier_stokes_solver.CreateSolver(self.model,settings["solver_settings"])

        print(self.model)
        self.fluid_model_part = self.model.GetModelPart(settings["solver_settings"]["model_part_name"].GetString())
        self.fluid_model_part.ProcessInfo.SetValue(CONVECTION_DIFFUSION_SETTINGS,thermal_settings)

    def setUpSolvers(self):

        if self.convection_diffusion_solver == 'bfecc':
            import KratosMultiphysics.ConvectionDiffusionApplication.bfecc_convection_diffusion_solver as thermal_solver
        elif self.convection_diffusion_solver == 'eulerian':
            import KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_solver as thermal_solver
        else:
            raise Exception("Unsupported convection-diffusion solver option: {0}".format(self.convection_diffusion_solver))

        self.fluid_solver.AddVariables()
        thermal_solver.AddVariables(self.fluid_model_part)

        self.fluid_solver.ImportModelPart()
        self.fluid_solver.PrepareModelPart()

        self.fluid_model_part.SetBufferSize(2)

        self.fluid_solver.AddDofs()
        thermal_solver.AddDofs(self.fluid_model_part)

        self.fluid_solver.Initialize()

        if self.convection_diffusion_solver == 'eulerian':
            # Duplicate model part

            thermal_model_part = self.model.CreateModelPart("Thermal")
            conv_diff_element = "EulerianConvDiff2D"
            conv_diff_condition = "LineCondition2D2N"

            MergeVariableListsUtility().Merge(self.fluid_model_part, thermal_model_part)

            modeler = ConnectivityPreserveModeler()
            modeler.GenerateModelPart(self.fluid_model_part,thermal_model_part,conv_diff_element,conv_diff_condition)

            # thermal solver
            self.thermal_solver = thermal_solver.ConvectionDiffusionSolver(thermal_model_part,self.domain_size)
        else:
            class SolverSettings:
                def __init__(self,domain_size):
                    self.domain_size = domain_size
            settings = SolverSettings(self.domain_size)
            self.thermal_solver = thermal_solver.CreateSolver(self.model_part,settings)
        self.thermal_solver.Initialize()


    def setUpProblem(self):
        xmin = 0.0
        xmax = 1.0
        ymin = 0.0
        ymax = 1.0

        ## For Ra~1e6
        g = 9.81     # accelertion of gravity m/s2
        T1 = 293.15  # Cold (reference) temperature K
        T2 = 303.15  # Hot temperature K
        rho = 1.2039 # (reference) density kg/m3
        c = 1004.84  # Specific heat J/kg K
        k = ( (rho*c)**2*g*(T2-T1)*(xmax-xmin)**3 / (1e6*T1*0.71) )**0.5 # Given Ra=1e6 & Pr=0.71
        mu = 0.71*k/c # For Prandlt = 0.71
        nu = mu/rho

        if self.thermal_expansion_coefficient is not None:
            parameter_string = '{ "gravity" : [ 0.0, -' + str(g) +', 0.0 ], "thermal_expansion_coefficient" : '+ str(self.thermal_expansion_coefficient) +' }'
        else:
            parameter_string = '{ "gravity" : [ 0.0, -' + str(g) +', 0.0 ] }'
        parameters = Parameters(parameter_string)
        self.buoyancy_process = BoussinesqForceProcess(self.fluid_model_part,parameters)

        ## Set initial and boundary conditions
        self.fluid_model_part.ProcessInfo.SetValue(AMBIENT_TEMPERATURE,T1)
        for node in self.fluid_model_part.Nodes:
            node.SetSolutionStepValue(DENSITY,rho)
            node.SetSolutionStepValue(VISCOSITY,nu)
            node.SetSolutionStepValue(CONDUCTIVITY,k)
            node.SetSolutionStepValue(SPECIFIC_HEAT,c)

            if node.X == xmin or node.X == xmax or node.Y == ymin or node.Y == ymax:
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                if node.X == xmin and node.Y == ymin:
                    node.Fix(PRESSURE)

            if node.X == xmin:
                node.Fix(TEMPERATURE)
                node.SetSolutionStepValue(TEMPERATURE,T2)
            elif node.X == xmax:
                node.Fix(TEMPERATURE)
                node.SetSolutionStepValue(TEMPERATURE,T1)
            else:
                T = T2 + (T1-T2)*(node.X-xmin)/(xmax-xmin)
                node.SetSolutionStepValue(TEMPERATURE,T)
                #node.SetSolutionStepValue(TEMPERATURE,T1)

        self.buoyancy_process.ExecuteInitialize()

    def runTest(self):
        time = 0.0

        for step in range(self.nsteps):
            time = self.fluid_solver.AdvanceInTime(time)
            self.buoyancy_process.ExecuteInitializeSolutionStep()

            self.fluid_solver.InitializeSolutionStep()
            self.fluid_solver.Predict()
            self.fluid_solver.SolveSolutionStep()
            self.fluid_solver.FinalizeSolutionStep()

            self.thermal_solver.Solve()

    def checkResults(self):

        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, VELOCITY_X, VELOCITY_Y, TEMPERATURE\n")
                for node in self.fluid_model_part.Nodes:
                    vel = node.GetSolutionStepValue(VELOCITY,0)
                    temp = node.GetSolutionStepValue(TEMPERATURE,0)
                    ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], temp))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()
                node_iter = self.fluid_model_part.Nodes

                for node in self.fluid_model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    node_id = values[0]
                    reference_vel_x = values[1]
                    reference_vel_y = values[2]
                    reference_temp = values[3]

                    velocity = node.GetSolutionStepValue(VELOCITY)
                    self.assertAlmostEqual(reference_vel_x, velocity[0], delta=self.check_tolerance)
                    self.assertAlmostEqual(reference_vel_y, velocity[1], delta=self.check_tolerance)
                    temperature = node.GetSolutionStepValue(TEMPERATURE)
                    self.assertAlmostEqual(reference_temp, temperature, delta=self.check_tolerance)

                    line = reference_file.readline()
                if line != '': # If we did not reach the end of the reference file
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

    def printOutput(self):
        gid_mode = GiDPostMode.GiD_PostBinary
        multifile = MultiFileFlag.SingleFile
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteElementsOnly
        gid_io = GidIO(self.input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)

        mesh_name = 0.0
        gid_io.InitializeMesh( mesh_name)
        gid_io.WriteMesh( self.fluid_model_part.GetMesh() )
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name,(self.fluid_model_part).GetMesh())

        label = self.fluid_model_part.ProcessInfo[TIME]
        gid_io.WriteNodalResults(VELOCITY,self.fluid_model_part.Nodes,label,0)
        gid_io.WriteNodalResults(PRESSURE,self.fluid_model_part.Nodes,label,0)
        gid_io.WriteNodalResults(TEMPERATURE,self.fluid_model_part.Nodes,label,0)
        gid_io.WriteNodalResults(DENSITY,self.fluid_model_part.Nodes,label,0)
        gid_io.WriteNodalResults(VISCOSITY,self.fluid_model_part.Nodes,label,0)
        gid_io.WriteNodalResults(CONDUCTIVITY,self.fluid_model_part.Nodes,label,0)
        gid_io.WriteNodalResults(SPECIFIC_HEAT,self.fluid_model_part.Nodes,label,0)
        gid_io.WriteNodalResults(BODY_FORCE,self.fluid_model_part.Nodes,label,0)

        gid_io.FinalizeResults()

if __name__ == '__main__':
    test = BuoyancyTest()
    test.setUp()
    #test.print_reference_values = True
    test.print_output = True
    #test.testEulerian()
    test.validationEulerian()
    test.tearDown()
