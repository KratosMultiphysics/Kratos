from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities

import KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic as navier_stokes_solver

class FluidElementTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.input_file = "cavity10"
        self.reference_file = "reference10_qasgs"
        self.work_folder = "FluidElementTest"
        self.element = "qsvms"
        self.is_time_integrated = False

        self.nsteps = 10

        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

        self.oss_switch = 0

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting(self.input_file+'.time')

    def runCavity(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            if self.is_time_integrated:
                self.setUpModelTimeIntegrated()
            else:
                self.setUpModel()
            self.setUpProblem()

            self.runTest()

            self.checkResults()


    def testCavityQSASGS(self):
        self.reference_file = "reference10_qsasgs"
        self.element = "qsvms"
        self.oss_switch = False
        self.runCavity()

    def testCavityQSOSS(self):
        self.reference_file = "reference10_qsoss"
        self.element = "qsvms"
        self.oss_switch = True
        self.runCavity()

    def testCavityDASGS(self):
        self.reference_file = "reference10_dasgs"
        self.element = "dvms"
        self.oss_switch = False
        self.runCavity()

    def testCavityDOSS(self):
        self.reference_file = "reference10_doss"
        self.element = "dvms"
        self.oss_switch = True
        self.runCavity()

    def testSymbolic(self):
        self.element = "symbolic"
        self.reference_file = "reference10_symbolic"
        self.is_time_integrated = True
        self.runCavity()

    def testTimeIntegratedQSVMS(self):
        self.element = "qsvms"
        self.reference_file = "reference10_time_integrated"
        self.is_time_integrated = True
        self.runCavity()

    def setUpModel(self):
        self.model = Model()

        with open("solver_settings.json","r") as settings_file:
            settings = Parameters(settings_file.read())

        settings["solver_settings"]["formulation"]["element_type"].SetString(self.element)
        settings["solver_settings"]["formulation"].AddEmptyValue("use_orthogonal_subscales")
        settings["solver_settings"]["formulation"]["use_orthogonal_subscales"].SetBool(self.oss_switch)

        self.fluid_solver = navier_stokes_solver.CreateSolver(self.model,settings["solver_settings"])
        self.fluid_solver.AddVariables()
        self.fluid_solver.ImportModelPart()
        self.fluid_solver.PrepareModelPart()
        self.fluid_solver.AddDofs()

        fluid_model_part_name = settings["solver_settings"]["model_part_name"].GetString() + "." + settings["solver_settings"]["volume_model_part_name"].GetString()
        self.fluid_model_part = self.model.GetModelPart(fluid_model_part_name)

    def setUpModelTimeIntegrated(self):
        self.model = Model()

        with open("solver_settings.json","r") as settings_file:
            settings = Parameters(settings_file.read())

        settings["solver_settings"]["formulation"]["element_type"].SetString(self.element)
        settings["solver_settings"]["formulation"].AddEmptyValue("dynamic_tau")
        settings["solver_settings"]["formulation"]["dynamic_tau"].SetDouble(0.0)
        if self.element == "qsvms":
            settings["solver_settings"]["formulation"].AddEmptyValue("use_orthogonal_subscales")
            settings["solver_settings"]["formulation"]["use_orthogonal_subscales"].SetBool(self.oss_switch)
            settings["solver_settings"]["formulation"].AddEmptyValue("element_manages_time_integration")
            settings["solver_settings"]["formulation"]["element_manages_time_integration"].SetBool(True)
        elif self.element == "symbolic":
            settings["solver_settings"]["formulation"].AddEmptyValue("sound_velocity")
            settings["solver_settings"]["formulation"]["sound_velocity"].SetDouble(1e12)

        settings["solver_settings"].AddEmptyValue("time_scheme")
        settings["solver_settings"]["time_scheme"].SetString("bdf2")

        self.fluid_solver = navier_stokes_solver.CreateSolver(self.model,settings["solver_settings"])
        self.fluid_solver.AddVariables()
        self.fluid_solver.ImportModelPart()
        self.fluid_solver.PrepareModelPart()
        self.fluid_solver.AddDofs()

        fluid_model_part_name = settings["solver_settings"]["model_part_name"].GetString() + "." + settings["solver_settings"]["volume_model_part_name"].GetString()
        self.fluid_model_part = self.model.GetModelPart(fluid_model_part_name)

    def setUpProblem(self):
        xmin = 0.0
        xmax = 1.0
        ymin = 0.0
        ymax = 1.0

        ux = 1.0

        if self.element == "symbolic":
            ## Set up consitutive law
            rho = 1.0
            mu = 0.01
            self.fluid_model_part.Properties[1].SetValue(DENSITY,rho)
            self.fluid_model_part.Properties[1].SetValue(DYNAMIC_VISCOSITY,mu)
            constitutive_law = Newtonian2DLaw()
            self.fluid_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW,constitutive_law)

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

        self.fluid_solver.Initialize()

        for step in range(self.nsteps):
            time = self.fluid_solver.AdvanceInTime(time)

            self.fluid_solver.InitializeSolutionStep()
            self.fluid_solver.Predict()
            self.fluid_solver.SolveSolutionStep()
            self.fluid_solver.FinalizeSolutionStep()

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
        self.gid_io.WriteNodalResults(VISCOSITY,self.fluid_model_part.Nodes,label,0)
        self.gid_io.WriteNodalResults(DENSITY,self.fluid_model_part.Nodes,label,0)
        self.gid_io.PrintOnGaussPoints(SUBSCALE_VELOCITY,self.fluid_model_part,label)
        self.gid_io.PrintOnGaussPoints(SUBSCALE_PRESSURE,self.fluid_model_part,label)

    def __FinializeOutput(self):
        self.gid_io.FinalizeResults()

if __name__ == '__main__':
    UnitTest.main()
