from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
try:
    import KratosMultiphysics.FSIApplication as FSIApplication
    import NonConformant_OneSideMap as ncosm
    have_fsi = True
except ImportError as e:
    have_fsi = False

import KratosMultiphysics.KratosUnittest as UnitTest

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(have_fsi,"Missing required application: FSIApplication")
class ThermalCouplingTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.left_input_file = "thermal_coupling_left"
        self.right_input_file = "thermal_coupling_right"
        self.reference_file = "coupling_reference"

        self.dt = 0.01 #1e10 # This is steady state test
        self.nsteps = 30
        self.theta = 1.0 # Note: I don't think the face condition supports Crank-Nicolson (JC)
        # Note: Crank-Nicolson (theta=0.5) won't converge in a single iteration (or at all, for huge dt)
        self.num_coupling_iterations = 10
        self.temperature_relaxation_factor = 0.7
        self.coupling_relative_tolerance = 1e-5

        self.rho = 1.0
        self.c = 1.0
        self.k = 1.0

        self.ux = 0.0
        self.source = 0.0

        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False
        self.calculate_reactions = True

    def tearDown(self):
        self.deleteOutFile(self.left_input_file+'.time')
        self.deleteOutFile(self.right_input_file+'.time')

    def deleteOutFile(self,filename):
        with WorkFolderScope("ThermalCouplingTest"):
            try:
                os.remove(filename)
            except FileNotFoundError as e:
                pass

    def testDirichletNeumann(self):
        current_model = KratosMultiphysics.Model()
        with WorkFolderScope("ThermalCouplingTest"):
            self.setUpModel(current_model)
            self.setUpSolvers()
            self.setUpMapper()

            self.setUpOuterBoundaryCondition(self.left_model_part,0.0,303.15)
            self.setUpOuterBoundaryCondition(self.right_model_part,1.0,293.15)
            self.setUpDirichletCouplingBoundary(self.right_model_part)
            #self.setUpDirichletCouplingBoundary(self.left_model_part)

            if self.print_output:
                num_left_nodes = len(self.left_model_part.Nodes)
                for node in self.right_model_part.Nodes:
                    node.Id = node.Id + num_left_nodes
                self.InitializeOutput()

            self.runTest()

            self.checkResults(self.mapper.str_interface.Nodes,KratosMultiphysics.TEMPERATURE)
            if self.print_output:
                self.FinalizeOutput()


    def setUpModel(self,current_model):

        self.left_model_part = current_model.CreateModelPart("LeftSide")
        self.right_model_part = current_model.CreateModelPart("RightSide")

        thermal_settings = KratosMultiphysics.ConvectionDiffusionSettings()
        thermal_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)
        thermal_settings.SetDensityVariable(KratosMultiphysics.DENSITY)
        thermal_settings.SetSpecificHeatVariable(KratosMultiphysics.SPECIFIC_HEAT)
        thermal_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)
        thermal_settings.SetVolumeSourceVariable(KratosMultiphysics.HEAT_FLUX)
        thermal_settings.SetSurfaceSourceVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        thermal_settings.SetVelocityVariable(KratosMultiphysics.VELOCITY)
        thermal_settings.SetMeshVelocityVariable(KratosMultiphysics.MESH_VELOCITY)
        #thermal_settings.SetProjectionVariable(KratosMultiphysics.PROJECTED_SCALAR1)
        thermal_settings.SetReactionVariable(KratosMultiphysics.REACTION_FLUX)

        self.left_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS,thermal_settings)
        self.right_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS,thermal_settings)

    def setUpSolvers(self):

        import convection_diffusion_solver as thermal_solver
        thermal_solver.AddVariables(self.left_model_part)
        thermal_solver.AddVariables(self.right_model_part)
        # Also add mapper variables
        ncosm.AddVariables(self.left_model_part,self.right_model_part)
        # auxiliary container for reaction->distributed flux conversion
        self.right_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)
        # Temporary container for un-relaxed temperature
        self.left_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_PAUX)

        model_part_io_left = KratosMultiphysics.ModelPartIO(self.left_input_file)
        model_part_io_left.ReadModelPart(self.left_model_part)

        model_part_io_right = KratosMultiphysics.ModelPartIO(self.right_input_file)
        model_part_io_right.ReadModelPart(self.right_model_part)

        self.left_model_part.SetBufferSize(2)
        self.right_model_part.SetBufferSize(2)
        thermal_solver.AddDofs(self.left_model_part)
        thermal_solver.AddDofs(self.right_model_part)

        self.left_solver = thermal_solver.ConvectionDiffusionSolver(self.left_model_part,self.domain_size)
        self.left_solver.calculate_reactions = self.calculate_reactions
        self.left_solver.theta = self.theta
        self.left_solver.echo_level = 0
        self.left_solver.Initialize()

        self.right_solver = thermal_solver.ConvectionDiffusionSolver(self.right_model_part,self.domain_size)
        self.right_solver.calculate_reactions = self.calculate_reactions
        self.right_solver.theta = self.theta
        self.right_solver.echo_level = 0
        self.right_solver.Initialize()

    def setUpOuterBoundaryCondition(self,model_part,boundary_x,boundary_value):

        velocity = KratosMultiphysics.Array3()
        velocity[0] = self.ux
        velocity[1] = 0.0
        velocity[2] = 0.0

        ## Set initial and boundary conditions
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY,self.rho)
            node.SetSolutionStepValue(KratosMultiphysics.CONDUCTIVITY,self.k)
            node.SetSolutionStepValue(KratosMultiphysics.SPECIFIC_HEAT,self.c)
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,self.source)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,velocity)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,boundary_value)

            if node.X == boundary_x:
                node.Fix(KratosMultiphysics.TEMPERATURE)

    def setUpDirichletCouplingBoundary(self,model_part):
        for cond in model_part.Conditions:
            for node in cond.GetNodes():
                node.Fix(KratosMultiphysics.TEMPERATURE)

    def setUpMapper(self):
        for cond in self.left_model_part.Conditions:
            for node in cond.GetNodes():
                node.Set(KratosMultiphysics.INTERFACE,True)

        for cond in self.right_model_part.Conditions:
            for node in cond.GetNodes():
                node.Set(KratosMultiphysics.INTERFACE,True)

        self.mapper = ncosm.NonConformant_OneSideMap(self.left_model_part,self.right_model_part,
                 search_radius_factor=2.0, it_max=50, tol=1e-5)

    def runTest(self):
        time = 0.0

        for step in range(self.nsteps):
            time = time+self.dt
            self.left_model_part.CloneTimeStep(time)
            self.right_model_part.CloneTimeStep(time)

            iter = 0
            while iter < self.num_coupling_iterations:

                # Solve Dirichlet side -> Get reactions
                self.right_solver.Solve()

                # Map reactions
                KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(
                    self.mapper.str_interface, KratosMultiphysics.REACTION_FLUX, KratosMultiphysics.NODAL_PAUX, 1e-5, 50)
                self.mapper.StructureToFluid_ScalarMap(KratosMultiphysics.NODAL_PAUX,KratosMultiphysics.FACE_HEAT_FLUX,False)

                # Solve Neumann side
                self.left_solver.Solve()

                # Get updated temperature
                self.mapper.FluidToStructure_ScalarMap(KratosMultiphysics.TEMPERATURE,KratosMultiphysics.NODAL_PAUX,True)
                temperature_difference = 0.0
                for node in self.mapper.str_interface.Nodes:
                    old_temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                    new_temperature = node.GetSolutionStepValue(KratosMultiphysics.NODAL_PAUX)
                    interpolated_temperature = (1.0-self.temperature_relaxation_factor)*old_temperature + self.temperature_relaxation_factor*new_temperature
                    temperature_difference += (old_temperature-new_temperature)**2
                    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, interpolated_temperature )

                iter += 1
                if (temperature_difference**0.5)/len(self.mapper.str_interface.Nodes) <= self.coupling_relative_tolerance:
                    #print("Convergence after {0} coupling iterations".format(iter))
                    break

            if self.print_output:
                self.PrintOutput()


    def checkResults(self,nodes,checked_variable):

        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, {0}\n".format(checked_variable.Name()))
                for node in nodes:
                    value = node.GetSolutionStepValue(checked_variable,0)
                    ref_file.write("{0}, {1}\n".format(node.Id, value))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()

                for node in nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    node_id = values[0]
                    reference_value = values[1]

                    value = node.GetSolutionStepValue(checked_variable)
                    self.assertAlmostEqual(reference_value, value, delta=self.check_tolerance)

                    line = reference_file.readline()
                if line != '': # If we did not reach the end of the reference file
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")


    def InitializeOutput(self):
        gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
        multifile = KratosMultiphysics.MultiFileFlag.SingleFile
        deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
        self.left_gid_io = KratosMultiphysics.GidIO(self.left_input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)
        self.right_gid_io = KratosMultiphysics.GidIO(self.right_input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)

        mesh_name = 0.0
        self.left_gid_io.InitializeMesh( mesh_name)
        self.left_gid_io.WriteMesh( self.left_model_part.GetMesh() )
        self.left_gid_io.FinalizeMesh()
        self.left_gid_io.InitializeResults(mesh_name,(self.left_model_part).GetMesh())


        self.right_gid_io.InitializeMesh( mesh_name)
        self.right_gid_io.WriteMesh( self.right_model_part.GetMesh() )
        self.right_gid_io.FinalizeMesh()
        self.right_gid_io.InitializeResults(mesh_name,(self.right_model_part).GetMesh())

    def FinalizeOutput(self):
        self.left_gid_io.FinalizeResults()
        self.right_gid_io.FinalizeResults()

    def PrintOutput(self):
        for model_part,gid_io in [[self.left_model_part,self.left_gid_io],[self.right_model_part,self.right_gid_io]]:
            label = model_part.ProcessInfo[KratosMultiphysics.TIME]
            gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(KratosMultiphysics.TEMPERATURE,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(KratosMultiphysics.DENSITY,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(KratosMultiphysics.CONDUCTIVITY,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(KratosMultiphysics.SPECIFIC_HEAT,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(KratosMultiphysics.HEAT_FLUX,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(KratosMultiphysics.FACE_HEAT_FLUX,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(KratosMultiphysics.REACTION_FLUX,model_part.Nodes,label,0)

if __name__ == '__main__':
    if have_fsi:
        test = ThermalCouplingTest()
        test.setUp()
        #test.print_reference_values = True
        test.print_output = True
        test.testDirichletNeumann()
        test.tearDown()
    else:
        UnitTest.skip("Missing required application: FSIApplication")
