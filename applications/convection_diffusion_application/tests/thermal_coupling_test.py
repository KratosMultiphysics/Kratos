from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
try:
    from KratosMultiphysics.FSIApplication import *
    import NonConformant_OneSideMap as ncosm
    have_fsi = True
except ImportError as e:
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
class ThermalCouplingTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.left_input_file = "thermal_coupling_left"
        self.right_input_file = "thermal_coupling_right"
        self.reference_file = "coupling_reference"

        self.dt = 0.01 #1e10 # This is steady state test
        self.nsteps = 50
        self.theta = 1.0 # Note: I don't think the face condition supports Crank-Nicolson (JC)
        # Note: Crank-Nicolson (theta=0.5) won't converge in a single iteration (or at all, for huge dt)
        self.num_coupling_iterations = 5

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
                remove(filename)
            except FileNotFoundError as e:
                pass

    def testDirichletNeumann(self):

        with WorkFolderScope("ThermalCouplingTest"):
            self.setUpModel()
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

            #self.checkResults(TEMPERATURE)
            if self.print_output:
                self.FinalizeOutput()


    def setUpModel(self):

        self.left_model_part = ModelPart("LeftSide")
        self.right_model_part = ModelPart("RightSide")

        thermal_settings = ConvectionDiffusionSettings()
        thermal_settings.SetUnknownVariable(TEMPERATURE)
        thermal_settings.SetDensityVariable(DENSITY)
        thermal_settings.SetSpecificHeatVariable(SPECIFIC_HEAT)
        thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
        thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
        thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
        thermal_settings.SetVelocityVariable(VELOCITY)
        thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
        #thermal_settings.SetProjectionVariable(PROJECTED_SCALAR1)
        thermal_settings.SetReactionVariable(REACTION_FLUX)

        self.left_model_part.ProcessInfo.SetValue(CONVECTION_DIFFUSION_SETTINGS,thermal_settings)
        self.right_model_part.ProcessInfo.SetValue(CONVECTION_DIFFUSION_SETTINGS,thermal_settings)

    def setUpSolvers(self):

        import convection_diffusion_solver as thermal_solver
        thermal_solver.AddVariables(self.left_model_part)
        thermal_solver.AddVariables(self.right_model_part)
        # Also add mapper variables
        ncosm.AddVariables(self.left_model_part,self.right_model_part)
        # auxiliary container for reaction->distributed flux conversion
        self.right_model_part.AddNodalSolutionStepVariable(NODAL_PAUX)

        model_part_io_left = ModelPartIO(self.left_input_file)
        model_part_io_left.ReadModelPart(self.left_model_part)

        model_part_io_right = ModelPartIO(self.right_input_file)
        model_part_io_right.ReadModelPart(self.right_model_part)

        self.left_model_part.SetBufferSize(2)
        self.right_model_part.SetBufferSize(2)
        thermal_solver.AddDofs(self.left_model_part)
        thermal_solver.AddDofs(self.right_model_part)

        self.left_solver = thermal_solver.ConvectionDiffusionSolver(self.left_model_part,self.domain_size)
        self.left_solver.calculate_reactions = self.calculate_reactions
        self.left_solver.theta = self.theta
        self.left_solver.Initialize()

        self.right_solver = thermal_solver.ConvectionDiffusionSolver(self.right_model_part,self.domain_size)
        self.right_solver.calculate_reactions = self.calculate_reactions
        self.right_solver.theta = self.theta
        self.right_solver.Initialize()

    def setUpOuterBoundaryCondition(self,model_part,boundary_x,boundary_value):

        velocity = Array3()
        velocity[0] = self.ux
        velocity[1] = 0.0
        velocity[2] = 0.0

        ## Set initial and boundary conditions
        for node in model_part.Nodes:
            node.SetSolutionStepValue(DENSITY,self.rho)
            node.SetSolutionStepValue(CONDUCTIVITY,self.k)
            node.SetSolutionStepValue(SPECIFIC_HEAT,self.c)
            node.SetSolutionStepValue(HEAT_FLUX,self.source)
            node.SetSolutionStepValue(VELOCITY,velocity)
            node.SetSolutionStepValue(TEMPERATURE,boundary_value)

            if node.X == boundary_x:
                node.Fix(TEMPERATURE)
    
    def setUpDirichletCouplingBoundary(self,model_part):
        for cond in model_part.Conditions:
            for node in cond.GetNodes():
                node.Fix(TEMPERATURE)

    def setUpMapper(self):
        for cond in self.left_model_part.Conditions:
            for node in cond.GetNodes():
                node.Set(INTERFACE,True)

        for cond in self.right_model_part.Conditions:
            for node in cond.GetNodes():
                node.Set(INTERFACE,True)

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
                VariableRedistributionUtility.DistributePointValues(
                    self.mapper.str_interface, REACTION_FLUX, NODAL_PAUX, 1e-5, 50)
                self.mapper.StructureToFluid_ScalarMap(NODAL_PAUX,FACE_HEAT_FLUX,False)

                # Solve Neumann side
                self.left_solver.Solve()

                # Get updated temperature
                self.mapper.FluidToStructure_ScalarMap(TEMPERATURE,TEMPERATURE,True)
                for node in self.mapper.str_interface.Nodes:
                    node.SetSolutionStepValue(TEMPERATURE,1,node.GetSolutionStepValue(TEMPERATURE))

                if self.print_output:
                    self.PrintOutput()

                iter += 1


    def checkResults(self,checked_variable):

        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, {0}\n".format(checked_variable.Name))
                for node in self.model_part.Nodes:
                    value = node.GetSolutionStepValue(checked_variable,0)
                    ref_file.write("{0}, {1}\n".format(node.Id, value))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()
                node_iter = self.model_part.Nodes

                for node in self.model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    node_id = values[0]
                    reference_value = values[1]

                    value = node.GetSolutionStepValue(checked_variable)
                    self.assertAlmostEqual(reference_value, value, delta=self.check_tolerance)

                    line = reference_file.readline()
                if line != '': # If we did not reach the end of the reference file
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")


    def InitializeOutput(self):
        gid_mode = GiDPostMode.GiD_PostBinary
        multifile = MultiFileFlag.SingleFile
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteElementsOnly
        self.left_gid_io = GidIO(self.left_input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)
        self.right_gid_io = GidIO(self.right_input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)

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
            label = model_part.ProcessInfo[TIME]
            gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(DENSITY,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(CONDUCTIVITY,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(SPECIFIC_HEAT,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(HEAT_FLUX,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(FACE_HEAT_FLUX,model_part.Nodes,label,0)
            gid_io.WriteNodalResults(REACTION_FLUX,model_part.Nodes,label,0)

if __name__ == '__main__':
    test = ThermalCouplingTest()
    test.setUp()
    test.print_reference_values = True
    test.print_output = True
    test.testDirichletNeumann()
    test.tearDown()