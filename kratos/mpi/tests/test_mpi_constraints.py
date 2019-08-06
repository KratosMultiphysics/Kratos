from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import os
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils
data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
import KratosMultiphysics.mpi as mpi #TODO remove once the test script accepts a command-line argument.

structural_mechanics_is_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
if structural_mechanics_is_available:
    import KratosMultiphysics.StructuralMechanicsApplication


def ReadModelPart(model_part, mdpa_file_name):
    import_flags = KratosMultiphysics.ModelPartIO.READ | KratosMultiphysics.ModelPartIO.SKIP_TIMER

    KratosMultiphysics.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

def ReadDistributedModelPart(model_part, mdpa_file_name):
    from KratosMultiphysics.mpi import distributed_import_model_part_utility
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

    importer_settings = KratosMultiphysics.Parameters("""{
        "model_import_settings": {
            "input_type": "mdpa",
            "input_filename": \"""" + mdpa_file_name + """\",
            "partition_in_memory" : true
        },
        "echo_level" : 0
    }""")

    model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(model_part, importer_settings)
    model_part_import_util.ImportModelPart()
    model_part_import_util.CreateCommunicators()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestLinearMultipointConstraints(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_variables(self):
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

    def _apply_material_properties(self, dim):
        #define properties
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS, 210e9)
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.THICKNESS, 1.0)
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 1.0)

        g = [0, 0, 0]
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,
                                       g)

        if dim == 2:
            cl = KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw(
            )
        else:
            cl = KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw(
            )
        self.mp.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, cl)        

    def _add_dofs(self):
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,
            self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,
            self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,
            self.mp)

        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_X, self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_Y, self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.VELOCITY_Z, self.mp)

        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_X, self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_Y, self.mp)
        KratosMultiphysics.VariableUtils().AddDof(
            KratosMultiphysics.ACCELERATION_Z, self.mp)

    def _apply_BCs(self):
        bcs = self.mp.GetSubModelPart("DISPLACEMENT_fixed_edge")
        KratosMultiphysics.VariableUtils().SetScalarVar(
            KratosMultiphysics.DISPLACEMENT_X, 0.0, bcs.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(
            KratosMultiphysics.DISPLACEMENT_Y, 0.0, bcs.Nodes)

        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_X, True, bcs.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_Y, True, bcs.Nodes)

        bcmn = self.mp.GetSubModelPart("DISPLACEMENT_moving_node")
        KratosMultiphysics.VariableUtils().SetScalarVar(
            KratosMultiphysics.DISPLACEMENT_X, 0.01, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(
            KratosMultiphysics.DISPLACEMENT_Y, 0.0, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_X, True, bcmn.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(
            KratosMultiphysics.DISPLACEMENT_Y, True, bcmn.Nodes)

    def _setup_solver(self):
        # #define a minimal newton raphson solver
        # self.linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        # self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(
        #     self.linear_solver)
        # self.scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(
        #     -0.01)
        # self.convergence_criterion = KratosMultiphysics.ResidualCriteria(
        #     1e-10, 1e-12)
        # self.convergence_criterion.SetEchoLevel(0)

        # max_iters = 100
        # compute_reactions = False
        # reform_step_dofs = True
        # move_mesh_flag = False
        # self.strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
        #     self.mp, self.scheme, self.linear_solver, self.convergence_criterion,
        #     self.builder_and_solver, max_iters, compute_reactions,
        #     reform_step_dofs, move_mesh_flag)
        # self.strategy.SetEchoLevel(0)
        # self.strategy.Initialize()

        # self.strategy.Check()
        pass #TODO: use all the trilinos app once they are implemented

    def _reset(self):
        # del self.strategy
        # del self.linear_solver
        # del self.builder_and_solver
        # del self.scheme
        # del self.convergence_criterion
        pass

    def _solve(self):
        #self.strategy.Solve()
        pass

    def _check_results(self):
        constrained_nodes = self.mp.GetSubModelPart("DISPLACEMENT_moving_edge").Nodes # these nodes are slaves
        for node in constrained_nodes:
            print("x : ", node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X))
            print("y : ", node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y))

    def _setup_model_part(self):
        mdpa_file_name = "mdpa_files_for_tests/mpi_constraints_test"
        self.input_file     = GetFilePath(mdpa_file_name)
        if data_comm.IsDistributed():
            ReadDistributedModelPart(self.mp, self.input_file)
        else:
            ReadModelPart(self.mp, self.input_file)
        self.mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 2 # needed for the partitioner!

    def _apply_mpc_constraints(self):
        fixed_node = self.mp.GetSubModelPart("DISPLACEMENT_moving_node").Nodes[1] # Has 0.1 x disp (master node)
        constrained_nodes = self.mp.GetSubModelPart("DISPLACEMENT_moving_edge").Nodes # these nodes are slaves
        constraint_id = 1
        for node in constrained_nodes:
            self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, fixed_node, KratosMultiphysics.DISPLACEMENT_X, node, KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)
            constraint_id += 1
            self.mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", constraint_id, fixed_node, KratosMultiphysics.DISPLACEMENT_Y, node, KratosMultiphysics.DISPLACEMENT_Y, 1.0, 0)
            constraint_id += 1

    def _apply_mpi_constraints(self):
        fixed_node_num = 1
        mpi_constraint_utility = mpi.MpiConstraintsUtility(self.mp)
        constrained_nodes = self.mp.GetSubModelPart("DISPLACEMENT_moving_edge").Nodes # these nodes are slaves
        # rVariable, ConstraintId, SlaveNodeId, MasterNodeId, Weight, Constant
        constraint_id = 1
        for node in constrained_nodes:
            mpi_constraint_utility.AddConstraint(KratosMultiphysics.DISPLACEMENT_X, constraint_id, node.Id, fixed_node_num, 1.0, 0)
            constraint_id += 1
            mpi_constraint_utility.AddConstraint(KratosMultiphysics.DISPLACEMENT_Y, constraint_id, node.Id, fixed_node_num, 1.0, 0)
            constraint_id += 1
        mpi_constraint_utility.SynchronizeAndCreateConstraints()

      
    def _set_and_fill_buffer(self, buffer_size, delta_time):
        # Set buffer size
        self.mp.SetBufferSize(buffer_size)
        # Fill buffer
        time = self.mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        self.mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size - 1)
            self.mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            self.mp.CloneTimeStep(time)

        self.mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    @KratosUnittest.skipUnless(structural_mechanics_is_available,"StructuralMechanicsApplication is not available")
    def test_MPI_MPC_Constraints(self):
        dim = 2
        current_model = KratosMultiphysics.Model()
        self.mp= current_model.CreateModelPart("MainModelPart")
        self._add_variables()
        self._apply_material_properties(dim)
        self._setup_model_part()
        self._add_dofs()

        #time integration parameters
        dt = 0.01
        time = 0.0
        end_time = 0.01
        step = 0

        self._set_and_fill_buffer(2, dt)
        # Applying boundary conditions
        self._apply_BCs()
        # Applying constraints
        #self._apply_mpc_constraints()
        self._apply_mpi_constraints()
        # Solving the system of equations
        self._setup_solver()

        while (time <= end_time):
            time = time + dt
            step = step + 1
            self.mp.CloneTimeStep(time)
            self._solve()
        # Checking the results
        self._check_results()
        self._reset()

if __name__ == '__main__':
    KratosUnittest.main()
