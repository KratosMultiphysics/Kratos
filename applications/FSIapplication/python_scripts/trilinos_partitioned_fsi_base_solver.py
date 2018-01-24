from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
import python_solvers_wrapper_structural       # Import the structure Python solvers wrapper
import convergence_accelerator_factory         # Import the FSI convergence accelerator factory

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications(
    "MetisApplication",
    "TrilinosApplication",
    "MappingApplication",
    "FSIApplication",
    "ALEApplication",
    "FluidDynamicsApplication",
    "StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.ALEApplication as KratosALE
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Import base class file
import partitioned_fsi_base_solver

def CreateSolver(structure_main_model_part, fluid_main_model_part, project_parameters):
    return TrilinosPartitionedFSIBaseSolver(structure_main_model_part, fluid_main_model_part, project_parameters)

class TrilinosPartitionedFSIBaseSolver(partitioned_fsi_base_solver.PartitionedFSIBaseSolver):
    def __init__(self, structure_main_model_part, fluid_main_model_part, project_parameters):

        if (KratosMPI.mpi.rank == 0) : print("** Calling the partitioned FSI Trilinos base solver constructor...")

        # Initial tests
        start_time_structure = project_parameters["structure_solver_settings"]["problem_data"]["start_time"].GetDouble()
        start_time_fluid = project_parameters["fluid_solver_settings"]["problem_data"]["start_time"].GetDouble()
        end_time_structure = project_parameters["structure_solver_settings"]["problem_data"]["end_time"].GetDouble()
        end_time_fluid = project_parameters["fluid_solver_settings"]["problem_data"]["end_time"].GetDouble()

        if start_time_structure != start_time_fluid:
            if (KratosMPI.mpi.rank == 0) : raise("ERROR: Different initial time among subdomains!")
        if end_time_structure != end_time_fluid:
            if (KratosMPI.mpi.rank == 0) : raise("ERROR: Different final time among subdomains!")

        self.structure_main_model_part = structure_main_model_part
        self.fluid_main_model_part = fluid_main_model_part

        # Time stepping checks (no sub-stepping between subdomains has been implemented yed)
        time_step_structure = project_parameters["structure_solver_settings"]["problem_data"]["time_step"].GetDouble()
        # If automatic time stepping has been selected in the fluid domain, deactivate it and use the structure time step
        if (project_parameters["fluid_solver_settings"]["solver_settings"]["time_stepping"]["automatic_time_step"].GetBool()):
            project_parameters["fluid_solver_settings"]["solver_settings"]["time_stepping"]["automatic_time_step"].SetBool(False)
            time_step_fluid = time_step_structure
            if (KratosMPI.mpi.rank == 0) : print("WARNING: Automatic fluid time stepping cannot be used. Setting structure time step as fluid time step.")
        else:
            time_step_fluid = project_parameters["fluid_solver_settings"]["solver_settings"]["time_stepping"]["time_step"].GetDouble()
            if time_step_structure != time_step_fluid:
                if (KratosMPI.mpi.rank == 0) : raise("ERROR: Different time step among subdomains! No sub-stepping implemented yet.")

        self.time_step = time_step_fluid

        # Take the each one of the solvers settings from the ProjectParameters
        # Note that the defaults check will be performed inside each field solver
        self.settings = KratosMultiphysics.Parameters("{}")
        self.settings.AddValue("structure_solver_settings",project_parameters["structure_solver_settings"]["solver_settings"])
        self.settings.AddValue("fluid_solver_settings",project_parameters["fluid_solver_settings"]["solver_settings"])
        self.settings.AddValue("coupling_solver_settings",project_parameters["coupling_solver_settings"]["solver_settings"])
        self.settings.AddValue("mapper_settings",project_parameters["coupling_solver_settings"]["mapper_settings"])

        # Auxiliar variables
        self.max_nl_it = self.settings["coupling_solver_settings"]["nl_max_it"].GetInt()
        self.nl_tol = self.settings["coupling_solver_settings"]["nl_tol"].GetDouble()
        self.solve_mesh_at_each_iteration = self.settings["coupling_solver_settings"]["solve_mesh_at_each_iteration"].GetBool()
        self.coupling_algorithm = self.settings["coupling_solver_settings"]["coupling_scheme"].GetString()
        self.fluid_interface_submodelpart_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][0].GetString()
        self.structure_interface_submodelpart_name = self.settings["coupling_solver_settings"]["structure_interfaces_list"][0].GetString()

        # Construct the structure solver
        self.structure_solver = python_solvers_wrapper_structural.CreateSolver(self.structure_main_model_part,
                                                                               project_parameters["structure_solver_settings"])
        if (KratosMPI.mpi.rank == 0) : print("* Structure solver constructed.")

        # Construct the fluid solver
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolver(self.fluid_main_model_part,
                                                                      project_parameters["fluid_solver_settings"])
        if (KratosMPI.mpi.rank == 0) : print("* Fluid solver constructed.")

        # Construct the ALE mesh solver
        mesh_solver_settings = KratosMultiphysics.Parameters("{}")
        self.mesh_solver_module = __import__(self.settings["coupling_solver_settings"]["mesh_solver"].GetString())
        self.mesh_solver = self.mesh_solver_module.CreateSolver(self.fluid_solver.main_model_part,
                                                                mesh_solver_settings)
        if (KratosMPI.mpi.rank == 0):
            print("* ALE mesh solver constructed.")
            print("** Partitioned FSI base solver constructed.")


    def AddVariables(self):
        ## Structure variables addition
        # Standard CSM variables addition
        self.structure_solver.AddVariables()

        ## Fluid variables addition
        # Standard CFD variables addition
        self.fluid_solver.AddVariables()
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_ACCELERATION) # TODO: This should be added in the mesh solvers
        # Mesh solver variables addition
        self.mesh_solver.AddVariables()

        ## FSIApplication variables addition
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_RESIDUAL)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_MESH_RESIDUAL)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)


    #######################################################################
    ##############          PRIVATE METHODS SECTION          ##############
    #######################################################################

    ### AUXILIAR METHODS ###
    def _GetPartitionedFSIUtilities(self):

        if (self.domain_size == 2):
            return KratosTrilinos.TrilinosPartitionedFSIUtilities2D(self.epetra_communicator)
        else:
            return KratosTrilinos.TrilinosPartitionedFSIUtilities3D(self.epetra_communicator)

    ### TODO: SUBSTITUTE IN THIS METHOD THE OLD MAPPER BY THE ONE IN THE FSI APPLICATION
    def _SetUpMapper(self):


        if (self.settings["mapper_settings"].size() == 1):
            fluid_submodelpart_name = self.settings["mapper_settings"][0]["fluid_interface_submodelpart_name"].GetString()
            structure_submodelpart_name = self.settings["mapper_settings"][0]["structure_interface_submodelpart_name"].GetString()

            project_parameters_mapper = KratosMultiphysics.Parameters("{}")
            project_parameters_mapper.AddEmptyValue("mapper_type").SetString("NearestElement")
            project_parameters_mapper.AddEmptyValue("interface_submodel_part_origin").SetString(fluid_submodelpart_name)
            project_parameters_mapper.AddEmptyValue("interface_submodel_part_destination").SetString(structure_submodelpart_name)

            self.interface_mapper = KratosMapping.MapperFactory(self.fluid_solver.main_model_part,
				                                                self.structure_solver.main_model_part,
				                                                project_parameters_mapper)

            self.double_faced_structure = False

        elif (self.settings["mapper_settings"].size() == 2):
            # Get the fluid interface faces submodelpart names
            for mapper_id in range(2):
                if (self.settings["mapper_settings"][mapper_id]["mapper_face"].GetString() == "Positive"):
                    pos_face_submodelpart_name = self.settings["mapper_settings"][mapper_id]["fluid_interface_submodelpart_name"].GetString()
                elif (self.settings["mapper_settings"][mapper_id]["mapper_face"].GetString() == "Negative"):
                    neg_face_submodelpart_name = self.settings["mapper_settings"][mapper_id]["fluid_interface_submodelpart_name"].GetString()
                else:
                    raise Exception("Unique mapper flag has been set but more than one mapper exist in mapper_settings.")
            # Get the structure submodelpart name
            structure_submodelpart_name = self.settings["mapper_settings"][0]["structure_interface_submodelpart_name"].GetString()

            # Grab the interface submodelparts
            pos_fluid_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(pos_face_submodelpart_name)
            neg_fluid_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(neg_face_submodelpart_name)
            structure_submodelpart = self.structure_solver.main_model_part.GetSubModelPart(structure_submodelpart_name)

            # TODO: SET THE DOUBLE SIDE SURFACE MAPPING

            # search_radius_factor = 2.0
            # mapper_max_iterations = 200
            # mapper_tolerance = 1e-12

            # self.interface_mapper = NonConformant_OneSideMap.NonConformantTwoFaces_OneSideMap(pos_fluid_submodelpart,
            #                                                                                   neg_fluid_submodelpart,
            #                                                                                   structure_submodelpart,
            #                                                                                   search_radius_factor,
            #                                                                                   mapper_max_iterations,
            #                                                                                   mapper_tolerance)

            # TODO: SET THE DOUBLE SIDE SURFACE MAPPING

            self.double_faced_structure = True

        else:
            raise Exception("Case with more than 2 mappers has not been implemented yet.\n \
                             Please, in case you are using single faced immersed bodies, set the skin entities in a unique submodelpart.\n \
                             In case you are considering double faced immersed bodies (shells or membranes), set all the positive faces \
                             in a unique submodelpart and all the negative ones in another submodelpart.")


    def _SetStructureNeumannCondition(self):

        structure_computational_submodelpart = self.structure_solver.GetComputingModelPart()

        # Get the maximum condition id
        max_cond_id = 0
        for condition in self.structure_solver.main_model_part.Conditions:
            max_cond_id = max(max_cond_id, condition.Id)

        max_cond_id = self.structure_solver.main_model_part.GetCommunicator().MaxAll(max_cond_id)

        # Set up the point load condition in the structure interface
        for i in range(self.settings["coupling_solver_settings"]["structure_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["structure_interfaces_list"][i].GetString()
            interface_submodelpart_i = self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name)

            # Get the number of conditions to be set in each processor
            local_nodes_number_accumulated = -1
            local_nodes_number = len(interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes)
            local_nodes_number_accumulated = interface_submodelpart_i.GetCommunicator().ScanSum(local_nodes_number, local_nodes_number_accumulated)

            # Create the point load condition
            aux_count = max_cond_id + local_nodes_number_accumulated
            if self.domain_size == 2:
                for node in interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes:
                    aux_count+=1
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition2D1N", int(aux_count), [node.Id], self.structure_solver.main_model_part.Properties[0])

            elif self.domain_size == 3:
                for node in interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes:
                    aux_count+=1
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition3D1N", int(aux_count), [node.Id], self.structure_solver.main_model_part.Properties[0])


    # TODO: GET IT FROM THE SERIAL BASE SOLVER ONCE THE MAPPER IN MAPPING APPLICATION IS USED IN SERIAL PROBLEMS AS WELL
    def _ComputeMeshPredictionSingleFaced(self):

            if (KratosMPI.mpi.rank == 0) : print("Computing time step ",self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]," prediction...")
            # Get the previous step fluid interface nodal fluxes
            self.interface_mapper.Map(KratosMultiphysics.REACTION,
                                      KratosStructural.POINT_LOAD,
                                      KratosMapping.MapperFactory.SWAP_SIGN |
                                      KratosMapping.MapperFactory.CONSERVATIVE)

            # Solve the current step structure problem with the previous step fluid interface nodal fluxes
            self.structure_solver.SolveSolutionStep()

            # Map the obtained structure displacement to the fluid interface
            self.interface_mapper.InverseMap(KratosMultiphysics.MESH_DISPLACEMENT,
                                             KratosMultiphysics.DISPLACEMENT)

            # Solve the mesh problem
            self.mesh_solver.Solve()

            if (KratosMPI.mpi.rank == 0) : print("Mesh prediction computed.")
