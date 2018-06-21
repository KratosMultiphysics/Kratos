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
    "MeshMovingApplication",
    "FluidDynamicsApplication",
    "StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Import base class file
import partitioned_fsi_base_solver

def CreateSolver(model, project_parameters):
    return TrilinosPartitionedFSIBaseSolver(model, project_parameters)

class TrilinosPartitionedFSIBaseSolver(partitioned_fsi_base_solver.PartitionedFSIBaseSolver):
    def __init__(self, model, project_parameters):
        super(TrilinosPartitionedFSIBaseSolver, self).__init__(model, project_parameters)

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
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VAUX_EQ_TRACTION)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_RESIDUAL)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_MESH_RESIDUAL)

        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VAUX_EQ_TRACTION)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE)

    def _PrintInfoOnRankZero(self, *args):
        KratosMPI.mpi.world.barrier()
        if KratosMPI.mpi.rank == 0:
            KratosMultiphysics.Logger.PrintInfo(" ".join(map(str, args)))

    def _PrintWarningOnRankZero(self, *args):
        KratosMPI.mpi.world.barrier()
        if KratosMPI.mpi.rank == 0:
            KratosMultiphysics.Logger.PrintWarning(" ".join(map(str, args)))

    def _GetPartitionedFSIUtilities(self):
        if (self.domain_size == 2):
            return KratosTrilinos.TrilinosPartitionedFSIUtilities2D(self.epetra_communicator)
        else:
            return KratosTrilinos.TrilinosPartitionedFSIUtilities3D(self.epetra_communicator)

    ### TODO: SUBSTITUTE IN THIS METHOD THE OLD MAPPER BY THE ONE IN THE FSI APPLICATION
    def _SetUpMapper(self):
        mapper_settings = self.settings["coupling_solver_settings"]["mapper_settings"]
        if (mapper_settings.size() == 1):
            fluid_submodelpart_name = mapper_settings[0]["fluid_interface_submodelpart_name"].GetString()
            structure_submodelpart_name = mapper_settings[0]["structure_interface_submodelpart_name"].GetString()

            mapper_project_parameters = KratosMultiphysics.Parameters("""{
                "mapper_type" : "",
                "interface_submodel_part_origin" : "",
                "interface_submodel_part_destination" : ""
            }""")
            mapper_project_parameters["mapper_type"].SetString("nearest_element")
            mapper_project_parameters["interface_submodel_part_origin"].SetString(fluid_submodelpart_name)
            mapper_project_parameters["interface_submodel_part_destination"].SetString(structure_submodelpart_name)

            self.interface_mapper = KratosMapping.MapperFactory.CreateMapper(self.fluid_solver.main_model_part,
                                                                             self.structure_solver.main_model_part,
                                                                             mapper_project_parameters)

            self.double_faced_structure = False

        elif (mapper_settings.size() == 2):
            # Get the fluid interface faces submodelpart names
            for mapper_id in range(2):
                if (mapper_settings[mapper_id]["mapper_face"].GetString() == "Positive"):
                    pos_face_submodelpart_name = mapper_settings[mapper_id]["fluid_interface_submodelpart_name"].GetString()
                elif (mapper_settings[mapper_id]["mapper_face"].GetString() == "Negative"):
                    neg_face_submodelpart_name = mapper_settings[mapper_id]["fluid_interface_submodelpart_name"].GetString()
                else:
                    raise Exception("Unique mapper flag has been set but more than one mapper exist in mapper_settings.")
            # Get the structure submodelpart name
            structure_submodelpart_name = mapper_settings[0]["structure_interface_submodelpart_name"].GetString()

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

    # TODO: GET IT FROM THE SERIAL BASE SOLVER ONCE THE MAPPER IN MAPPING APPLICATION IS USED IN SERIAL PROBLEMS AS WELL
    def _ComputeMeshPredictionSingleFaced(self):
        step = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        self._PrintInfoOnRankZero("Computing time step ",str(step)," prediction...")
        
        # Set the redistribution settings
        redistribution_tolerance = 1e-8
        redistribution_max_iters = 50

        # Convert the nodal reaction to traction loads before transfering
        KratosFSI.VariableRedistributionUtility.DistributePointValues(
            self._GetFluidInterfaceSubmodelPart(),
            KratosMultiphysics.REACTION,
            KratosMultiphysics.VAUX_EQ_TRACTION,
            redistribution_tolerance,
            redistribution_max_iters)

        # Transfer fluid tractions to the structure interface
        self.interface_mapper.Map(KratosMultiphysics.VAUX_EQ_TRACTION,
                                  KratosMultiphysics.VAUX_EQ_TRACTION,
                                  KratosMapping.Mapper.SWAP_SIGN)

        # Convert the transferred traction loads to point loads
        KratosFSI.VariableRedistributionUtility.ConvertDistributedValuesToPoint(
            self._GetStructureInterfaceSubmodelPart(),
            KratosMultiphysics.VAUX_EQ_TRACTION,
            KratosStructural.POINT_LOAD)

        # Solve the current step structure problem with the previous step fluid interface nodal fluxes
        is_converged = self.structure_solver.SolveSolutionStep()
        if not is_converged:
            self._PrintWarningOnRankZero("Mesh prediction structure solver did not converge.")
        
        # Map the obtained structure displacement to the fluid interface
        self.interface_mapper.InverseMap(KratosMultiphysics.MESH_DISPLACEMENT, KratosMultiphysics.DISPLACEMENT)

        # Solve the mesh problem
        self.mesh_solver.Solve()

        self._PrintInfoOnRankZero("Mesh prediction computed.")

    def _ComputeMeshPredictionDoubleFaced(self):
        raise Exception("Trilinos _ComputeMeshPredictionDoubleFaced() is to be implemented!")
