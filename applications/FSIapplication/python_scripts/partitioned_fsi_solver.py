from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
import NonConformant_OneSideMap                # Import non-conformant mapper
import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as KratosALE
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(structure_main_model_part, fluid_main_model_part, project_parameters):
    return PartitionedFSISolver(structure_main_model_part, fluid_main_model_part, project_parameters)


class PartitionedFSISolver:
    def __init__(self, structure_main_model_part, fluid_main_model_part, project_parameters):

        # Initial tests
        start_time_structure = project_parameters["structure_solver_settings"]["problem_data"]["start_time"].GetDouble()
        start_time_fluid = project_parameters["fluid_solver_settings"]["problem_data"]["start_step"].GetDouble()
        end_time_structure = project_parameters["structure_solver_settings"]["problem_data"]["end_time"].GetDouble()
        end_time_fluid = project_parameters["fluid_solver_settings"]["problem_data"]["end_time"].GetDouble()

        if start_time_structure != start_time_fluid:
            raise("ERROR: Different initial time among subdomains!")
        if end_time_structure != end_time_fluid:
            raise("ERROR: Different final time among subdomains!")

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.structure_main_model_part = structure_main_model_part
        self.fluid_main_model_part = fluid_main_model_part

        # Settings string in JSON format
        default_settings = KratosMultiphysics.Parameters("""
        {
        "structure_solver_settings":
            {
            "solver_type": "solid_mechanics_implicit_dynamic_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "echo_level": 0,
            "time_integration_method": "Implicit",
            "analysis_type": "nonlinear",
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "solution_type": "Dynamic",
            "scheme_type": "Newmark",
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance" : 1.0e-3,
            "displacement_absolute_tolerance" : 1.0e-5,
            "residual_relative_tolerance"     : 1.0e-3,
            "residual_absolute_tolerance"     : 1.0e-5,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type"   : "SuperLUSolver",
                "max_iteration" : 200,
                "tolerance"     : 1e-7,
                "scaling"       : false,
                "verbosity"     : 1
            },
            "processes_sub_model_part_list": [""],
            "problem_domain_sub_model_part_list": ["solid_model_part"]
            },
        "fluid_solver_settings":
            {
            "solver_type": "navier_stokes_solver_vmsmonolithic",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "dynamic_tau" : 0.0,
            "oss_switch"  : 0,
            "echo_level"  : 0,
            "consider_periodic_conditions" : false,
            "compute_reactions"            : true,
            "divergence_clearance_steps"   : 0,
            "reform_dofs_at_each_step"     : true,
            "relative_velocity_tolerance"  : 1e-3,
            "absolute_velocity_tolerance"  : 1e-5,
            "relative_pressure_tolerance"  : 1e-3,
            "absolute_pressure_tolerance"  : 1e-5,
            "linear_solver_settings"        : {
                "solver_type"         : "AMGCL",
                "max_iteration"       : 200,
                "tolerance"           : 1e-9,
                "provide_coordinates" : true,
                "smoother_type"       : "ilu0",
                "krylov_type"         : "gmres",
                "coarsening_type"     : "aggregation",
                "scaling"             : true,
                "verbosity"           : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts"             : [""],
            "no_skin_parts"          : [""],
            "time_stepping"          : {
                "automatic_time_step" : false,
                "time_step"           : 0.1
            },
            "alpha"                  :-0.3,
            "move_mesh_strategy"     : 0,
            "periodic"               : "periodic",
            "move_mesh_flag"         : false,
            "turbulence_model"       : "None"
            },
        "coupling_solver_settings":
            {
            "coupling_scheme"                : "DirichletNeumann",
            "solver_type"                    : "partitioned_fsi_solver",
            "nl_tol"                         : 1e-5,
            "nl_max_it"                      : 50,
            "move_interface"                 : true,
            "mesh_prediction"                : true,
            "solve_mesh_at_each_iteration"   : true,
            "coupling_strategy" : {
                "solver_type"       : "Relaxation",
                "acceleration_type" : "Aitken",
                "w_0"               : 0.825
                },
            "mesh_solver"                    : "mesh_solver_structural_similarity",
            "mesh_reform_dofs_each_step"     : false,
            "structure_interfaces_list"      : [""],
            "fluid_interfaces_list"          : [""]
            }
        }
        """)

        # Time stepping checks (no sub-stepping between subdomains has been implemented yed)
        time_step_structure = project_parameters["structure_solver_settings"]["problem_data"]["time_step"].GetDouble()
        # If automatic time stepping has been selected in the fluid domain, deactivate it and use the structure time step
        if (project_parameters["fluid_solver_settings"]["solver_settings"]["time_stepping"]["automatic_time_step"].GetBool()):
            project_parameters["fluid_solver_settings"]["solver_settings"]["time_stepping"]["automatic_time_step"].SetBool(False)
            time_step_fluid = time_step_structure
            print("WARNING: Automatic fluid time stepping cannot be used. Setting structure time step as fluid time step.")
        else:
            time_step_fluid = project_parameters["fluid_solver_settings"]["solver_settings"]["time_stepping"]["time_step"].GetDouble()
            if time_step_structure != time_step_fluid:
                raise("ERROR: Different time step among subdomains! No sub-stepping implemented yet.")

        self.time_step = time_step_fluid

        # Take the each one of the solvers settings from the ProjectParameters
        self.settings = KratosMultiphysics.Parameters("{}")
        self.settings.AddValue("structure_solver_settings",project_parameters["structure_solver_settings"]["solver_settings"])
        self.settings.AddValue("fluid_solver_settings",project_parameters["fluid_solver_settings"]["solver_settings"])
        self.settings.AddValue("coupling_solver_settings",project_parameters["coupling_solver_settings"]["solver_settings"])

        # Overwrite the default settings with user-provided parameters
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)

        # Auxiliar variables
        self.max_nl_it = self.settings["coupling_solver_settings"]["nl_max_it"].GetInt()
        self.nl_tol = self.settings["coupling_solver_settings"]["nl_tol"].GetDouble()
        self.solve_mesh_at_each_iteration = self.settings["coupling_solver_settings"]["solve_mesh_at_each_iteration"].GetBool()
        # self.move_interface = self.settings["coupling_solver_settings"]["move_interface"].GetBool()
        # self.mesh_prediction = self.settings["coupling_solver_settings"]["mesh_prediction"].GetBool()
        self.coupling_algorithm = self.settings["coupling_solver_settings"]["coupling_scheme"].GetString()
        self.fluid_interface_submodelpart_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][0].GetString()
        self.structure_interface_submodelpart_name = self.settings["coupling_solver_settings"]["structure_interfaces_list"][0].GetString()
        coupling_utility_parameters = self.settings["coupling_solver_settings"]["coupling_strategy"]


        print("*** Partitioned FSI solver construction starts...")

        # Construct the structure solver
        structure_solver_module = __import__(self.settings["structure_solver_settings"]["solver_type"].GetString())
        self.structure_solver = structure_solver_module.CreateSolver(self.structure_main_model_part,
                                                                     self.settings["structure_solver_settings"])
        print("* Structure solver constructed.")

        # Construct the fluid solver
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolver(self.fluid_main_model_part,
                                                                      project_parameters["fluid_solver_settings"])
        print("* Fluid solver constructed.")

        # Construct the coupling partitioned strategy
        import convergence_accelerator_factory
        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(coupling_utility_parameters)
        print("* Coupling strategy constructed.")

        # Construct the ALE mesh solver
        mesh_solver_settings = KratosMultiphysics.Parameters("{}")
        mesh_solver_settings.AddValue("mesh_reform_dofs_each_step",self.settings["coupling_solver_settings"]["mesh_reform_dofs_each_step"])

        self.mesh_solver_module = __import__(self.settings["coupling_solver_settings"]["mesh_solver"].GetString())
        self.mesh_solver = self.mesh_solver_module.CreateSolver(self.fluid_solver.main_model_part,
                                                                mesh_solver_settings)
        print("* ALE mesh solver constructed.")

        print("*** FSI partitioned solver construction finished.")


    def GetMinimumBufferSize(self):
        # Get structure buffer size
        buffer_structure = self.structure_solver.GetMinimumBufferSize()
        # Get fluid buffer size
        buffer_fluid = self.fluid_solver.GetMinimumBufferSize()

        return min(buffer_structure,buffer_fluid)


    def AddVariables(self):
        ## Structure variables addition
        # Standard CSM variables addition
        self.structure_solver.AddVariables()

        ## Fluid variables addition
        # Standard CFD variables addition
        self.fluid_solver.AddVariables()
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        # Mesh solver variables addition
        self.mesh_solver.AddVariables()

        ## FSIApplication variables addition
        NonConformant_OneSideMap.AddVariables(self.fluid_solver.main_model_part,self.structure_solver.main_model_part)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.VECTOR_PROJECTED)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.FSI_INTERFACE_RESIDUAL)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.FSI_INTERFACE_MESH_RESIDUAL)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.VECTOR_PROJECTED)


    def ImportModelPart(self):
        # Import structure model part
        self.structure_solver.ImportModelPart()

        # Import fluid model part
        self.fluid_solver.ImportModelPart()


    def AddDofs(self):
        # Add DOFs structure
        self.structure_solver.AddDofs()

        # Add DOFs fluid
        self.fluid_solver.AddDofs()
        self.mesh_solver.AddDofs()


    def Initialize(self):

        # Get the domain size
        self.domain_size = self._GetDomainSize()

        # Get the partitioned FSI utilities
        self.partitioned_fsi_utilities = self._GetPartitionedFSIUtilities()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        if self.coupling_algorithm == "DirichletNeumann":
            # Ensure that the fluid reaction fluxes are computed if D-N scheme is considered
            if self.fluid_solver.settings["compute_reactions"].GetBool() == False:
                self.fluid_solver.settings["compute_reactions"].SetBool(True)
            # In the D-N scheme the interface correction is done over the velocity
            self.correction_over_velocity = True

        elif self.coupling_algorithm == "NeumannNeumann":
            # In the N-N scheme the interface correction is done over the interface fluxes
            self.correction_over_velocity = False

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Python mesh solver initialization
        self.mesh_solver.Initialize()

        # Construct the interface mapper
        # Recall, to set the INTERFACE flag in both the fluid and solid interface before the mapper construction
        # Currently this is done with the FSI application Python process set_interface_process.py
        search_radius_factor = 2.0
        mapper_max_iterations = 100
        mapper_tolerance = 0.001*self.nl_tol
        self.interface_mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.fluid_solver.main_model_part,
                                                                                  self.structure_solver.main_model_part,
                                                                                  search_radius_factor,
                                                                                  mapper_max_iterations,
                                                                                  mapper_tolerance)

        # Set the Neumann B.C. in the structure interface
        self._SetStructureNeumannCondition() #TODO change when the interface is able to correctly transfer distributed forces

        # Set the Neumann B.C. in the fluid interface
        if self.coupling_algorithm == "NeumannNeumann":
            self._SetFluidNeumannCondition()

        # Note that the FSI problem is defined in terms of the fluid interface
        # Initialize the iteration value for the residual computation
        fluid_interface_residual_size = self.partitioned_fsi_utilities.GetFluidInterfaceResidualSize()
        self.iteration_value = KratosMultiphysics.Vector(fluid_interface_residual_size)     # Interface solution guess (it might be velocity or fluxes depending on the type of coupling)
        for i in range(0,fluid_interface_residual_size):
            self.iteration_value[i] = 0.0001

        # Compute the fluid domain NODAL_AREA values (required as weight in the residual norm computation)
        KratosMultiphysics.CalculateNodalAreaProcess(self.fluid_solver.GetComputingModelPart(), self.domain_size).Execute()

        # Strategies initialization
        self.fluid_solver.SolverInitialize()
        self.structure_solver.SolverInitialize()
        self.coupling_utility.Initialize()


    def GetComputingModelPart(self):
        pass


    def GetOutputVariables(self):
        pass


    def ComputeDeltaTime(self):
        return self.time_step


    def SaveRestart(self):
        pass


    def Solve(self):

        self.fluid_solver.SolverInitializeSolutionStep()
        self.structure_solver.SolverInitializeSolutionStep()
        self.coupling_utility.InitializeSolutionStep()

        self.fluid_solver.SolverPredict()
        self.structure_solver.SolverPredict()

        ## Compute mesh prediction ##
        self._ComputeMeshPrediction()

        ## Non-Linear interface coupling iteration ##
        for nl_it in range(1,self.max_nl_it+1):

            print("     NL-ITERATION ",nl_it,"STARTS.")
            self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            self.coupling_utility.InitializeNonLinearIteration()

            # Update the fluid interface according to the previous iteration
            # If the correction is done over the velocity, the interface displacement must be done according the corrected velocity.
            if self.correction_over_velocity == True:
                self._ComputeCorrectedInterfacePosition()

            # The interface movement can be directly done with the obtained values
            elif self.correction_over_velocity == False:
                    # for mapper in self.list_of_mappers:
                    #     mapper.InverseMap(KratosALE.MESH_DISPLACEMENT, KratosMultiphysics.DISPLACEMENT)
                    keep_sign = True
                    distribute_load = False
                    self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                                     KratosALE.MESH_DISPLACEMENT,
                                                                     keep_sign,
                                                                     distribute_load)

            if (self.solve_mesh_at_each_iteration == True):
                self.mesh_solver.Solve()
            else:
                self.mesh_solver.MoveNodes()

            # Residual computation
            print("     Residual computation starts...")

            if self.coupling_algorithm == "DirichletNeumann":
                vel_residual = self._ComputeDirichletNeumannResidual()

            elif self.coupling_algorithm == "NeumannNeumann":
                vel_residual = self._ComputeNeumannNeumannResidual()

            nl_res_norm = self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.FSI_INTERFACE_RESIDUAL_NORM]
            FluidInterfaceArea = self.partitioned_fsi_utilities.GetFluidInterfaceArea()
            print("     Residual computation finished. |res|/A =", nl_res_norm/FluidInterfaceArea)

            # Check convergence
            if nl_res_norm/FluidInterfaceArea < self.nl_tol:
                print("     NON-LINEAR ITERATION CONVERGENCE ACHIEVED")
                print("     Total non-linear iterations: ",nl_it," |res|/A = ",nl_res_norm/FluidInterfaceArea)
                break

            else:
                # If convergence is not achieved, perform the correction of the prediction
                print("     Performing non-linear iteration ",nl_it," correction.")
                self.coupling_utility.UpdateSolution(vel_residual, self.iteration_value)
                self.coupling_utility.FinalizeNonLinearIteration()

        ## Compute the mesh residual
        self.partitioned_fsi_utilities.ComputeFluidInterfaceMeshVelocityResidualNorm()
        mesh_res_norm = self.fluid_solver.main_model_part.ProcessInfo.GetValue(KratosFSI.FSI_INTERFACE_MESH_RESIDUAL_NORM)
        print("     NL residual norm: ", nl_res_norm)
        print("     Mesh residual norm: ", mesh_res_norm)

        ## Finalize solution step
        self.fluid_solver.SolverFinalizeSolutionStep()
        self.structure_solver.SolverFinalizeSolutionStep()
        self.coupling_utility.FinalizeSolutionStep()


    def SetEchoLevel(self, structure_echo_level, fluid_echo_level):
        self.structure_solver.SetEchoLevel(self, structure_echo_level)
        self.fluid_solver.SetEchoLevel(self, fluid_echo_level)


    def SetTimeStep(self, step):
        self.fluid_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME_STEPS, step)
        self.structure_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME_STEPS, step)


    def Clear(self):
        self.fluid_solver.Clear()
        self.structure_solver.Clear()


    def Check(self):
        self.fluid_solver.Check()
        self.structure_solver.Check()


    #######################################################################
    ##############          PRIVATE METHODS SECTION          ##############
    #######################################################################

    ### AUXILIAR METHODS ###
    def _GetFluidInterfaceSubmodelPart(self):
        return self.fluid_solver.main_model_part.GetSubModelPart(self.fluid_interface_submodelpart_name)


    def _GetStructureInterfaceSubmodelPart(self):
        return self.structure_solver.main_model_part.GetSubModelPart(self.structure_interface_submodelpart_name)


    def _GetDomainSize(self):

        fluid_domain_size = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        structure_domain_size = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        if fluid_domain_size !=structure_domain_size:
            raise("ERROR: Solid domain size and fluid domain size are not equal!")

        return fluid_domain_size


    def _GetPartitionedFSIUtilities(self):

        if (self.domain_size == 2):
            return KratosFSI.PartitionedFSIUtilities2D(self._GetFluidInterfaceSubmodelPart(), self._GetStructureInterfaceSubmodelPart())
        else:
            return KratosFSI.PartitionedFSIUtilities3D(self._GetFluidInterfaceSubmodelPart(), self._GetStructureInterfaceSubmodelPart())


    def _SetStructureNeumannCondition(self):

        structure_computational_submodelpart = self.structure_solver.GetComputingModelPart()

        aux_count = 0
        for cond in self.structure_solver.main_model_part.Conditions:
            if(cond.Id > aux_count):
                aux_count = cond.Id

        for i in range(self.settings["coupling_solver_settings"]["structure_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["structure_interfaces_list"][i].GetString()
            interface_submodelpart_i = self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name)
            # NOTE: In this manner, two interface submodelparts cannot share a node (it would be repeated in the pointload conditions...)

            # Create the point load condition
            if self.domain_size == 2:
                for node in interface_submodelpart_i.Nodes:
                    aux_count+=1
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition2D1N",aux_count,[node.Id],self.structure_solver.main_model_part.Properties[0])

            elif self.domain_size == 3:
                for node in interface_submodelpart_i.Nodes:
                    aux_count+=1
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition3D1N",aux_count,[node.Id],self.structure_solver.main_model_part.Properties[0])


    # TODO: This function must be checked as soon as the fluid Neumann BC has been implemented.
    def _SetFluidNeumannCondition(self):

        fluid_computational_volume_submodelpart = self.fluid_solver.GetComputingModelPart()

        aux_count = len(self.fluid_solver.main_model_part.Conditions)       # Get the last existing condition numbering
        aux_count += 1
        print("max aux_count",aux_count)
        aux_count = 0
        for cond in self.fluid_solver.main_model_part.Conditions:
            if(cond.Id > aux_count):
                aux_count = cond.Id
        aux_count += 1
        print("max aux_count",aux_count)


        for i in range(self.settings["coupling_solver_settings"]["fluid_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][i].GetString()
            interface_submodelpart_i = self.fluid_solver.main_model_part.GetSubModelPart(interface_submodelpart_name)
            # NOTE: In this manner, two interface submodelparts cannot share a node (it would be repeated in the pointload conditions...)
            # DO CreateNewCondition CHECK IF THERE EXIST A CONDITION IN A NODE?

            for node in interface_submodelpart_i.Nodes:

                # NOTE: THIS CONDITION REMAINS TO BE IMPLEMENTED IN THE FluidDynamicsApplication, DECIDE WHAT TO DO.
                # Create the fluid load condition
                if self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                    fluid_computational_volume_submodelpart.CreateNewCondition("PointForce2Dfluid",aux_count,[node.Id],self.fluid_solver.main_model_part.Properties[0])
                elif self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
                    fluid_computational_volume_submodelpart.CreateNewCondition("PointForce3Dfluid",aux_count,[node.Id],self.fluid_solver.main_model_part.Properties[0])

                aux_count+=1


    def _ComputeMeshPrediction(self):

            print("Computing time step ",self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS]," prediction...")
            # Get the previous step fluid interface nodal fluxes
            keep_sign = False
            distribute_load = True
            self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                             KratosSolid.POINT_LOAD,
                                                             keep_sign,
                                                             distribute_load)

            # Solve the current step structure problem with the previous step fluid interface nodal fluxes
            self.structure_solver.SolverSolveSolutionStep()

            # Map the obtained structure displacement to the fluid interface
            keep_sign = True
            distribute_load = False
            self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                             KratosALE.MESH_DISPLACEMENT,
                                                             keep_sign,
                                                             distribute_load)

            # Solve the mesh problem
            self.mesh_solver.Solve()

            print("Mesh prediction computed.")

    ### RESIDUALS ###

    # The residual schemes have been implemented such that they retur a "vel_residual" vector.
    # "vel_residual" is a vector containing the nodal values of the residual at the fluid interface.
    # The residual computation starts from the iteration prediction stored in "self.iteration_value"

    # Dirichlet-Neumann scheme interface velocity residual
    def _ComputeDirichletNeumannResidual(self):

        # Fluid domain velocities imposition
        self.partitioned_fsi_utilities.SetAndFixFluidInterfaceVectorVariable(KratosMultiphysics.VELOCITY, True, self.iteration_value)
        self.partitioned_fsi_utilities.SetAndFixFluidInterfaceVectorVariable(KratosMultiphysics.MESH_VELOCITY, False, self.iteration_value)

        # Solve fluid problem
        self.fluid_solver.SolverSolveSolutionStep()

        # Transfer fluid reaction to solid interface
        keep_sign = False
        distribute_load = True
        self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                         KratosSolid.POINT_LOAD,
                                                         keep_sign,
                                                         distribute_load)

        # Solve structure problem
        self.structure_solver.SolverSolveSolutionStep()

        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.VELOCITY,
                                                         KratosFSI.VECTOR_PROJECTED,
                                                         keep_sign,
                                                         distribute_load)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        vel_residual = KratosMultiphysics.Vector(self.partitioned_fsi_utilities.GetFluidInterfaceResidualSize())
        self.partitioned_fsi_utilities.ComputeFluidInterfaceVelocityResidual(vel_residual)

        return vel_residual


    # Neumann-Neumann scheme interface velocity residual
    def _ComputeNeumannNeumannResidual(self):

        # Fluid domain interface flux imposition
        self.partitioned_fsi_utilities.SetAndFixFluidInterfaceVectorVariable(KratosMultiphysics.FORCE, False, self.iteration_value)

        # Solve the fluid problem
        self.fluid_solver.SolverSolveSolutionStep()

        # Flux prediction is defined in terms of the fluid interface. Transfer it to the structure interface.
        keep_sign = False
        distribute_load = True
        self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.FORCE,
                                                         KratosSolid.POINT_LOAD,
                                                         keep_sign,
                                                         distribute_load)

        # Solve structure problem
        self.structure_solver.SolverSolveSolutionStep()

        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.VELOCITY,
                                                         KratosFSI.VECTOR_PROJECTED,
                                                         keep_sign,
                                                         distribute_load)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        vel_residual = KratosMultiphysics.Vector(self.partitioned_fsi_utilities.GetFluidInterfaceResidualSize())
        self.partitioned_fsi_utilities.ComputeFluidInterfaceVelocityResidual(vel_residual)

        return vel_residual


    ### INTERFACE MOVEMENT UTILITY ###

    # Function to update the position of the interface during iterations
    def _ComputeCorrectedInterfacePosition(self):

        # Bossak parameters
        alpha = -1/3
        gamma = 0.5*(1-2*alpha)
        beta = ((1-alpha)**2)/4

        i = 0
        if self.domain_size == 2:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                u_n = node.GetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,1)
                v_n = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
                a_n = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION,1)

                v_n1 = KratosMultiphysics.Vector(3)
                v_n1[0] = self.iteration_value[i]
                v_n1[1] = self.iteration_value[i+1]
                v_n1[2] = 0.0
                i+=2

                # Compute the current acceleration associated to the corrected interface velocity with the Bossak formulaes
                a_n1 = KratosMultiphysics.Vector(3)
                a_n1[0] = (v_n1[0] - v_n[0] - self.time_step*(1-gamma*(alpha-1))*a_n[0]) / ((1-alpha)*self.time_step*gamma)
                a_n1[1] = (v_n1[1] - v_n[1] - self.time_step*(1-gamma*(alpha-1))*a_n[1]) / ((1-alpha)*self.time_step*gamma)
                a_n1[2] = 0.0

                # Compute the current displacement associated to the corrected interface velocity with the Bossak formulaes
                u_n1 = KratosMultiphysics.Vector(3)
                u_n1[0] = u_n[0] + self.time_step*v_n[0] + (self.time_step**2)*(0.5-beta)*a_n[0] + (self.time_step**2)*beta*((1-alpha)*a_n1[0] + alpha*a_n[0])
                u_n1[1] = u_n[1] + self.time_step*v_n[1] + (self.time_step**2)*(0.5-beta)*a_n[1] + (self.time_step**2)*beta*((1-alpha)*a_n1[1] + alpha*a_n[1])
                u_n1[2] = 0.0

                # Set the obtained corrected interface displacement
                node.SetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,0,u_n1)

        elif self.domain_size == 3:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                u_n = node.GetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,1)
                v_n = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
                a_n = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION,1)

                v_n1 = KratosMultiphysics.Vector(3)
                v_n1[0] = self.iteration_value[i]
                v_n1[1] = self.iteration_value[i+1]
                v_n1[2] = self.iteration_value[i+2]
                i+=3

                # Compute the current acceleration associated to the corrected interface velocity with the Bossak formulaes
                a_n1 = KratosMultiphysics.Vector(3)
                a_n1[0] = (v_n1[0] - v_n[0] - self.time_step*(1-gamma*(alpha-1))*a_n[0]) / ((1-alpha)*self.time_step*gamma)
                a_n1[1] = (v_n1[1] - v_n[1] - self.time_step*(1-gamma*(alpha-1))*a_n[1]) / ((1-alpha)*self.time_step*gamma)
                a_n1[2] = (v_n1[2] - v_n[2] - self.time_step*(1-gamma*(alpha-1))*a_n[2]) / ((1-alpha)*self.time_step*gamma)

                # Compute the current displacement associated to the corrected interface velocity with the Bossak formulaes
                u_n1 = KratosMultiphysics.Vector(3)
                u_n1[0] = u_n[0] + self.time_step*v_n[0] + (self.time_step**2)*(0.5-beta)*a_n[0] + (self.time_step**2)*beta*((1-alpha)*a_n1[0] + alpha*a_n[0])
                u_n1[1] = u_n[1] + self.time_step*v_n[1] + (self.time_step**2)*(0.5-beta)*a_n[1] + (self.time_step**2)*beta*((1-alpha)*a_n1[1] + alpha*a_n[1])
                u_n1[2] = u_n[2] + self.time_step*v_n[2] + (self.time_step**2)*(0.5-beta)*a_n[2] + (self.time_step**2)*beta*((1-alpha)*a_n1[2] + alpha*a_n[2])

                # Set the obtained corrected interface displacement
                node.SetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,0,u_n1)
