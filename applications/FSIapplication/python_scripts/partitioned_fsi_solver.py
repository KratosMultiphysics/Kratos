from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
# import NonConformant_OneSideMap                # Import non-conformant mapper

# Import libraries
#~ import time as timemodule                   # Import time library as timemodule (avoid interferences with "time" var)
#~ import json                                 # Encoding library (for data exchange)

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as KratosALE
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Define the UblasSparseSpace to compute auxilar array operations.
space = KratosMultiphysics.UblasSparseSpace()

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(structure_main_model_part, fluid_main_model_part, project_parameters):
    return PartitionedFSISolver(structure_main_model_part, fluid_main_model_part, project_parameters)


class PartitionedFSISolver:
    def __init__(self, structure_main_model_part, fluid_main_model_part, project_parameters):

        # Initial tests
        self.time_step_structure = project_parameters["structure_solver_settings"]["problem_data"]["time_step"].GetDouble()
        self.time_step_fluid = project_parameters["fluid_solver_settings"]["problem_data"]["time_step"].GetDouble()
        start_time_structure = project_parameters["structure_solver_settings"]["problem_data"]["start_time"].GetDouble()
        start_time_fluid = project_parameters["fluid_solver_settings"]["problem_data"]["start_step"].GetDouble()
        end_time_structure = project_parameters["structure_solver_settings"]["problem_data"]["end_time"].GetDouble()
        end_time_fluid = project_parameters["fluid_solver_settings"]["problem_data"]["end_time"].GetDouble()

        if self.time_step_structure != self.time_step_fluid:
            raise("ERROR: Different time step among subdomains! No sub-stepping implemented yet.")
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
            "displacement_relative_tolerance": 1.0e-3,
            "displacement_absolute_tolerance": 1.0e-5,
            "residual_relative_tolerance": 1.0e-3,
            "residual_absolute_tolerance": 1.0e-5,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 200,
                "tolerance": 1e-7,
                "scaling": false,
                "verbosity": 1
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
            "dynamic_tau": 0.0,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": true,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_step": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                "solver_type" : "AMGCL_NS_Solver",
                "krylov_type" : "lgmres",
                "velocity_block_preconditioner" : {
                    "tolerance" : 1e-3,
                    "precondioner_type" : "spai0"
                },
                "pressure_block_preconditioner" : {
                    "tolerance" : 1e-2,
                    "precondioner_type" : "spai0"
                },
                "tolerance" : 1e-7,
                "gmres_krylov_space_dimension": 50,
                "max_iteration": 50,
                "verbosity" : 0,
                "scaling": true,
                "coarse_enough" : 5000
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts": [""],
            "alpha":-0.3,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "regularization_coef": 1000,
            "MoveMeshFlag": false,
            "use_slip_conditions": false,
            "turbulence_model": "None",
            "use_spalart_allmaras": false
            },
        "coupling_solver_settings":
            {
            "coupling_scheme"   : "DirichletNeumann",
            "solver_type"       : "partitioned_fsi_solver",
            "nl_tol"            : 1e-5,
            "nl_max_it"         : 50,
            "move_interface"    : true,
            "mesh_prediction"   : true,
            "coupling_strategy" : {
                "solver_type"       : "Relaxation",
                "acceleration_type" : "Aitken",
                "w_0"               : 0.825
                },
            "mesh_solver"               : "mesh_solver_structural_similarity",
            "mesh_reform_dofs_each_step": false,
            "structure_interfaces_list" : [""],
            "fluid_interfaces_list" : [""],
        	"mappers_settings": [
                    {
        		        "mapper_type": "NearestNeighbor",
        		        "interface_submodel_part_origin": "FluidInterface",
        		        "interface_submodel_part_destination": "StructureInterface"
        	        }
                ]
            }
        }
        """)

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
        self.move_interface = self.settings["coupling_solver_settings"]["move_interface"].GetBool()
        self.mesh_prediction = self.settings["coupling_solver_settings"]["mesh_prediction"].GetBool()
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
        fluid_solver_module = __import__(self.settings["fluid_solver_settings"]["solver_type"].GetString())
        self.fluid_solver = fluid_solver_module.CreateSolver(self.fluid_main_model_part,
                                                             self.settings["fluid_solver_settings"])
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
        # NonConformant_OneSideMap.AddVariables(self.fluid_solver.main_model_part,self.structure_solver.main_model_part)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosFSI.VECTOR_PROJECTED)
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
        # Initialize structure solver
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

        # Initialize fluid solver
        self.fluid_solver.Initialize()

        # Mesh solver initialization
        self.mesh_solver.Initialize()

        # Get the domain size
        self.domain_size = self._GetDomainSize()

        # Construct the interface mapper
        # Recall, to set the INTERFACE flag in both the fluid and solid interface before the mapper construction
        # Currently this is done with the FSI application Python process set_interface_process.py
        # search_radius_factor = 2.0
        # mapper_max_iterations = 50
        # mapper_tolerance = 1e-5
        # self.interface_mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.fluid_solver.main_model_part,
        #                                                                           self.structure_solver.main_model_part,
        #                                                                           search_radius_factor,
        #                                                                           mapper_max_iterations,
        #                                                                           mapper_tolerance)

        mappers_settings = self.settings["coupling_solver_settings"]["mappers_settings"]

        self.list_of_mappers = []

        for mapper_id in range(0,mappers_settings.size()):
            self.list_of_mappers.append(KratosMapping.MapperFactory(self.fluid_solver.main_model_part,      # Origin model part
                                                                    self.structure_solver.main_model_part,  # Destination model part
                                                                    mappers_settings[mapper_id]))

            print("Mapper ",mapper_id," from MappingApplication constructed.")

        # Set the Neumann B.C. in the structure interface
        self._SetStructureNeumannCondition() #TODO change when the interface is able to correctly transfer distributed forces

        # Set the Neumann B.C. in the fluid interface
        if self.coupling_algorithm == "NeumannNeumann":
            self._SetFluidNeumannCondition()

        # Get interface problem sizes
        interface_problem_sizes = self._GetInterfaceProblemSizes()
        fluid_interface_problem_size = interface_problem_sizes[1]

        # Note that the FSI problem is defined in terms of the fluid interface
        self.fluid_interface_residual_size = fluid_interface_problem_size*self.domain_size
        self.iteration_value = KratosMultiphysics.Vector(self.fluid_interface_residual_size)     # Interface solution guess (it might be velocity or fluxes depending on the type of coupling)

        for i in range(0,self.fluid_interface_residual_size):
            self.iteration_value[i] = 0.0001

        self.fluid_solver.SolverInitialize()
        self.structure_solver.SolverInitialize()
        self.coupling_utility.Initialize()


    def GetComputingModelPart(self):
        pass


    def GetOutputVariables(self):
        pass


    def ComputeDeltaTime(self):
        pass


    def SaveRestart(self):
        pass


    def Solve(self):

        self.fluid_solver.SolverInitializeSolutionStep()
        self.structure_solver.SolverInitializeSolutionStep()
        self.coupling_utility.InitializeSolutionStep()

        self.fluid_solver.SolverPredict()
        self.structure_solver.SolverPredict()

        ## Compute mesh prediction ##
        if self.mesh_prediction == True:
            self._ComputeMeshPrediction()

        ## Non-Linear interface coupling iteration ##
        for nl_it in range(1,self.max_nl_it+1):

            self.coupling_utility.InitializeNonLinearIteration()

            print("     NL-ITERATION ",nl_it,"STARTS.")
            self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            # Residual computation
            print("     Residual computation starts...")

            if self.coupling_algorithm == "DirichletNeumann":
                vel_residual = self._ComputeDirichletNeumannResidual()

            elif self.coupling_algorithm == "NeumannNeumann":
                vel_residual = self._ComputeNeumannNeumannResidual()

            nl_res_norm = self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.FSI_INTERFACE_RESIDUAL_NORM]
            print("     Residual computation finished. |res|=",nl_res_norm)

            # Check convergence
            if nl_res_norm < self.nl_tol:
                print("     NON-LINEAR ITERATION CONVERGENCE ACHIEVED")
                print("     Total non-linear iterations: ",nl_it," NL residual norm: ",nl_res_norm)
                break

            else:
                # If convergence is not achieved, perform the correction of the prediction
                print("     Performing non-linear iteration ",nl_it," correction.")

                self.coupling_utility.UpdateSolution(vel_residual,
                                                     self.iteration_value)

                # Move interface nodes
                # If the correction is done over the velocity, the interface displacement must be done according the corrected velocity.
                if self.correction_over_velocity == True:
                    self._ComputeCorrectedInterfacePosition()

                # The interface movement can be directly done with the obtained values
                elif self.correction_over_velocity == False:
                    if self.move_interface == True:
                        for mapper in self.list_of_mappers:
                            mapper.InverseMap(KratosALE.MESH_DISPLACEMENT, KratosMultiphysics.DISPLACEMENT)
                        # keep_sign = True
                        # distribute_load = False
                        # self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                        #                                                  KratosALE.MESH_DISPLACEMENT,
                        #                                                  keep_sign,
                        #                                                  distribute_load)            # Project the structure interface displacement onto the fluid interface

                self.mesh_solver.MoveNodes()

                self.coupling_utility.FinalizeNonLinearIteration()

        ## Mesh update
        for mapper in self.list_of_mappers:
            mapper.InverseMap(KratosALE.MESH_DISPLACEMENT, KratosMultiphysics.DISPLACEMENT)
        # keep_sign = True
        # distribute_load = False
        # self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
        #                                                  KratosALE.MESH_DISPLACEMENT,
        #                                                  keep_sign,
        #                                                  distribute_load)                       # Project the structure interface displacement onto the fluid interface
        self.mesh_solver.Solve()                                                                # Solve the mesh problem

        ## Compute the mesh residual
        mesh_res_norm = self._ComputeMeshResidual()
        print("     Mesh residual norm: ",mesh_res_norm)

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


    def _GetInterfaceProblemSizes(self):

        # Get solid interface problem size
        structure_interface_pb_size = 0
        for i in range(self.settings["coupling_solver_settings"]["structure_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["structure_interfaces_list"][i].GetString()
            structure_interface_pb_size += len(self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name).Nodes)

        # Get fluid interface problem size
        fluid_interface_pb_size = 0
        for i in range(self.settings["coupling_solver_settings"]["fluid_interfaces_list"].size()):
            interface_submodelpart_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][i].GetString()
            fluid_interface_pb_size += len(self.fluid_solver.main_model_part.GetSubModelPart(interface_submodelpart_name).Nodes)

        return (structure_interface_pb_size, fluid_interface_pb_size)


    def _GetDomainSize(self):

        fluid_domain_size = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        structure_domain_size = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        if fluid_domain_size !=structure_domain_size:
            raise("ERROR: Solid domain size and fluid domain size are not equal!")

        return fluid_domain_size


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
            # keep_sign = False
            # distribute_load = True
            # self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.REACTION,
            #                                                  KratosSolid.POINT_LOAD,
            #                                                  keep_sign,
            #                                                  distribute_load)
            for mapper in self.list_of_mappers:
                mapper.Map(KratosMultiphysics.REACTION, KratosSolid.POINT_LOAD, KratosMapping.MapperFactory.SWAP_SIGN | KratosMapping.MapperFactory.CONSERVATIVE)

            # Solve the current step structure problem with the previous step fluid interface nodal fluxes
            self.structure_solver.SolverSolveSolutionStep()

            # Map the obtained structure displacement to the fluid interface
            # keep_sign = True
            # distribute_load = False
            # self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
            #                                                  KratosALE.MESH_DISPLACEMENT,
            #                                                  keep_sign,
            #                                                  distribute_load)
            for mapper in self.list_of_mappers:
                mapper.InverseMap(KratosALE.MESH_DISPLACEMENT, KratosMultiphysics.DISPLACEMENT)

            # Solve the mesh problem
            self.mesh_solver.Solve()

            # Map the obtained structure velocity to the fluid interface
            #~ self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.VELOCITY,
                                                             #~ KratosMultiphysics.VELOCITY,
                                                             #~ True,False)

            # Solve the fluid problem with the updated mesh and the solid interface velocity as Dirichlet B.C.
            #~ self.fluid_solver.SolveSolutionStep()
            print("Mesh prediction computed.")


    ### RESIDUALS ###

    # The residual schemes have been implemented such that they have to modify "self.vel_residual"
    # "self.vel_residual" is a vector containing the nodal residual at the fluid interface
    # The residual computation starts from the iteration prediction stored in "self.iteration_value"

    # Dirichlet-Neumann scheme interface velocity residual
    def _ComputeDirichletNeumannResidual(self):

        # Fluid domain velocity imposition
        i = 0
        if self.domain_size == 2:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                velocity = KratosMultiphysics.Vector(3)
                velocity[0] = self.iteration_value[i]
                velocity[1] = self.iteration_value[i+1]
                velocity[2] = 0.0

                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)

                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0,velocity)
                i+=2

        elif self.domain_size == 3:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                velocity = KratosMultiphysics.Vector(3)
                velocity[0] = self.iteration_value[i]
                velocity[1] = self.iteration_value[i+1]
                velocity[2] = self.iteration_value[i+2]

                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)

                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0,velocity)
                i+=3

        # Solve fluid problem
        self.fluid_solver.SolverSolveSolutionStep()

        # Transfer fluid reaction to solid interface
        # keep_sign = False
        # distribute_load = True
        # self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.REACTION,
        #                                                  KratosSolid.POINT_LOAD,
        #                                                  keep_sign,
        #                                                  distribute_load)
        for mapper in self.list_of_mappers:
            mapper.Map(KratosMultiphysics.REACTION, KratosSolid.POINT_LOAD, KratosMapping.MapperFactory.SWAP_SIGN | KratosMapping.MapperFactory.CONSERVATIVE)

        # Solve structure problem
        self.structure_solver.SolverSolveSolutionStep()

        # Project the structure velocity onto the fluid interface
        # keep_sign = True
        # distribute_load = False
        # self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.VELOCITY,
        #                                                  KratosFSI.VECTOR_PROJECTED,
        #                                                  keep_sign,
        #                                                  distribute_load)
        for mapper in self.list_of_mappers:
            mapper.InverseMap(KratosFSI.VECTOR_PROJECTED, KratosMultiphysics.VELOCITY)

        # Compute the fluid interface residual by means of the VECTOR_PROJECTED variable
        vel_residual = KratosMultiphysics.Vector(self.fluid_interface_residual_size)

        i = 0
        if self.domain_size == 2:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                vector_projected = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED,0)

                vel_residual[i] = vector_projected[0] - self.iteration_value[i]
                vel_residual[i+1] = vector_projected[1] - self.iteration_value[i+1]
                i+=2

        elif self.domain_size == 3:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                vector_projected = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED,0)

                vel_residual[i] = vector_projected[0] - self.iteration_value[i]
                vel_residual[i+1] = vector_projected[1] - self.iteration_value[i+1]
                vel_residual[i+2] = vector_projected[2] - self.iteration_value[i+2]
                i+=3

        # Compute the residual norm and store it in the ProcessInfo
        res_norm = space.TwoNorm(vel_residual)
        self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.FSI_INTERFACE_RESIDUAL_NORM] = res_norm

        return vel_residual


    # Neumann-Neumann scheme interface velocity residual
    def _ComputeNeumannNeumannResidual(self):

        # Fluid domain interface flux imposition
        i = 0
        if self.domain_size == 2:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                force = KratosMultiphysics.Vector(3)
                force[0] = self.iteration_value[i]
                force[1] = self.iteration_value[i+1]
                force[2] = 0.0

                node.SetSolutionStepValue(KratosMultiphysics.FORCE,0,force)
                i+=2

        elif self.domain_size == 3:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                force = KratosMultiphysics.Vector(3)
                force[0] = self.iteration_value[i]
                force[1] = self.iteration_value[i+1]
                force[2] = self.iteration_value[i+2]

                node.SetSolutionStepValue(KratosMultiphysics.FORCE,0,force)
                i+=3

        # Solve the fluid problem
        self.fluid_solver.SolverSolveSolutionStep()

        # Flux prediction is defined in terms of the fluid interface. Transfer it to the structure interface.
        # keep_sign = False
        # distribute_load = True
        # self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.FORCE,
        #                                                  KratosSolid.POINT_LOAD,
        #                                                  keep_sign,
        #                                                  distribute_load)
        for mapper in self.list_of_mappers:
            mapper.Map(KratosMultiphysics.FORCE, KratosSolid.POINT_LOAD, KratosMapping.MapperFactory.SWAP_SIGN | KratosMapping.MapperFactory.CONSERVATIVE)

        # Solve structure problem
        self.structure_solver.SolverSolveSolutionStep()

        # Project the structure velocity onto the fluid interface
        # keep_sign = True
        # distribute_load = False
        # self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.VELOCITY,
        #                                                  KratosFSI.VECTOR_PROJECTED,
        #                                                  keep_sign,
        #                                                  distribute_load)
        for mapper in self.list_of_mappers:
            mapper.InverseMap(KratosFSI.VECTOR_PROJECTED, KratosMultiphysics.VELOCITY)

        # Compute the fluid interface residual by means of the VECTOR_PROJECTED variable
        vel_residual = KratosMultiphysics.Vector(self.fluid_interface_residual_size)

        i = 0
        if self.domain_size == 2:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                velocity_fluid = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                vector_projected = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED,0)

                vel_residual[i] = velocity_fluid[0] - vector_projected[0]
                vel_residual[i+1] = velocity_fluid[1] - vector_projected[1]
                i+=2

        elif self.domain_size == 3:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                velocity_fluid = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                vector_projected = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED,0)

                vel_residual[i] = velocity_fluid[0] - vector_projected[0]
                vel_residual[i+1] = velocity_fluid[1] - vector_projected[1]
                vel_residual[i+2] = velocity_fluid[2] - vector_projected[2]
                i+=3

        # Compute the residual norm and store it in the ProcessInfo
        res_norm = space.TwoNorm(vel_residual)
        self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.FSI_INTERFACE_RESIDUAL_NORM] = res_norm

        return vel_residual


    # Auxiliar function to compute the L2 norm of the difference between the fluid velocity and the mesh velocity at the interface
    def _ComputeMeshResidual(self):

        i = 0
        mesh_res_norm = 0.0

        if self.domain_size == 2:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                velocity_fluid = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                mesh_velocity = node.GetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY,0)

                mesh_res_norm += (velocity_fluid[0]-mesh_velocity[0])**2 + (velocity_fluid[1]-mesh_velocity[1])**2

        elif self.domain_size == 3:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                velocity_fluid = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                mesh_velocity = node.GetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY,0)

                mesh_res_norm += (velocity_fluid[0]-mesh_velocity[0])**2 + (velocity_fluid[1]-mesh_velocity[1])**2 + (velocity_fluid[2]-mesh_velocity[2])**2

        return mesh_res_norm**0.5


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
                a_n1[0] = (v_n1[0] - v_n[0] - self.time_step_fluid*(1-gamma*(alpha-1))*a_n[0]) / ((1-alpha)*self.time_step_fluid*gamma)
                a_n1[1] = (v_n1[1] - v_n[1] - self.time_step_fluid*(1-gamma*(alpha-1))*a_n[1]) / ((1-alpha)*self.time_step_fluid*gamma)
                a_n1[2] = 0.0

                # Compute the current displacement associated to the corrected interface velocity with the Bossak formulaes
                u_n1 = KratosMultiphysics.Vector(3)
                u_n1[0] = u_n[0] + self.time_step_fluid*v_n[0] + (self.time_step_fluid**2)*(0.5-beta)*a_n[0] + (self.time_step_fluid**2)*beta*((1-alpha)*a_n1[0] + alpha*a_n[0])
                u_n1[1] = u_n[1] + self.time_step_fluid*v_n[1] + (self.time_step_fluid**2)*(0.5-beta)*a_n[1] + (self.time_step_fluid**2)*beta*((1-alpha)*a_n1[1] + alpha*a_n[1])
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
                a_n1[0] = (v_n1[0] - v_n[0] - self.time_step_fluid*(1-gamma*(alpha-1))*a_n[0]) / ((1-alpha)*self.time_step_fluid*gamma)
                a_n1[1] = (v_n1[1] - v_n[1] - self.time_step_fluid*(1-gamma*(alpha-1))*a_n[1]) / ((1-alpha)*self.time_step_fluid*gamma)
                a_n1[2] = (v_n1[2] - v_n[2] - self.time_step_fluid*(1-gamma*(alpha-1))*a_n[2]) / ((1-alpha)*self.time_step_fluid*gamma)

                # Compute the current displacement associated to the corrected interface velocity with the Bossak formulaes
                u_n1 = KratosMultiphysics.Vector(3)
                u_n1[0] = u_n[0] + self.time_step_fluid*v_n[0] + (self.time_step_fluid**2)*(0.5-beta)*a_n[0] + (self.time_step_fluid**2)*beta*((1-alpha)*a_n1[0] + alpha*a_n[0])
                u_n1[1] = u_n[1] + self.time_step_fluid*v_n[1] + (self.time_step_fluid**2)*(0.5-beta)*a_n[1] + (self.time_step_fluid**2)*beta*((1-alpha)*a_n1[1] + alpha*a_n[1])
                u_n1[2] = u_n[2] + self.time_step_fluid*v_n[2] + (self.time_step_fluid**2)*(0.5-beta)*a_n[2] + (self.time_step_fluid**2)*beta*((1-alpha)*a_n1[2] + alpha*a_n[2])

                # Set the obtained corrected interface displacement
                node.SetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,0,u_n1)
