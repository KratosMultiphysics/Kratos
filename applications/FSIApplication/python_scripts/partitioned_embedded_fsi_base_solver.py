from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from math import sqrt   # Import the square root from python library

# Import utilities
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural # Import the structure Python solvers wrapper
from KratosMultiphysics.FSIApplication import fsi_coupling_interface                            # Import the FSI coupling interface utility
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory                   # Import the FSI convergence accelerator factory

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

def CreateSolver(model, project_parameters):
    return PartitionedEmbeddedFSIBaseSolver(model, project_parameters)

#TODO: AT THIS POINT, THE EMBEDDED SOLVER COULD DERIVE FROM THE BODY FITTED BASE ONE
class PartitionedEmbeddedFSIBaseSolver(PythonSolver):

    def __init__(self, model, project_parameters):
        # TODO: Remove this as soon as the MPCs are implemented in MPI
        # This has to be done prior to the defaults check to avoid the structural solver to throw an error in MPI
        if not project_parameters["structure_solver_settings"].Has("multi_point_constraints_used"):
            project_parameters["structure_solver_settings"].AddEmptyValue("multi_point_constraints_used")
            project_parameters["structure_solver_settings"]["multi_point_constraints_used"].SetBool(False)

        # Call the base Python solver constructor
        super().__init__(model, project_parameters)

        # Auxiliar variables
        self.parallel_type = self.settings["parallel_type"].GetString()
        coupling_settings = self.settings["coupling_settings"]
        self.max_nl_it = coupling_settings["nl_max_it"].GetInt()
        self.nl_tol = coupling_settings["nl_tol"].GetDouble()
        self.structure_interface_submodelpart_name = coupling_settings["structure_interfaces_list"][0].GetString()

        # Construct the structure solver
        self.structure_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structure_solver_settings"], self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Structure solver construction finished')

        # Construct the fluid solver
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"], self.parallel_type)
        self.level_set_type = self.settings["fluid_solver_settings"]["formulation"]["level_set_type"].GetString()

        # First call to create the embedded intersections model part
        self.__GetEmbedddedSkinUtilityModelPart()

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Fluid solver construction finished')
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Partitioned embedded FSI base solver construction finished')

    @classmethod
    def GetDefaultParameters(cls):

        # Note that only the coupling settings are validated
        # The subdomain solver settings will be validated while instantiating these
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level": 0,
            "parallel_type": "OpenMP",
            "solver_type": "partitioned_embedded",
            "coupling_scheme": "dirichlet_neumann",
            "structure_solver_settings": {
            },
            "fluid_solver_settings":{
            },
            "coupling_settings":{
                "coupling_strategy_settings": {
                    "abs_cut_off_tol": 1e-06,
                    "solver_type": "MVQN",
                    "w_0": 0.5
                },
                "nl_max_it": 30,
                "nl_tol": 1e-07,
                "structure_interfaces_list": []
            }
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def ValidateSettings(self):
        default_settings = self.GetDefaultParameters()

        ## Base class settings validation
        super().ValidateSettings()

        ## Validate coupling settings
        self.settings["coupling_settings"].ValidateAndAssignDefaults(default_settings["coupling_settings"])

    def GetMinimumBufferSize(self):
        buffer_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_structure = self.structure_solver.GetMinimumBufferSize()
        return max(buffer_structure,buffer_fluid)

    def AddVariables(self):
        # Fluid and structure solvers variables addition
        self.fluid_solver.AddVariables()
        self.structure_solver.AddVariables()

    def ImportModelPart(self):
        # Fluid and structure solvers ImportModelPart() call
        self.fluid_solver.ImportModelPart()
        self.structure_solver.ImportModelPart()

    def PrepareModelPart(self):
        # Fluid and structure solvers PrepareModelPart() call
        self.fluid_solver.PrepareModelPart()
        self.structure_solver.PrepareModelPart()

        # FSI interface coupling interfaces initialization
        # The __GetFSICouplingInterfaceStructure is supposed to construct the FSI coupling structure interface in here
        self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart()
        self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart()

    def AddDofs(self):
        # Add DOFs structure
        self.structure_solver.AddDofs()
        # Add DOFs fluid
        self.fluid_solver.AddDofs()

    def Initialize(self):
        # Coupling utility initialization
        # The _GetConvergenceAccelerator is supposed to construct the convergence accelerator in here
        self._GetConvergenceAccelerator().Initialize()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Compute the fluid domain NODAL_AREA values
        # Required by the parallel distance calculator if the distance has to be extended
        if (self.level_set_type == "continuous"):
            KratosMultiphysics.CalculateNodalAreaProcess(self.GetFluidComputingModelPart(), self.__GetDomainSize()).Execute()

        # Initialize the Dirichlet-Neumann interface
        self.__InitializeFSIInterfaces()

        # Initialize the distance field
        update_distance_process = True
        self.__GetDistanceToSkinProcess(update_distance_process).Execute()
        if (self.level_set_type == "continuous"):
            self.__ExtendLevelSet()

        # Initialize the embedded skin utility
        self.__GetEmbeddedSkinUtility()

        with open("FSI_iterations.txt",'w') as f:
            f.write("{}\t{}\t{}\n".format("Step","It.", "err_u"))
            f.close()

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', "Finished initialization.")

    def AdvanceInTime(self, current_time):
        fluid_new_time = self.fluid_solver.AdvanceInTime(current_time)
        structure_new_time = self.structure_solver.AdvanceInTime(current_time)

        self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart().GetRootModelPart().CloneTimeStep(fluid_new_time)
        self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart().GetRootModelPart().CloneTimeStep(structure_new_time)

        if abs(fluid_new_time - structure_new_time) > 1e-12:
            err_msg =  'Fluid new time is: ' + str(fluid_new_time) + '\n'
            err_msg += 'Structure new time is: ' + str(structure_new_time) + '\n'
            err_msg += 'No substepping has been implemented yet. Fluid and structure time step must coincide.'
            raise Exception(err_msg)

        return fluid_new_time

    def InitializeSolutionStep(self):
        # Initialize solution step of fluid, structure and coupling solvers
        self.fluid_solver.InitializeSolutionStep()
        self.structure_solver.InitializeSolutionStep()
        self._GetConvergenceAccelerator().InitializeSolutionStep()

    def Predict(self):
        # Structure solver prediction. It is important to firstly perform the structure
        # prediction to update the current buffer position before the FM-ALE operations.
        # Otherwise position 0 and 1 of the buffer coincide since the advance in time
        # has been already performed but no update has been done yet. Besides, this will
        # give a better approximation of the level-set position at the end of step.
        self.structure_solver.Predict()

        # Update the level set position with the structure prediction
        self.__UpdateLevelSet()

        # Correct the updated level set
        self.fluid_solver.GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Fluid solver prediction
        self.fluid_solver.Predict()

        # Restore the fluid node fixity to its original status
        self.fluid_solver.GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

    def GetComputingModelPart(self):
        err_msg =  'Calling GetComputingModelPart() method in a partitioned solver.\n'
        err_msg += 'Specify the domain of interest by calling:\n'
        err_msg += '\t- GetFluidComputingModelPart()\n'
        err_msg += '\t- GetStructureComputingModelPart()\n'
        raise Exception(err_msg)

    def GetFluidComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def GetStructureComputingModelPart(self):
        return self.structure_solver.GetComputingModelPart()

    def GetStructureSkinModelPart(self):
        return self.model.GetModelPart(self.__GetStructureInterfaceModelPartName())

    def GetStructureSkinElementBasedModelPart(self):
        # Create an auxiliar model part to save the element based skin
        element_based_skin_model_part_name = self.__GetStructureInterfaceModelPartName() + "ElementBased"
        if self.model.HasModelPart(element_based_skin_model_part_name):
            self.model.DeleteModelPart(element_based_skin_model_part_name)
        self.element_based_skin_model_part = self.model.CreateModelPart(element_based_skin_model_part_name)

        # Copy the skin model part conditions to an auxiliar model part elements.
        # This is required for the computation of the distance function, which
        # takes the elements of the second modelpart as skin. If this operation
        # is not performed, no elements are found, yielding a wrong level set.
        self.__GetPartitionedFSIUtilities().CopySkinToElements(
            self.GetStructureSkinModelPart(),
            self.element_based_skin_model_part)

        return self.element_based_skin_model_part

    def GetStructureIntersectionsModelPart(self):
        if not hasattr(self, '_embedded_intersections_model_part'):
            embedded_intersections_root_part = self.model.CreateModelPart("EmbeddedIntersectionsModelPart")
            self._embedded_intersections_model_part = embedded_intersections_root_part.CreateSubModelPart("SkinEmbeddedIntersectionsModelPart")
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)
        return self._embedded_intersections_model_part

    ##TODO: BODY FITTED COPIED AND PASTED
    def SolveSolutionStep(self):
        ## Safe ward to avoid coupling if fluid problem is not initialized
        if not self.fluid_solver._TimeBufferIsInitialized():
            return True

        ## Initialize residual
        dis_residual_norm = self.__ComputeInitialResidual()

        ## FSI nonlinear iteration loop
        nl_it = 0
        is_converged = False
        while (nl_it < self.max_nl_it):
            KratosMultiphysics.Logger.PrintInfo('PartitionedFSIBaseSolver', 'FSI non-linear iteration = {0}'.format(nl_it))
            # Check convergence
            if self.__CheckFSIConvergence(dis_residual_norm):
                is_converged = True
                break

            # Perform the displacement convergence accelerator update
            self._GetConvergenceAccelerator().InitializeNonLinearIteration()
            if not self._GetConvergenceAccelerator().IsBlockNewton():
                self.__GetFSICouplingInterfaceStructure().Update()
            else:
                self.__GetFSICouplingInterfaceStructure().UpdateDisplacement()

            # Update the structure interface position
            self.__GetFSICouplingInterfaceStructure().UpdatePosition(KratosMultiphysics.RELAXED_DISPLACEMENT)

            # Map the RELAXED_DISP from the structure FSI coupling interface to fluid FSI coupling interface
            # Note that we take advance of the fact that the coupling interfaces coincide in the embedded case
            # Then update the fluid FSI coupling interface position
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.RELAXED_DISPLACEMENT,
                KratosMultiphysics.DISPLACEMENT,
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                0)
            self.__GetFSICouplingInterfaceFluid().UpdatePosition(KratosMultiphysics.DISPLACEMENT)

            # Update the EMBEDDED_VELOCITY and solve the fluid problem
            self.__SolveFluid()

            # Transfer the fluid load from the fluid computational mesh to the fluid FSI coupling interface
            # This operation contains all the operations required to interpolate from the embedded CFD mesh
            self.__CalculateFluidInterfaceTraction()

            # Transfer the fluid load to the structure FSI coupling interface
            self.__MapFluidInterfaceTraction()     

            # Compute the current iteration traction
            if not self._GetConvergenceAccelerator().IsBlockNewton():
                # Directly send the map load from the structure FSI coupling interface to the parent one
                self.__GetFSICouplingInterfaceStructure().TransferValuesToFatherModelPart(self.__GetTractionVariable())
            else:
                # Perform the traction convergence accelerator update
                self.__GetFSICouplingInterfaceStructure().UpdateTraction()
                KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                    KratosMultiphysics.RELAXED_TRACTION,
                    self.__GetTractionVariable(),
                    self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                    self.__GetFSICouplingInterfaceStructure().GetFatherModelPart(),
                    0)
            self._GetConvergenceAccelerator().FinalizeNonLinearIteration()

            # Solve the structure problem
            self.__SolveStructure()

            # Compute the residual vector
            dis_residual_norm = self.__GetFSICouplingInterfaceStructure().ComputeResidualVector()
            
            # Check and update iterator counter
            if nl_it == self.max_nl_it:
                KratosMultiphysics.Logger.PrintInfo('PartitionedFSIBaseSolver', 'FSI non-linear converged not achieved in {0} iterations'.format(self.max_nl_it))
            else:
                nl_it += 1

        with open("FSI_iterations.txt",'a') as f:
            f.write("{0}\t{1}\t{2}\n".format(self.fluid_solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP],nl_it, dis_residual_norm))
            f.close()

        return is_converged

    # #TODO: OLD EMBEDDED SOLVER
    # def SolveSolutionStep(self):
    #     ## Safe ward to avoid coupling if fluid problem is not initialized
    #     if not self.fluid_solver._TimeBufferIsInitialized():
    #         return True

    #     ## FSI nonlinear iteration loop
    #     nl_it = 0
    #     while (nl_it < self.max_nl_it and self.fluid_solver._TimeBufferIsInitialized()):
    #         KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\tFSI non-linear iteration = {0}'.format(nl_it))


    #         #FIXME: I THINK THAT THE BLOCK BETWEEN TODO MUST BE DONE ONLY IF nl_it > 0
    #         #TODO: CHECK WHAT HAPPENS FOR nl_it == 0
    #         #TODO: CHECK WHAT HAPPENS FOR nl_it == 0
    #         # Update the structure interface position #FIXME:  THIS CopyVariable WAS ONLY REQUIRED TO DO THE UPDATE. WITH THIS NEW VERSION SHOULD BE EQUIVALENT AND MORE EFFICIENT
    #         # KratosMultiphysics.VariableUtils().CopyVariable(
    #         #     KratosMultiphysics.RELAXED_DISPLACEMENT,
    #         #     KratosMultiphysics.DISPLACEMENT,
    #         #     self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart().Nodes)
    #         # self.__GetFSICouplingInterfaceStructure().UpdatePosition()
    #         self.__GetFSICouplingInterfaceStructure().UpdatePosition(RELAXED_DISPLACEMENT) 

    #         # Map the RELAXED_DISP from the structure FSI coupling interface to fluid FSI coupling interface
    #         # Note that we take advance of the fact that the coupling interfaces coincide in the embedded case
    #         # Then update the fluid FSI coupling interface position
    #         KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
    #             KratosMultiphysics.RELAXED_DISPLACEMENT,
    #             KratosMultiphysics.DISPLACEMENT,
    #             self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
    #             self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
    #             0)
    #         self.__GetFSICouplingInterfaceFluid().UpdatePosition()
    #         #TODO: CHECK WHAT HAPPENS FOR nl_it == 0
    #         #TODO: CHECK WHAT HAPPENS FOR nl_it == 0

    #         # Update the EMBEDDED_VELOCITY and solve the fluid problem
    #         self.__SolveFluid()

    #         # Perform the pressure convergence accelerator update
    #         p_residual_norm = 0 #TODO: REMOVE WHEN REMOVING THE TEMPORARY ITERATION COUNTER
    #         if (self.traction_update and nl_it > 0):
    #             # Compute the pressure residual vector
    #             p_residual_norm = self.__GetFSICouplingInterfaceFluid().ComputeResidualVector()
    #             KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\t|p_res| = {}'.format(p_residual_norm))

    #             # Call the fluid FSI coupling interface update
    #             self.__GetTractionConvergenceAccelerator().InitializeNonLinearIteration()
    #             self.__GetFSICouplingInterfaceFluid().Update()
    #             self.__GetTractionConvergenceAccelerator().FinalizeNonLinearIteration()

    #             # Update the fluid interface position
    #             #FIXME: THIS SHOULD BE COMPATIBLE WITH TWO INTERFACES
    #             if (self.level_set_type == "continuous"):
    #                 KratosMultiphysics.VariableUtils().CopyVariable(
    #                 KratosMultiphysics.RELAXED_SCALAR,
    #                 KratosMultiphysics.POSITIVE_FACE_PRESSURE,
    #                 self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart().Nodes)
    #             elif (self.level_set_type == "discontinuous"):
    #                 #TODO: THINK ABOUT THE UPDATE OF POSITIVE AND NEGATIVE FACE PRESSURES
    #                 raise Exception("NOT IMPLEMENTED YET!!!!")
    #             else:
    #                 err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
    #                 raise Exception(err_msg)

    #         # Transfer the fluid load to the structure FSI coupling interface
    #         self.__MapFluidInterfaceTraction()                

    #         # Solve the structure problem
    #         self.__SolveStructure()

    #         # Compute the residual vector
    #         dis_residual_norm = self.__GetFSICouplingInterfaceStructure().ComputeResidualVector()

    #         # Check convergence
    #         if self.__CheckFSIConvergence(dis_residual_norm):
    #             with open("FSI_iterations.txt",'a') as f:
    #                 f.write("{0:<33}{1:^33}{2:^33}{3:>33}\n".format(self.fluid_solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP],nl_it, dis_residual_norm, p_residual_norm))
    #                 f.close()
    #             return True
    #         else:
    #             # Perform the displacement convergence accelerator update
    #             self._GetConvergenceAccelerator().InitializeNonLinearIteration() #TODO: This call must be done within the fsi_coupling_interface
    #             self.__GetFSICouplingInterfaceStructure().Update()
    #             self._GetConvergenceAccelerator().FinalizeNonLinearIteration() #TODO: This call must be done within the fsi_coupling_interface
                
    #             # Update iterator counter
    #             nl_it += 1

    #     with open("FSI_iterations.txt",'a') as f:
    #         f.write("{0:<33}{1:^33}{2:^33}{3:>33}\n".format(self.fluid_solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP],nl_it, dis_residual_norm, p_residual_norm))
    #         f.close()

    #     # Maximum iterations reached without convergence
    #     return False

    def FinalizeSolutionStep(self):
        # Finalize solution step
        self.fluid_solver.FinalizeSolutionStep()
        self.structure_solver.FinalizeSolutionStep()
        self._GetConvergenceAccelerator().FinalizeSolutionStep()

    def Finalize(self):
        self.fluid_solver.Finalize()
        self.structure_solver.Finalize()
        self._GetConvergenceAccelerator().Finalize()

    def SetEchoLevel(self, structure_echo_level, fluid_echo_level):
        self.fluid_solver.SetEchoLevel(self, fluid_echo_level)
        self.structure_solver.SetEchoLevel(self, structure_echo_level)

    def Clear(self):
        self.fluid_solver.Clear()
        self.structure_solver.Clear()

    def Check(self):
        self.fluid_solver.Check()
        self.structure_solver.Check()

    #######################################################################
    ##############          PRIVATE METHODS SECTION          ##############
    #######################################################################

    def __GetDistanceToSkinProcess(self, update_distance_process = False):
        if update_distance_process:
            self._distance_to_skin_process = self.__CreateDistanceToSkinProcess()
        return self._distance_to_skin_process

    def __CreateDistanceToSkinProcess(self):
        # Set the distance computation process
        if (self.level_set_type == "continuous"):
            raycasting_relative_tolerance = 1.0e-10
            if self.__GetDomainSize() == 2:
                return KratosMultiphysics.CalculateDistanceToSkinProcess2D(
                    self.GetFluidComputingModelPart(),
                    self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                    raycasting_relative_tolerance)
            elif self.__GetDomainSize() == 3:
                return KratosMultiphysics.CalculateDistanceToSkinProcess3D(
                    self.GetFluidComputingModelPart(),
                    self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                    raycasting_relative_tolerance)
            else:
                raise Exception("Domain size expected to be 2 or 3. Got " + str(self.__GetDomainSize()))
        elif (self.level_set_type == "discontinuous"):
            if self.__GetDomainSize() == 2:
                return KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess2D(
                    self.GetFluidComputingModelPart(),
                    self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart())
            elif self.__GetDomainSize() == 3:
                return KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess3D(
                    self.GetFluidComputingModelPart(),
                    self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart())
            else:
                raise Exception("Domain size expected to be 2 or 3. Got " + str(self.__GetDomainSize()))
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    def __GetParallelDistanceCalculator(self):
        if not hasattr(self, '_parallel_distance_calculator'):
            self._parallel_distance_calculator = self.__CreateParallelDistanceCalculator()
        return self._parallel_distance_calculator

    def __CreateParallelDistanceCalculator(self):
        if self.__GetDomainSize() == 2:
            return KratosMultiphysics.ParallelDistanceCalculator2D()
        elif self.__GetDomainSize() == 3:
            return KratosMultiphysics.ParallelDistanceCalculator3D()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.__GetDomainSize()))

    def __GetEmbeddedSkinUtility(self):
        if not hasattr(self, '_embedded_skin_utility'):
            self._embedded_skin_utility = self.__CreateEmbeddedSkinUtility()
        return self._embedded_skin_utility

    def __CreateEmbeddedSkinUtility(self):
        if self.__GetDomainSize() == 2:
            return KratosMultiphysics.EmbeddedSkinUtility2D(
                self.GetFluidComputingModelPart(),
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.level_set_type)
        elif self.__GetDomainSize() == 3:
            return KratosMultiphysics.EmbeddedSkinUtility3D(
                self.GetFluidComputingModelPart(),
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.level_set_type)
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.__GetDomainSize()))

    def __GetEmbedddedSkinUtilityModelPart(self):
        if not hasattr(self, '_embedded_skin_utility_model_part'):
            self._embedded_skin_utility_model_part = self.__CreateEmbeddedSkinUtilityModelPart()
        return self._embedded_skin_utility_model_part

    def __CreateEmbeddedSkinUtilityModelPart(self):
        embedded_skin_utility_skin_model_part_name = "EmbeddedSkinUtilityModelPart"
        embedded_skin_utility_skin_model_part =self.model.CreateModelPart(embedded_skin_utility_skin_model_part_name)
        embedded_skin_utility_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        if (self.level_set_type == "continuous"):
            embedded_skin_utility_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        elif (self.level_set_type == "discontinuous"):
            embedded_skin_utility_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
            embedded_skin_utility_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)
        return embedded_skin_utility_skin_model_part

    def __GetStructureInterfaceModelPartName(self):
        str_int_list = self.settings["coupling_settings"]["structure_interfaces_list"]
        if (str_int_list.size() != 1):
            raise Exception("FSI embedded solver structure skin must be contained in a unique model part")
        return str_int_list[0].GetString()

    #TODO: CHECK WHY THIS IS NEEDED...
    #TODO: THINK WE SHOULD ADD THE CONDITIONS IN HERE... RIGHT NOW WE ASSUME THIS ARE PROVIDED IN THE MDPA
    def __InitializeFSIInterfaces(self):
        # Initialize Neumann structure interface
        str_interface_submodelpart = self.model.GetModelPart(self.__GetStructureInterfaceModelPartName())

        # Set the INTERFACE flag to the structure skin
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, str_interface_submodelpart.Nodes)

    def __UpdateLevelSet(self):
        # Recompute the distance field with the obtained solution
        self.__GetDistanceToSkinProcess().Execute()

        # Extend the level set to the first layer of non-intersected elements
        # This is required in case the distance modification process moves the level set
        # to a non-intersected element to prevent almost empty fluid elements.
        # Note that non-intersected elements have a large default distance value, which might
        # alter the zero isosurface when the distance modification avoids almost empty elements.
        if (self.level_set_type == "continuous"):
            self.__ExtendLevelSet()

    def __ExtendLevelSet(self):
        max_layers = 2
        max_distance = 1.0e+12
        self.__GetParallelDistanceCalculator().CalculateDistances(
            self.GetFluidComputingModelPart(),
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.NODAL_AREA,
            max_layers,
            max_distance)

    def __SolveFluid(self):
        # Update the current iteration level-set position
        self.__UpdateLevelSet()

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep() # This contains the FM-ALE operations and the level set correction

    def __CalculateFluidInterfaceTraction(self):
        if (self.level_set_type == "continuous"):
            # Interpolate the pressure to the fluid FSI coupling interface
            self.__GetPartitionedFSIUtilities().EmbeddedPressureToPositiveFacePressureInterpolator(
                self.GetFluidComputingModelPart(),
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart())

        elif (self.level_set_type == "discontinuous"):
            # Generate the intersections skin to map from
            self.__GetEmbeddedSkinUtility().GenerateSkin()

            # Interpolate POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the background mesh
            self.__GetEmbeddedSkinUtility().InterpolateDiscontinuousMeshVariableToSkin(
                KratosMultiphysics.PRESSURE, KratosMultiphysics.POSITIVE_FACE_PRESSURE, "positive")
            self.__GetEmbeddedSkinUtility().InterpolateDiscontinuousMeshVariableToSkin(
                KratosMultiphysics.PRESSURE, KratosMultiphysics.NEGATIVE_FACE_PRESSURE, "negative")

        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    def __MapFluidInterfaceTraction(self):
        if (self.level_set_type == "continuous"):
            # Map PRESSURE from fluid FSI coupling interface to structure FSI coupling interface
            # Note that in here we take advantage of the fact that the coupling interfaces coincide in the embedded case
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)

            # Convert the pressure scalar load to a traction vector one
            self.__GetPartitionedFSIUtilities().CalculateTractionFromPressureValues(
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                self.__GetTractionVariable())

        elif (self.level_set_type == "discontinuous"):
            # Map the POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the auxiliary embedded skin model part,
            # which is created from the elemental level set intersections, to the fluid FSI coupling interface
            mapper_params = KratosMultiphysics.Parameters("""{
                "mapper_type": "nearest_element",
                "echo_level" : 0
            }""")
            mapper = KratosMapping.MapperFactory.CreateMapper(
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                mapper_params)
            mapper.Map(KratosMultiphysics.POSITIVE_FACE_PRESSURE, KratosMultiphysics.POSITIVE_FACE_PRESSURE)
            mapper.Map(KratosMultiphysics.NEGATIVE_FACE_PRESSURE, KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

            # Transfer POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the fluid coupling interface to the structure one
            # Note that in here we take advantage of the fact that the coupling interfaces coincide in the embedded case
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.NEGATIVE_FACE_PRESSURE,
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)

            # Convert the pressure scalar load to a traction vector one
            self.__GetPartitionedFSIUtilities().CalculateTractionFromPressureValues(
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                KratosMultiphysics.NEGATIVE_FACE_PRESSURE,
                self.__GetTractionVariable())

        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    #TODO: CHECK IF WE COULD DERIVE THIS FROM THE BASE PARTITIONED SOLVER AND USE IT THIS METHOD IN HERE (IS EQUAL IN BODY FITTED)
    #TODO: THE UNIQUE DIFFERENCE IS THAT THE BODY-FITTED ONE DOES THE TRACTION TO POINT LOAD CONVERSION
    def __SolveStructure(self):
        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()

        # Transfer the obtained DISPLACEMENT to the structure FSI coupling interface
        # Note that this are the current non-linear iteration unrelaxed values (\tilde{u}^{k+1})
        # These values will be employed to calculate the interface residual vector
        # #TODO: REMOVE THIS. IT IS ALREADY DONE IN THE FSI COUPLING INTERFACE
        # self.__GetFSICouplingInterfaceStructure().GetValuesFromFatherModelPart(KratosMultiphysics.DISPLACEMENT)

    #TODO: MOVE THIS TO A BASE CLASS
    def __CheckFSIConvergence(self, residual_norm):
        interface_dofs = self.__GetPartitionedFSIUtilities().GetInterfaceResidualSize(self.__GetStructureInterfaceSubmodelPart())
        normalised_residual = residual_norm/sqrt(interface_dofs)
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\t|res|/sqrt(nDOFS) = ' + str(normalised_residual))
        return normalised_residual < self.nl_tol

    # This method returns the convergence accelerator.
    # If it is not created yet, it calls the __CreateConvergenceAccelerator first
    def _GetConvergenceAccelerator(self):
        if not hasattr(self, '_convergence_accelerator'):
            self._convergence_accelerator = self.__CreateConvergenceAccelerator()
        return self._convergence_accelerator

    # This method returns the fluid interface convergence accelerator.
    # If it is not created yet, it calls the __CreateConvergenceAccelerator first
    def __GetTractionConvergenceAccelerator(self):
        if not hasattr(self, '_traction_convergence_accelerator'):
            self._traction_convergence_accelerator = self.__CreateConvergenceAccelerator()
        return self._traction_convergence_accelerator

    # This method constructs the convergence accelerator coupling utility
    def __CreateConvergenceAccelerator(self):
        convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self.settings["coupling_settings"]["coupling_strategy_settings"])
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Coupling strategy construction finished')
        return convergence_accelerator

    def __GetFSICouplingInterfaceStructure(self):
        if not hasattr(self, '_fsi_coupling_interface_structure'):
            self._fsi_coupling_interface_structure = self.__CreateFSICouplingInterfaceStructure()
        return self._fsi_coupling_interface_structure

    def __CreateFSICouplingInterfaceStructure(self):
        # Set auxiliary settings
        if (self.level_set_type == "continuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceStructure",
                "parent_model_part_name": "",
                "input_variable_list": [],
                "output_variable_list": ["DISPLACEMENT"],
                "auxiliary_variable_list": ["POSITIVE_FACE_PRESSURE"]
            }""")
        elif (self.level_set_type == "discontinuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceStructure",
                "parent_model_part_name": "",
                "input_variable_list": [],
                "output_variable_list": ["DISPLACEMENT"],
                "auxiliary_variable_list": ["POSITIVE_FACE_PRESSURE","NEGATIVE_FACE_PRESSURE"]
            }""")
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

        aux_settings["parent_model_part_name"].SetString(self.structure_interface_submodelpart_name)
        aux_settings["input_variable_list"].Append(self.__GetTractionVariable().Name())

        # Construct the FSI coupling interface
        fsi_coupling_interface_structure = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings,
            self._GetConvergenceAccelerator())

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Structure FSI coupling interface created')

        return fsi_coupling_interface_structure

    def __GetFSICouplingInterfaceFluid(self):
        if not hasattr(self, '_fsi_coupling_interface_fluid'):
            self._fsi_coupling_interface_fluid = self.__CreateFSICouplingInterfaceFluid()
        return self._fsi_coupling_interface_fluid

    def __CreateFSICouplingInterfaceFluid(self):
        # Set auxiliary settings
        # Note that in the embedded case, the fluid interface is identical to the structure one
        # This is intentionally done, since this copy will be used in the level set computation
        if (self.level_set_type == "continuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceFluid",
                "parent_model_part_name": "",
                "input_variable_list": ["DISPLACEMENT"],
                "output_variable_list": ["POSITIVE_FACE_PRESSURE"]
            }""")
        elif (self.level_set_type == "discontinuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceFluid",
                "parent_model_part_name": "",
                "input_variable_list": ["DISPLACEMENT"],
                "output_variable_list": ["POSITIVE_FACE_PRESSURE","NEGATIVE_FACE_PRESSURE"]
            }""")
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

        aux_settings["parent_model_part_name"].SetString(self.structure_interface_submodelpart_name)

        # Construct the FSI coupling interface
        fsi_coupling_interface_fluid = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings)

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Fluid FSI coupling interface created')

        return fsi_coupling_interface_fluid

    #TODO: THIS SHOULD BE IMPLEMENTED IN AN EVENTUAL BASE CLASS
    def __GetStructureInterfaceSubmodelPart(self):
        # Returns the structure interface submodelpart that will be used in the residual minimization
        return self.model.GetModelPart(self.structure_interface_submodelpart_name)

    #TODO: MOVE TO EVENTUAL BASE SOLVER
    #TODO: SHOULDN'T ADD THIS TO THE BASE PYTHON SOLVER
    def __GetDomainSize(self):
        if not hasattr(self, 'domain_size'):
            fluid_domain_size = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            structure_domain_size = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            if fluid_domain_size !=structure_domain_size:
                raise("ERROR: Solid domain size and fluid domain size are not equal!")
            self.domain_size = fluid_domain_size
        return self.domain_size

    def __ComputeInitialResidual(self):
        # Save as RELAXED_DISPLACEMENT the DISPLACEMENT coming from the structure Predict()
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            KratosMultiphysics.DISPLACEMENT,
            KratosMultiphysics.RELAXED_DISPLACEMENT,
            self.__GetStructureInterfaceSubmodelPart(),
            self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
            0)

        # Update the structure interface position with the DISPLACEMENT values from the predict
        self.__GetFSICouplingInterfaceStructure().UpdatePosition(KratosMultiphysics.RELAXED_DISPLACEMENT)

        # Map the RELAXED_DISP from the structure FSI coupling interface to fluid FSI coupling interface
        # Note that we take advance of the fact that the coupling interfaces coincide in the embedded case
        # Then update the fluid FSI coupling interface position
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            KratosMultiphysics.RELAXED_DISPLACEMENT,
            KratosMultiphysics.DISPLACEMENT,
            self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
            self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
            0)
        self.__GetFSICouplingInterfaceFluid().UpdatePosition()

        # Update the EMBEDDED_VELOCITY and solve the fluid problem
        self.__SolveFluid()

        # Transfer the fluid load from the fluid computational mesh to the fluid FSI coupling interface
        # This operation contains all the operations required to interpolate from the embedded CFD mesh
        self.__CalculateFluidInterfaceTraction()

        # Transfer the fluid load to the structure FSI coupling interface
        self.__MapFluidInterfaceTraction()

        # Save as RELAXED_TRATION the TRACTION coming from the fluid
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            self.__GetTractionVariable(),
            KratosMultiphysics.RELAXED_TRACTION,
            self.__GetStructureInterfaceSubmodelPart(),
            self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
            0)

        # Directly send the map load from the structure FSI coupling interface to the parent one
        self.__GetFSICouplingInterfaceStructure().TransferValuesToFatherModelPart(self.__GetTractionVariable())                    

        # Solve the structure problem
        self.__SolveStructure()

        # Compute the residual vector
        dis_residual_norm = self.__GetFSICouplingInterfaceStructure().ComputeResidualVector()

        return dis_residual_norm

    def __GetPartitionedFSIUtilities(self):
        if not hasattr(self, '_partitioned_fsi_utilities'):
            self._partitioned_fsi_utilities = self.__CreatePartitionedFSIUtilities()
        return self._partitioned_fsi_utilities

    def __CreatePartitionedFSIUtilities(self):
        if self.__GetDomainSize() == 2:
            return KratosFSI.PartitionedFSIUtilitiesArray2D()
        elif self.__GetDomainSize() == 3:
            return KratosFSI.PartitionedFSIUtilitiesArray3D()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.__GetDomainSize()))

    def __GetTractionVariable(self):
        if self.__GetDomainSize() == 2:
            return KratosStructural.LINE_LOAD
        elif self.__GetDomainSize() == 3:
            return KratosStructural.SURFACE_LOAD
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.__GetDomainSize()))
        