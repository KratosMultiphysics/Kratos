# Kratos imports
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC
from KratosMultiphysics.restart_utility import RestartUtility

# RigidBody imports
from . import rigid_body_input_check as input_check

# Other imports
import numpy as np
import json
import os


class RigidBodySolver:
    """
    @details This class implements a Rigid Body within Kratos. It can be seen as a combination of 6 single-degree-of-freedom
    (SDOF) solvers, solving independently the 3 displacements and 3 rotations of a rigid body. It is meant to be used
    with CoSimulationApplication (e.g. rigid object submerged in a certain flow). For standalone usage, it can also
    be called directly, since the class acts at the same time as a solver and as an analysis.

    Two sub model parts, each of them with one node in (0,0,0), are directly generated here: "RigidBody" (representing
    the body itself) and "RootPoint" (representing the attachment point). Here is a summary of the variables that each
    sub model part has (note that the sketch represents each of the DOFs):

     ---------------        Sub model part: RigidBody                          \ 
    |               |       Node ID: 1                                         |
    |    Node 1     |       Specific variables:                                |
    |       X       |           FORCE, MOMENT,                                 |
    |    (0,0,0)    |           IMPOSED_FORCE, IMPOSED_MOMENT                  |
    |               |           EFFECTIVE_FORCE, EFFECTIVE_MOMENT              |    Model part: Main
     ---------------            BODY_FORCE, BODY_MOMENT                        |    Nodes IDs: 1, 2
        |       |                                                               >   General variables:
       <_      _|_                                                             |        DISPLACEMENT, ROTATION,
        _>    |___|                                                            |        VELOCITY, ANGULAR_VELOCITY
       <_     | | |                                                            |        ACCELERATION, ANGULAR_ACCELERATION
         >      |           Sub model part: RootPoint                          |
        |       |           Node ID: 2                                         |
    /////// X ////////      Specific variables:                                |
         Node 2                 REACTION, REACTION_MOMENT                      |
         (0,0,0)                IMPOSED_DISPLACEMENT, IMPOSED_ROTATION         /

    @note It is possible to apply forces to the rigid body as well as displacements to the reference point (e.g. to be used
    as a TMD). The "IMPOSED_*" variables have the same behaviour as their original version but are necessary to avoid
    overwriting data in some cases (e.g. when a force comes from another solver with CoSimulation but an extra force must
    be prescribed directly from the project parameters).
    """

    def __init__(self, model, project_parameters):

        # Check that the input parameters are in the right format
        if not isinstance(model, KM.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")
        if not isinstance(project_parameters, KM.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        # Basic check to see if project_parameters has all the necessary data to set up the problem
        input_check._CheckMandatoryInputParameters(project_parameters)

        # Read the basic problem data and store it in class variables
        self.problem_name = project_parameters["problem_data"]["problem_name"].GetString()
        self.start_time = project_parameters["problem_data"]["start_time"].GetDouble()
        self.end_time = project_parameters["problem_data"]["end_time"].GetDouble()
        
        # Check the solver settings, fill empty fields with the default values and save them in class variables
        solver_settings = input_check._ValidateAndAssignRigidBodySolverDefaults(project_parameters["solver_settings"])
        self._InitializeSolutionVariables(solver_settings)
        buffer_size = solver_settings["buffer_size"].GetInt() # we only need it in the constructor

        # The degrees of freedom (DOFs) that can be activated from the parameters depend on the domain size
        # The system size, quantity of linear variables (displacement, force...)
        # and angular variables (rotation, moment...) needs to be adjusted
        if self.domain_size == 3:
            self.available_dofs = ['displacement_x', 'displacement_y', 'displacement_z',
                'rotation_x', 'rotation_y', 'rotation_z']
        elif self.domain_size == 2:
            self.available_dofs = ['displacement_x', 'displacement_y', 'rotation_z']
        self.system_size = len(self.available_dofs)
        self.linear_size = int(np.ceil(self.system_size/2)) # 3 or 2 with a 3D or 2D problem respectively
        self.angular_size = int(np.floor(self.system_size/2)) # 3 or 1 with a 3D or 2D problem respectively

        # Check the settings related to each DOF, fill empty fields and save them in class variables
        dof_settings, self.active_dofs = input_check._ValidateAndAssignDofDefaults(solver_settings["active_dofs"], self.available_dofs)
        self._InitializeDofsVariables(dof_settings)

        # Prepare the coefficients to apply the generalized-alpha method
        self._InitializeGeneralizedAlphaParameters()

        # Create the Kratos model and generate/import the sub model parts with all the necessary variables
        self.model = model
        self.main_model_part = self.model.CreateModelPart("Main")
        # If it is not a restart, the sub model parts must be created from scratch
        if solver_settings["model_import_settings"]["input_type"].GetString() == "none":
            self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.domain_size
            self.rigid_body_model_part = self.main_model_part.CreateSubModelPart("RigidBody")
            self.root_point_model_part = self.main_model_part.CreateSubModelPart("RootPoint")
            self.AddVariables()
            self.rigid_body_model_part.CreateNewNode(1,0.0,0.0,0.0)
            self.root_point_model_part.CreateNewNode(2,0.0,0.0,0.0)
            self.main_model_part.SetBufferSize(buffer_size)
        # If it is a restart, the sub model parts can be imported from the restart file
        elif solver_settings["model_import_settings"]["input_type"].GetString() == "rest":
            model_import_settings = solver_settings["model_import_settings"]
            model_import_settings.RemoveValue("input_type")
            RestartUtility(self.main_model_part, model_import_settings).LoadRestart()
            self.rigid_body_model_part = self.main_model_part.GetSubModelPart("RigidBody")
            self.root_point_model_part = self.main_model_part.GetSubModelPart("RootPoint")
            self.start_time = self.main_model_part.ProcessInfo[KM.TIME]
        
        # Create all the processes stated in the project parameters
        self._list_of_processes = input_check._CreateListOfProcesses(self.model, project_parameters, self.main_model_part)
        self._list_of_output_processes = input_check._CreateListOfOutputProcesses(self.model, project_parameters)
        self._list_of_processes.extend(self._list_of_output_processes)


    def _InitializeSolutionVariables(self, solver_settings):
        
        # Save all the data that does not depend on the degree of freedom
        self.domain_size = solver_settings["domain_size"].GetInt()
        self.echo_level = solver_settings["echo_level"].GetInt()
        self.rho_inf = solver_settings["time_integration_parameters"]["rho_inf"].GetDouble()
        self.delta_t = solver_settings["time_integration_parameters"]["time_step"].GetDouble()


    def _InitializeDofsVariables(self, dof_settings):

        # Initialize variables that depend on the degree of freedom
        self.is_constrained = {}
        self.M = np.zeros((self.system_size,self.system_size)) # Mass matrix
        self.C = np.zeros((self.system_size,self.system_size)) # Damping matrix
        self.K = np.zeros((self.system_size,self.system_size)) # Stiffness matrix

        # Fill the initialised values with the ones from the parameters
        for index, dof in enumerate(self.available_dofs):
            self.is_constrained[dof] = dof_settings[dof]["constrained"].GetBool()
            self.M[index][index] = dof_settings[dof]['system_parameters']['mass'].GetDouble()
            self.C[index][index] = dof_settings[dof]['system_parameters']['damping'].GetDouble()
            self.K[index][index] = dof_settings[dof]['system_parameters']['stiffness'].GetDouble()


    def _InitializeGeneralizedAlphaParameters(self):

        # Parameters of the method (see Chung, 1993)
        self.alpha_f = self.rho_inf / (self.rho_inf + 1)
        self.alpha_m = (2*self.rho_inf - 1) / (self.rho_inf + 1)
        self.beta = 0.25 * (1- self.alpha_m + self.alpha_f)**2
        self.gamma =  0.50 - self.alpha_m + self.alpha_f
        
        # Coefficients for LHS
        self.a1h = (1.0 - self.alpha_m) / (self.beta * self.delta_t**2)
        self.a2h = (1.0 - self.alpha_f) * self.gamma / (self.beta * self.delta_t)
        self.a3h = 1.0 - self.alpha_f
        # Since the input can only be linear, the Left-Hand Side is always
        # the same and does not need to be calculated before each solution step.
        self.LHS = self.a1h * self.M + self.a2h * self.C + self.a3h * self.K

        # Coefficients for the RHS
        # for the mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.delta_t
        self.a3m = (1.0 - self.alpha_m - 2.0 * self.beta) / (2.0 * self.beta)
        # for the damping
        self.a1b = (1.0 - self.alpha_f) * self.gamma / (self.beta * self.delta_t)
        self.a2b = (1.0 - self.alpha_f) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alpha_f) * (0.5 * self.gamma / self.beta - 1.0) * self.delta_t
        # for the stiffness
        self.a1k = -1.0 * self.alpha_f

        # Coefficients to update the velocity
        self.a1v = self.gamma / (self.beta * self.delta_t)
        self.a2v = 1.0 - self.gamma / self.beta
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.delta_t

        # Coefficients to update the acceleration
        self.a1a = self.a1v / (self.delta_t * self.gamma)
        self.a2a = -1.0 / (self.beta * self.delta_t)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)
        

    def AddVariables(self):

        # Kinematic variables (work with both RigidBody and RootPoint model parts)
        self.main_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ANGULAR_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ANGULAR_ACCELERATION)

        # Specific variables for the model part RigidBody
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KM.FORCE)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KM.MOMENT)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KMC.IMPOSED_FORCE)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KMC.IMPOSED_MOMENT)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KMC.EFFECTIVE_FORCE)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KMC.EFFECTIVE_MOMENT)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KM.BODY_FORCE)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KM.BODY_MOMENT)

        # Specific variables for the model part RootPoint
        self.root_point_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.root_point_model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        self.root_point_model_part.AddNodalSolutionStepVariable(KMC.IMPOSED_DISPLACEMENT)
        self.root_point_model_part.AddNodalSolutionStepVariable(KMC.IMPOSED_ROTATION)


    def Run(self):

        # This carries out the whole simulation
        # (standard method, see Kratos analysis stage)
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()


    def Initialize(self):

        # Initialize the time only if it is not a restart.
        # If it is, the simulation will beginn at the time saved in the restart file
        if not self.main_model_part.ProcessInfo[KM.IS_RESTARTED]:
            self.main_model_part.ProcessInfo[KM.TIME] = self.start_time

        # Initialize other variables that are used in SolveSolutionStep()
        self.total_root_point_displ = np.zeros(self.system_size)
        self.total_load = np.zeros(self.system_size)

        # Let the processes do their initialization tasks
        for process in self._list_of_processes:
            process.ExecuteInitialize()

        # The solution loop starts after this method ends, 
        # so let the processes do their necessary tasks
        for process in self._list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        # Standard Kratos output
        if self.echo_level > 0:
            KM.Logger.PrintInfo("Rigid Body Solver", "Analysis -START- ")


    def RunSolutionLoop(self):

        # Solve each time step until the end time is reached
        # (standard method, see Kratos analysis stage)
        while self.main_model_part.ProcessInfo[KM.TIME] < self.end_time:
            self.AdvanceInTime(self.main_model_part.ProcessInfo[KM.TIME])
            self.InitializeSolutionStep()
            self.Predict()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    
    def Finalize(self):
        
        # Let the process do their finalization tasks
        for process in self._list_of_processes:
            process.ExecuteFinalize()

        # Standard Kratos output
        if self.echo_level > 0:
            KM.Logger.PrintInfo("Rigid Body Solver", "Analysis -END- ")


    def AdvanceInTime(self, current_time):

        # Update the step and time of the simulation
        time = current_time + self.delta_t
        # This updates KM.TIME and rotates slides the buffer values
        self.main_model_part.CloneTimeStep(time)
        self.main_model_part.ProcessInfo[KM.STEP] += 1

        # Since the time step is cloned, the variables still have the previous values
        # It is necessary to set to zero the variables that might not be overwritten
        self._ResetExternalVariables()

        # It is necessary to return the time to use the solver with CoSim
        return time

    
    def InitializeSolutionStep(self):

        # Standard Kratos output
        if self.echo_level > 0:
            KM.Logger.PrintInfo("Rigid Body Solver", "STEP: ", self.main_model_part.ProcessInfo[KM.STEP])
            KM.Logger.PrintInfo("Rigid Body Solver", "TIME: ", self.main_model_part.ProcessInfo[KM.TIME])
        
        # Let the processes do their tasks before solving the step
        for process in self._list_of_processes:
            process.ExecuteInitializeSolutionStep()

    
    def Predict(self):
        # In case it is implemented in the future
        pass


    def SolveSolutionStep(self):
        
        # Calculate the effective load which considers both external/imposed
        # loads and the equivalent load from the root point displacement
        self._CalculateEffectiveLoad()
        eff_load = self._GetCompleteVector("rigid_body", KMC.EFFECTIVE_FORCE, KMC.EFFECTIVE_MOMENT)
        eff_load_prev = self._GetCompleteVector("rigid_body", KMC.EFFECTIVE_FORCE, KMC.EFFECTIVE_MOMENT, buffer=1)

        # Calculate the gen-alpha load for the construction of the RHS
        F = (1.0 - self.alpha_f) * eff_load + self.alpha_f * eff_load_prev

        # Construct the right-hand-side (RHS) of the system according to the gen-alpha method
        x_prev, v_prev, a_prev = self._GetKinematics("rigid_body", buffer=1)
        RHS = np.dot(self.M, (self.a1m * x_prev + self.a2m * v_prev + self.a3m * a_prev))
        RHS += np.dot(self.C, (self.a1b * x_prev + self.a2b * v_prev + self.a3b * a_prev))
        RHS += np.dot(self.a1k * self.K, x_prev) + F

        # Force to zero the corrresponding values of RHS so the inactive dofs are not excited
        for index, dof in enumerate(self.available_dofs):
            if dof not in self.active_dofs:
                RHS[index] = 0

        # Solve the solution step and find the new displacements
        x = np.linalg.solve(self.LHS, RHS)

        # Force the constrained dofs to have the same rigid body displacement as in the root point
        for index, dof in enumerate(self.available_dofs):
            if self.is_constrained[dof]:
                x[index] = self.total_root_point_displ[index]

        # Update velocity and acceleration according to the gen-alpha method
        self._UpdateKinematics("rigid_body", x)

        # Update the reaction now that the new displacements have been found
        reaction = self._CalculateReaction()
        self._SetCompleteVector("root_point", KM.REACTION, KM.REACTION_MOMENT, reaction)

    
    def FinalizeSolutionStep(self):
        
        # Let the processes do their tasks after solving the step
        for process in self._list_of_processes:
            process.ExecuteFinalizeSolutionStep()


    def OutputSolutionStep(self):
            
        # Ensure that execute is called only once or never (see Kratos analysis stage)
        execute_was_called = False
        for output_process in self._list_of_output_processes:
            if output_process.IsOutputStep():

                # Let regular processes do their tasks before outputting the step
                if not execute_was_called:
                    for process in self._list_of_processes:
                        process.ExecuteBeforeOutputStep()
                    execute_was_called = True
                
                # Output the step
                output_process.PrintOutput()

        # Let regular processes do their tasks after outputting the step
        if execute_was_called:
            for process in self._list_of_processes:
                process.ExecuteAfterOutputStep()
            
    
    def _ResetExternalVariables(self):

        # Set to zero variables that might not be overwritten each time step
        zero_vector = np.zeros(self.system_size)
        self._SetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT, zero_vector, broadcast=False)
        self._SetCompleteVector("rigid_body", KMC.IMPOSED_FORCE, KMC.IMPOSED_MOMENT, zero_vector, broadcast=False)
        self._SetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION, zero_vector, broadcast=False)
        self._SetCompleteVector("root_point", KMC.IMPOSED_DISPLACEMENT, KMC.IMPOSED_ROTATION, zero_vector, broadcast=False)


    def _CalculateEffectiveLoad(self):

        # Calculate the total load
        self_weight = self._GetCompleteVector("rigid_body", KM.BODY_FORCE, KM.BODY_MOMENT)
        external_load = self._GetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT)
        imposed_load = self._GetCompleteVector("rigid_body", KMC.IMPOSED_FORCE, KMC.IMPOSED_MOMENT)
        self.total_load = external_load + imposed_load + self_weight

        # Calculate the total root point displacement and the equivalent force it generates
        external_root_point_displ = self._GetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION)
        imposed_root_point_displ = self._GetCompleteVector("root_point", KMC.IMPOSED_DISPLACEMENT, KMC.IMPOSED_ROTATION)
        self.total_root_point_displ = external_root_point_displ + imposed_root_point_displ
        self._UpdateKinematics("root_point", self.total_root_point_displ)
        root_point_force = self._CalculateEquivalentForceFromRootPointDisplacement()
        
        # Sum up both loads and save the resulting effective load
        effective_load = self.total_load + root_point_force
        self._SetCompleteVector("rigid_body", KMC.EFFECTIVE_FORCE, KMC.EFFECTIVE_MOMENT, effective_load)

    
    def _GetKinematics(self, model_part_name, buffer=0):

        # Obtain displacements, velocities and accelerations of all the DOFs
        x = self._GetCompleteVector(model_part_name, KM.DISPLACEMENT, KM.ROTATION, buffer=buffer)
        v = self._GetCompleteVector(model_part_name, KM.VELOCITY, KM.ANGULAR_VELOCITY, buffer=buffer)
        a = self._GetCompleteVector(model_part_name, KM.ACCELERATION, KM.ANGULAR_ACCELERATION, buffer=buffer)

        return x, v, a


    def _UpdateKinematics(self, model_part_name, x):

        # The previous' step kinematics are necessary to calculate the new velocity and acceleration
        x_prev, v_prev, a_prev = self._GetKinematics(model_part_name, buffer=1)

        # Calculate the velocity and acceleration according to the gen-alpha method
        v = self.a1v * (x - x_prev) + self.a2v * v_prev + self.a3v * a_prev
        a = self.a1a * (x - x_prev) + self.a2a * v_prev + self.a3a * a_prev

        # Update the newly calculated kinematics
        self._SetCompleteVector(model_part_name, KM.DISPLACEMENT, KM.ROTATION, x)
        self._SetCompleteVector(model_part_name, KM.VELOCITY, KM.ANGULAR_VELOCITY, v)
        self._SetCompleteVector(model_part_name, KM.ACCELERATION, KM.ANGULAR_ACCELERATION, a)


    def _CalculateReaction(self, buffer=0):

        # Get displacement, velocity and acceleration of both the rigid body and the root point
        x, v, a = self._GetKinematics("rigid_body", buffer=buffer)
        #x_root, v_root, a_root = self._GetKinematics("root_point", buffer=buffer)

        # Calculate the reaction
        # It is important not to calculate from the stiffness,
        # since constrained dofs will have an inf*O factor there.
        reaction = -self.M.dot(a) + self.total_load
        #reaction = self.C.dot(v - v_root) + self.K.dot(x - x_root)

        # Inactive dofs should have a zero reaction, in case they were coupled by mistake
        for index, dof in enumerate(self.available_dofs):
            if dof not in self.active_dofs:
                reaction[index] = 0

        return reaction


    def _CalculateEquivalentForceFromRootPointDisplacement(self):

        # Transform the movement of the root point into an equivalent force
        x_root, v_root, a_root = self._GetKinematics("root_point")
        equivalent_force = self.K.dot(x_root) + self.C.dot(v_root)

        return equivalent_force
    

    def _GetCompleteVector(self, model_part_name, linear_variable, angular_variable, buffer=0):
        # Method to transform the linear and angular variables stored in the model part to a single vector.
        # For example, the displacement vector will be the combination of KM.DISPLACEMENT and KM.ROTATION.
        
        # Obtain the linear and angular data that are going to form the output
        if model_part_name == "rigid_body":
            linear_values = self.rigid_body_model_part.Nodes[1].GetSolutionStepValue(linear_variable, buffer)
            angular_values = self.rigid_body_model_part.Nodes[1].GetSolutionStepValue(angular_variable, buffer)
        elif model_part_name == "root_point":
            linear_values = self.root_point_model_part.Nodes[2].GetSolutionStepValue(linear_variable, buffer)
            angular_values = self.root_point_model_part.Nodes[2].GetSolutionStepValue(angular_variable, buffer)
        else:
            raise Exception('model_part_name should be "rigid_body" or "root_point".')

        # Join all the data in a single vector so the solver can handle it more easily
        return np.array(list(linear_values) + list(angular_values))

    
    def _SetCompleteVector(self, model_part_name, linear_variable, angular_variable, values, buffer=0, broadcast=True):
        # Method to transform a single vector into the linear and angular variables stored in the model part.
        # For example, the displacement vector will splitted and saved in KM.DISPLACEMENT and KM.ROTATION.

        # Decompose the original vector into its linear and angular components
        linear_values = list(values[:self.linear_size])
        angular_values = list(values[-self.angular_size:])

        # Save both of them separately
        if model_part_name == "rigid_body":
            self.rigid_body_model_part.Nodes[1].SetSolutionStepValue(linear_variable, buffer, linear_values)
            self.rigid_body_model_part.Nodes[1].SetSolutionStepValue(angular_variable, buffer, angular_values)
        elif model_part_name == "root_point":
            self.root_point_model_part.Nodes[2].SetSolutionStepValue(linear_variable, buffer, linear_values)
            self.root_point_model_part.Nodes[2].SetSolutionStepValue(angular_variable, buffer, angular_values)
        else:
            raise Exception('model_part_name should be "rigid_body" or "root_point".')


