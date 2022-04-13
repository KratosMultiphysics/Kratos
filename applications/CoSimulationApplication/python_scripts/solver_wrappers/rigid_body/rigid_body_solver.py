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


class RigidBodySolver(object):
    """
    This class implements a Rigid Body solver independent of Kratos.
    Several types of load applications are available, and they can be applyed to each degree of freedom.
    """

    def __init__(self, model, project_parameters):

        if not isinstance(model, KM.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")
        if not isinstance(project_parameters, KM.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        input_check._CheckMandatoryInputParameters(project_parameters)
        self.problem_name = project_parameters["problem_data"]["problem_name"].GetString()
        self.start_time = project_parameters["problem_data"]["start_time"].GetDouble()
        self.end_time = project_parameters["problem_data"]["end_time"].GetDouble()

        solver_settings = input_check._ValidateAndAssignRigidBodySolverDefaults(project_parameters["solver_settings"])
        self.domain_size = solver_settings["domain_size"].GetInt()
        buffer_size = solver_settings["buffer_size"].GetInt()
        self.echo_level = solver_settings["echo_level"].GetInt()

        self.model = model
        self.main_model_part = self.model.CreateModelPart("Main")
        if solver_settings["model_import_settings"]["input_type"].GetString() == "none":
            self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.domain_size
            self.rigid_body_model_part = self.main_model_part.CreateSubModelPart("RigidBody")
            self.root_point_model_part = self.main_model_part.CreateSubModelPart("RootPoint")
            self.AddVariables()
            self.rigid_body_model_part.CreateNewNode(1,0.0,0.0,0.0)
            self.root_point_model_part.CreateNewNode(2,0.0,0.0,0.0)
            self.main_model_part.SetBufferSize(buffer_size)
        elif solver_settings["model_import_settings"]["input_type"].GetString() == "rest":
            model_import_settings = solver_settings["model_import_settings"]
            model_import_settings.RemoveValue("input_type")
            RestartUtility(self.main_model_part, model_import_settings).LoadRestart()
            self.rigid_body_model_part = self.main_model_part.GetSubModelPart("RigidBody")
            self.root_point_model_part = self.main_model_part.GetSubModelPart("RootPoint")
            self.start_time = self.main_model_part.ProcessInfo[KM.TIME]

        # Degrees of freedom that can be activated from the parameters
        if self.domain_size == 3:
            self.available_dofs = ['displacement_x', 'displacement_y', 'displacement_z',
                'rotation_x', 'rotation_y', 'rotation_z']
        elif self.domain_size == 2:
            self.available_dofs = ['displacement_x', 'displacement_y', 'rotation_z']

        self.system_size = len(self.available_dofs)
        self.linear_size = int(np.ceil(self.system_size/2))
        self.angular_size = int(np.floor(self.system_size/2))

        # Fill with defaults and check that mandatory fields are given
        dof_settings, self.active_dofs = input_check._ValidateAndAssignDofDefaults(solver_settings["active_dofs"], self.available_dofs)
        
        # Create all the processes stated in the project parameters
        self._list_of_processes = input_check._CreateListOfProcesses(self.model, project_parameters, self.main_model_part)
        self._list_of_output_processes = input_check._CreateListOfOutputProcesses(self.model, project_parameters)

        # Safe all the filled data in their respective class variables
        self._InitializeSolutionVariables(solver_settings)
        self._InitializeDofsVariables(dof_settings)

        # Prepare the parameters for the generalized-alpha method
        self._InitializeGeneralizedAlphaParameters()
        

    def AddVariables(self):

        self.main_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ANGULAR_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ANGULAR_ACCELERATION)

        self.rigid_body_model_part.AddNodalSolutionStepVariable(KM.FORCE)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KM.MOMENT)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KMC.PRESCRIBED_FORCE)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KMC.PRESCRIBED_MOMENT)

        self.root_point_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.root_point_model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        self.root_point_model_part.AddNodalSolutionStepVariable(KMC.PRESCRIBED_DISPLACEMENT)
        self.root_point_model_part.AddNodalSolutionStepVariable(KMC.PRESCRIBED_ROTATION)


    def _InitializeSolutionVariables(self, solver_settings):
        
        # Save all the data that does not depend on the degree of freedom
        self.rho_inf = solver_settings["time_integration_parameters"]["rho_inf"].GetDouble()
        self.delta_t = solver_settings["time_integration_parameters"]["time_step"].GetDouble()
        self.buffer_size = solver_settings["buffer_size"].GetInt()


    def _InitializeDofsVariables(self, dof_settings):

        # Initialize variables that depend on the degree of freedom
        self.is_constrained = {}
        self.M = np.zeros((self.system_size,self.system_size)) # Mass matrix
        self.C = np.zeros((self.system_size,self.system_size)) # Damping matrix
        self.K = np.zeros((self.system_size,self.system_size)) # Stiffness matrix
        self.modulus_self_weight = np.zeros(self.system_size) # Gravity acceleration

        # Fill the initialised values with the ones from the parameters
        for index, dof in enumerate(self.available_dofs):
            self.is_constrained[dof] = dof_settings[dof]["constrained"].GetBool()
            self.M[index][index] = dof_settings[dof]['system_parameters']['mass'].GetDouble()
            self.C[index][index] = dof_settings[dof]['system_parameters']['damping'].GetDouble()
            self.K[index][index] = dof_settings[dof]['system_parameters']['stiffness'].GetDouble()
            self.modulus_self_weight[index] = dof_settings[dof]['system_parameters']['modulus_self_weight'].GetDouble()


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


    def Run(self):

        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()


    def Initialize(self):

        # Initialize the time
        if not self.main_model_part.ProcessInfo[KM.IS_RESTARTED]:
            self.main_model_part.ProcessInfo[KM.TIME] = self.start_time

        # Other variables that are used in SolveSolutionStep()
        self.total_root_point_displ = np.zeros(self.system_size)
        self.total_load = np.zeros(self.system_size)
        self.effective_load = np.zeros((self.system_size, self.buffer_size))

        for process in self._list_of_processes:
            process.ExecuteInitialize()

        for process in self._list_of_processes:
            process.ExecuteBeforeSolutionLoop()


    def RunSolutionLoop(self):

        while self.main_model_part.ProcessInfo[KM.TIME] < self.end_time:
            self.AdvanceInTime(self.main_model_part.ProcessInfo[KM.TIME])
            self.InitializeSolutionStep()
            self.Predict()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    
    def Finalize(self):

        for process in self._list_of_processes:
            process.ExecuteFinalize()


    def AdvanceInTime(self, current_time):    
        # Similar to the Kratos CloneTimeStep function
        # Advances values along the buffer axis (so rolling columns) using numpy's roll
        # Column 0 is the current time step, column 1 the previous one...

        # Variables whith buffer. Column 0 will be overwriten later
        #self.total_root_point_displ = np.roll(self.total_root_point_displ,1,axis=1)
        self.effective_load = np.roll(self.effective_load,1,axis=1)

        # Variables that need to be reseted. They might not be overwriten later so they
        # need to be zero to avoid values from previous time steps ar not continously used.
        self.total_load = np.zeros(self.system_size)

        # Update the time of the simulation
        time = current_time + self.delta_t
        self.main_model_part.CloneTimeStep(time)
        self.main_model_part.ProcessInfo[KM.STEP] += 1
        self._ResetExternalVariables()

        return time

    
    def InitializeSolutionStep(self):
        
        for process in self._list_of_processes:
            process.ExecuteInitializeSolutionStep()

    
    def Predict(self):
        pass


    def SolveSolutionStep(self):
        
        # Calculate the effective load which considers both actual loads
        # and the equivalent load from the root point displacement
        self.effective_load[:,0] = self._CalculateEffectiveLoad()

        # Calculate the gen-alpha load for the construction of the RHS
        F = (1.0 - self.alpha_f) * self.effective_load[:,0] + self.alpha_f * self.effective_load[:,1]

        x_prev, v_prev, a_prev = self._GetKinematics("rigid_body", buffer=1)

        # Creation of the RHS according to the gen-alpha method
        RHS = np.dot(self.M, (self.a1m * x_prev + self.a2m * v_prev + self.a3m * a_prev))
        RHS += np.dot(self.C, (self.a1b * x_prev + self.a2b * v_prev + self.a3b * a_prev))
        RHS += np.dot(self.a1k * self.K, x_prev) + F

        # Make zero the corrresponding values of RHS so the inactive dofs are not excited
        for index, dof in enumerate(self.available_dofs):
            if dof not in self.active_dofs:
                RHS[index] = 0

        # Solve the solution step and find the new displacements
        x = np.linalg.solve(self.LHS, RHS)

        # constrained dofs will have the root_point_displacement as a total displacement
        for index, dof in enumerate(self.available_dofs):
            if self.is_constrained[dof]:
                x[index] = self.total_root_point_displ[index]

        # Update velocity and acceleration according to the gen-alpha method
        self._UpdateDisplacement("rigid_body", x)

        reaction = self._CalculateReaction()
        self._SetCompleteVector("root_point", KM.REACTION, KM.REACTION_MOMENT, reaction)

    
    def FinalizeSolutionStep(self):
        
        for process in self._list_of_processes:
            process.ExecuteFinalizeSolutionStep()


    def OutputSolutionStep(self):
        """This function writes output files after the solution of a step, exactly as in AnalysisStage()
        """
        # Carry out this task only in the first rank (for parallel computing)
        # TODO: This might not be necessary anymore, since it doesn't run in MPI
        data_comm = KM.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:

            execute_was_called = False
            for output_process in self._list_of_output_processes:
                if output_process.IsOutputStep():
                    if not execute_was_called:
                        for process in self._list_of_processes:
                            process.ExecuteBeforeOutputStep()
                        execute_was_called = True

                    output_process.PrintOutput()

            if execute_was_called:
                for process in self._list_of_processes:
                    process.ExecuteAfterOutputStep()
            
    
    def _ResetExternalVariables(self):
        zero_vector = np.zeros(self.system_size)
        self._SetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT, zero_vector)
        self._SetCompleteVector("rigid_body", KMC.PRESCRIBED_FORCE, KMC.PRESCRIBED_MOMENT, zero_vector)
        self._SetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION, zero_vector)
        self._SetCompleteVector("root_point", KMC.PRESCRIBED_DISPLACEMENT, KMC.PRESCRIBED_ROTATION, zero_vector)


    def _CalculateEffectiveLoad(self):

        # Calculate the total load
        self_weight = self._CalculateSelfWeight()
        external_load = self._GetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT)
        prescribed_load = self._GetCompleteVector("rigid_body", KMC.PRESCRIBED_FORCE, KMC.PRESCRIBED_MOMENT)
        self.total_load = external_load + prescribed_load + self_weight

        # Calculate the total root point displacement and the equivalent force it generates
        external_root_point_displ = self._GetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION)
        prescribed_root_point_displ = self._GetCompleteVector("root_point", KMC.PRESCRIBED_DISPLACEMENT, KMC.PRESCRIBED_ROTATION)
        self.total_root_point_displ = external_root_point_displ + prescribed_root_point_displ
        self._UpdateDisplacement("root_point", self.total_root_point_displ)
        root_point_force = self._CalculateEquivalentForceFromRootPointDisplacement()
        
        # Sum up both loads
        effective_load = self.total_load + root_point_force

        return effective_load

    
    def _GetKinematics(self, model_part_name, buffer=0):

        x = self._GetCompleteVector(model_part_name, KM.DISPLACEMENT, KM.ROTATION, buffer=buffer)
        v = self._GetCompleteVector(model_part_name, KM.VELOCITY, KM.ANGULAR_VELOCITY, buffer=buffer)
        a = self._GetCompleteVector(model_part_name, KM.ACCELERATION, KM.ANGULAR_ACCELERATION, buffer=buffer)

        return x, v, a


    def _UpdateDisplacement(self, model_part_name, x):

        x_prev, v_prev, a_prev = self._GetKinematics(model_part_name, buffer=1)

        # Calculate the velocity and acceleration according to the gen-alpha method
        v = self.a1v * (x - x_prev) + self.a2v * v_prev + self.a3v * a_prev
        a = self.a1a * (x - x_prev) + self.a2a * v_prev + self.a3a * a_prev

        self._SetCompleteVector(model_part_name, KM.DISPLACEMENT, KM.ROTATION, x)
        self._SetCompleteVector(model_part_name, KM.VELOCITY, KM.ANGULAR_VELOCITY, v)
        self._SetCompleteVector(model_part_name, KM.ACCELERATION, KM.ANGULAR_ACCELERATION, a)


    def _CalculateReaction(self, buffer=0):

        x, v, a = self._GetKinematics("rigid_body", buffer=buffer)
        x_root, v_root, a_root = self._GetKinematics("root_point", buffer=buffer)
        # TODO: Check how does it work with constrained DOFs
        # This is a very big TODO!
        reaction = self.C.dot(v - v_root) + self.K.dot(x - x_root)

        return reaction


    def _CalculateSelfWeight(self):

        self_weight = self.M.dot(self.modulus_self_weight)
        
        return self_weight


    def _CalculateEquivalentForceFromRootPointDisplacement(self):

        # Transform the movement of the root point into an equivalent force
        x_root, v_root, a_root = self._GetKinematics("root_point")
        equivalent_force = self.K.dot(x_root) + self.C.dot(v_root)

        return equivalent_force
    

    def _GetCompleteVector(self, model_part_name, linear_variable, angular_variable, buffer=0):

        if model_part_name == "rigid_body":
            linear_values = self.rigid_body_model_part.Nodes[1].GetSolutionStepValue(linear_variable, buffer)
            angular_values = self.rigid_body_model_part.Nodes[1].GetSolutionStepValue(angular_variable, buffer)
        elif model_part_name == "root_point":
            linear_values = self.root_point_model_part.Nodes[2].GetSolutionStepValue(linear_variable, buffer)
            angular_values = self.root_point_model_part.Nodes[2].GetSolutionStepValue(angular_variable, buffer)
        else:
            raise Exception('model_part_name should be "rigid_body" or "root_point".')

        return np.array(list(linear_values) + list(angular_values))

    
    def _SetCompleteVector(self, model_part_name, linear_variable, angular_variable, values, buffer=0):

        linear_values = list(values[:self.linear_size])
        angular_values = list(values[-self.angular_size:])

        if model_part_name == "rigid_body":
            self.rigid_body_model_part.Nodes[1].SetSolutionStepValue(linear_variable, buffer, linear_values)
            self.rigid_body_model_part.Nodes[1].SetSolutionStepValue(angular_variable, buffer, angular_values)
        elif model_part_name == "root_point":
            self.root_point_model_part.Nodes[2].SetSolutionStepValue(linear_variable, buffer, linear_values)
            self.root_point_model_part.Nodes[2].SetSolutionStepValue(angular_variable, buffer, angular_values)
        else:
            raise Exception('model_part_name should be "rigid_body" or "root_point".')


