# Kratos imports
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# RigidBody imports
from . import rigid_body_process
from . import input_check

# Other imports
import numpy as np
import json
import os
from importlib import import_module

class RigidBodySolver(object):
    """
    This class implements a Rigid Body solver independent of Kratos.
    Several types of load applications are available, and they can be applyed to each degree of freedom.
    """

    def __init__(self, model, parameters):

        # TODO: Centralized input check
        #self.parameters = input_check.FullInputCheck(parameters)

        self.model = model
        self.main_model_part = self.model.CreateModelPart("Main")
        self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        self.rigid_body_model_part = self.main_model_part.CreateSubModelPart("RigidBody")
        self.root_point_model_part = self.main_model_part.CreateSubModelPart("RootPoint")

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
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KMC.PRESCRIBED_DISPLACEMENT)
        self.rigid_body_model_part.AddNodalSolutionStepVariable(KMC.PRESCRIBED_ROTATION)

        self.rigid_body_model_part.CreateNewNode(1,0.0,0.0,0.0)
        self.root_point_model_part.CreateNewNode(2,0.0,0.0,0.0)

        # Degrees of freedom that can be activated from the parameters
        self.available_dofs = ['displacement_x', 'displacement_y', 'displacement_z',
            'rotation_x', 'rotation_y', 'rotation_z']
        self.system_size = len(self.available_dofs)

        # How many dofs are linear (displacement/force) and how many are angular (rotation/moment)?
        # TODO: For future implementation of a 2D version. It will only be needed to change self.available_dofs
        # CAN BE ERASED
        self.linear_size = int(np.ceil(self.system_size/2))
        self.angular_size = int(np.floor(self.system_size/2))

        # Note which variable labels are linear/angular for later inpu/output use
        # CAN BE ERASED
        self.linear_variables = ["DISPLACEMENT", "ROOT_POINT_DISPLACEMENT", "FORCE", "REACTION"]
        self.angular_variables = ["ROTATION", "ROOT_POINT_ROTATION", "MOMENT", "REACTION_MOMENT"]

        # Check that the activated dofs are not repeated and are among the available ones
        self.active_dofs = input_check._CheckActiveDofs(parameters, self.available_dofs)

        # Fill with defaults and check that mandatory fields are given
        input_check._CheckMandatoryInputParameters(parameters)
        dof_params, sol_params = input_check._ValidateAndAssignRigidBodySolverDefaults(parameters, self.available_dofs)

        # Create all the processes stated in the project parameters
        parameters_json = json.loads(parameters.WriteJsonString())
        if "processes" in parameters_json:
            self.process_list = []
            for process in parameters_json["processes"]:
                python_module = process["python_module"]
                kratos_module = process["kratos_module"]
                process_module = import_module(kratos_module + "." + python_module)
                process_settings = KM.Parameters(json.dumps(process["Parameters"]))
                self.process_list.append(process_module.Factory(self, process_settings))
        
        # Safe all the filled data in their respective class variables
        self._InitializeDofsVariables(dof_params)
        self._InitializeSolutionVariables(sol_params)

        # Prepare the parameters for the generalized-alpha method
        self._InitializeGeneralizedAlphaParameters()

        # Set buffer size
        self.main_model_part.SetBufferSize(self.buffer_size)


    def _InitializeSolutionVariables(self, sol_params):
        
        # Save all the data that does not depend on the degree of freedom
        self.rho_inf = sol_params["time_integration_parameters"]["rho_inf"].GetDouble()
        self.delta_t = sol_params["time_integration_parameters"]["time_step"].GetDouble()
        self.start_time = sol_params["time_integration_parameters"]["start_time"].GetDouble()
        self.buffer_size = sol_params["solver_parameters"]["buffer_size"].GetInt()
        self.output_file_path = sol_params["output_parameters"]["file_path"].GetString()
        self.write_output_file = sol_params["output_parameters"]["write_output_files"].GetBool()


    def _InitializeDofsVariables(self, dof_params):

        # Initialize variables that depend on the degree of freedom
        self.is_constrained = {}
        self.M = np.zeros((self.system_size,self.system_size)) # Mass matrix
        self.C = np.zeros((self.system_size,self.system_size)) # Damping matrix
        self.K = np.zeros((self.system_size,self.system_size)) # Stiffness matrix
        self.modulus_self_weight = np.zeros(self.system_size) # Gravity acceleration
        self.initial_displacement = np.zeros(self.system_size)
        self.initial_velocity = np.zeros(self.system_size)
        self.initial_acceleration = np.zeros(self.system_size)
        self.load_impulse = np.zeros(self.system_size)

        # Fill the initialised values with the ones from the parameters
        for index, dof in enumerate(self.available_dofs):
            self.is_constrained[dof] = dof_params[dof]["constrained"].GetBool()
            self.M[index][index] = dof_params[dof]['system_parameters']['mass'].GetDouble()
            self.C[index][index] = dof_params[dof]['system_parameters']['damping'].GetDouble()
            self.K[index][index] = dof_params[dof]['system_parameters']['stiffness'].GetDouble()
            self.modulus_self_weight[index] = dof_params[dof]['system_parameters']['modulus_self_weight'].GetDouble()
            self.initial_displacement[index] = dof_params[dof]["initial_conditions"]["displacement"].GetDouble()
            self.initial_velocity[index] = dof_params[dof]["initial_conditions"]["velocity"].GetDouble()
            self.load_impulse[index] = dof_params[dof]["initial_conditions"]["load_impulse"].GetDouble()
        
        # Calculate the initial acceleration from the initial conditions
        factor = self.load_impulse - self.K.dot(self.initial_displacement)
        self.initial_acceleration = np.dot(np.linalg.inv(self.M), factor)


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
        

    def Initialize(self):

        # Initialize the time
        self.time = self.start_time

        # Initialize total displacement, velocity and acceleration
        #self.x = np.zeros((self.system_size, self.buffer_size))
        #self.v = np.zeros((self.system_size, self.buffer_size))
        #self.a = np.zeros((self.system_size, self.buffer_size))

        # Initialize the displacement, velocity and accelerration of the root point
        #self.x_root = np.zeros((self.system_size, self.buffer_size))
        #self.v_root = np.zeros((self.system_size, self.buffer_size))
        #self.a_root = np.zeros((self.system_size, self.buffer_size))

        # Other variables that are used in SolveSolutionStep()
        self.total_root_point_displ = np.zeros((self.system_size, self.buffer_size))
        self.total_load = np.zeros(self.system_size)
        #self.prescribed_load = np.zeros(self.system_size)
        #self.prescribed_root_point_displ = np.zeros(self.system_size)
        #self.external_load = np.zeros(self.system_size)
        #self.external_root_point_displ = np.zeros(self.system_size)
        self.effective_load = np.zeros((self.system_size, self.buffer_size))

        # Apply initial conditions
        initial_disp = list(self.initial_displacement[:self.linear_size])
        initial_rot = list(self.initial_displacement[-self.angular_size:])
        initial_vel = list(self.initial_velocity[:self.linear_size])
        initial_ang_vel = list(self.initial_velocity[-self.angular_size:])
        initial_acc = list(self.initial_acceleration[:self.linear_size])
        initial_ang_acc = list(self.initial_acceleration[-self.angular_size:])
        self._SetCompleteVector("rigid_body", KM.DISPLACEMENT, KM.ROTATION, initial_disp+initial_rot)
        self._SetCompleteVector("rigid_body", KM.VELOCITY, KM.ANGULAR_VELOCITY, initial_disp+initial_rot)
        self._SetCompleteVector("rigid_body", KM.ACCELERATION, KM.ANGULAR_ACCELERATION, initial_disp+initial_rot)
        #self.x[:,0] = self.initial_displacement
        #self.v[:,0] = self.initial_velocity
        #self.a[:,0] = self.initial_acceleration

        # Apply external load as an initial impulse
        self.total_load = self.load_impulse
        self.effective_load[:,0] = self.total_load

        # Create output file and wrrite the time=start_time step
        if self.write_output_file:
            self.InitializeOutput()


    def InitializeOutput(self):

        # Carry out this task only in the first rank (for parallel computing)
        data_comm = KM.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:

            # Create the directory where the results will be stored
            if not os.path.exists(self.output_file_path):
                os.makedirs(self.output_file_path)
            
            # File paths for each DOF will be stored here
            self.output_file_name = {}

            # Loop through all the possible DOFs
            for dof in self.available_dofs:

                # Save output file path for this dof
                dof_file_name = os.path.join(self.output_file_path, dof + '.dat')

                # If the DOF is active, create the results file and write the headers
                if dof in self.active_dofs:
                    self.output_file_name[dof] = dof_file_name
                    with open(self.output_file_name[dof], "w") as results_rigid_body:
                        results_rigid_body.write("time"+ " " +
                                                "displacement" + " " +
                                                "velocity" + " " +
                                                "acceleration" + " " +
                                                "root_point_displacement" + " " +
                                                "root_point_velocity" + " " +
                                                "root_point_acceleration" + " " +
                                                "relative_displacement" + " " +
                                                "relative_velocity" + " " +
                                                "relative_accleration" + " " +
                                                "reaction" + "\n")
                
                # If the DOF is not active, erase results files from other runs
                else:
                    if os.path.isfile(dof_file_name):
                        os.remove(dof_file_name)
            
            # Write initial time step output
            self.OutputSolutionStep()


    def OutputSolutionStep(self):

        # Carry out this task only in the first rank (for parallel computing)
        data_comm = KM.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:

            # Only write if the output is enabled
            if self.write_output_file:

                # Calculate the reaction for the output
                # TODO: this sould be calculated elsewhere or not outputed
                reaction = self.CalculateReaction()

                x, v, a = self._GetKinematics("rigid_body")
                x_root, v_root, a_root = self._GetKinematics("root_point")

                # Write the output only for the active DOFs
                for index, dof in enumerate(self.available_dofs):
                    if dof in self.active_dofs:
                        with open(self.output_file_name[dof], "a") as results_rigid_body:
                            results_rigid_body.write(str(self.time) + " " +
                                                    str(x[index]) + " " +
                                                    str(v[index]) + " " +
                                                    str(a[index]) + " " +
                                                    str(x_root[index]) + " " +
                                                    str(v_root[index]) + " " +
                                                    str(a_root[index]) + " " +
                                                    str(x[index] - x_root[index]) + " " +
                                                    str(v[index] - v_root[index]) + " " +
                                                    str(a[index] - a_root[index]) + " " +
                                                    str(reaction[index]) + "\n")


    def AdvanceInTime(self, current_time):
        # Similar to the Kratos CloneTimeStep function
        # Advances values along the buffer axis (so rolling columns) using numpy's roll
        # Column 0 is the current time step, column 1 the previous one...

        # Variables whith buffer. Column 0 will be overwriten later
        #self.x = np.roll(self.x,1,axis=1)
        #self.x_root = np.roll(self.x_root,1,axis=1)
        #self.v = np.roll(self.v,1,axis=1)
        #self.v_root = np.roll(self.v_root,1,axis=1)
        #self.a = np.roll(self.a,1,axis=1)
        #self.a_root = np.roll(self.a_root,1,axis=1)
        self.total_root_point_displ = np.roll(self.total_root_point_displ,1,axis=1)
        self.effective_load = np.roll(self.effective_load,1,axis=1)

        # Variables that need to be reseted. They might not be overwriten later so they
        # need to be zero to avoid values from previous time steps ar not continously used.
        self.total_load = np.zeros(self.system_size)
        #self.prescribed_load = np.zeros(self.system_size)
        #self.prescribed_root_point_displ = np.zeros(self.system_size)
        #self.external_load = np.zeros(self.system_size)
        #self.external_root_point_displ = np.zeros(self.system_size)

        # Update the time of the simulation
        self.time = current_time + self.delta_t
        self.main_model_part.CloneTimeStep(self.time)
        self.main_model_part.ProcessInfo[KM.STEP] += 1
        self._ResetExternalVariables()

        return self.time


    def SolveSolutionStep(self):
        
        # Execute all generated processes
        # TODO: Call the corresponding methods in other places
        for process in self.process_list:
            process.ExecuteBeforeSolutionLoop()
        
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
                x[index] = self.total_root_point_displ[index,0]

        # Update velocity and acceleration according to the gen-alpha method
        self._UpdateDisplacement("rigid_body", x)
        #self.v[:,0] = self.a1v * (x - x_prev) + self.a2v * \
        #                self.v[:,1] + self.a3v * self.a[:,1]
        #self.a[:,0] = self.a1a * (x - x_prev) + self.a2a * \
        #                self.v[:,1] + self.a3a * self.a[:,1]
         
        #self._SetCompleteVector("rigid_body", KM.DISPLACEMENT, KM.ROTATION, x)


    def _UpdateRootPointDisplacement(self, displ):
        # Save root point displacement
        self.x_root[:,0] = displ
        # Update the velocity and acceleration according to the gen-alpha method
        self.v_root[:,0] = self.a1v * (self.x_root[:,0] - self.x_root[:,1]) + self.a2v * \
            self.v_root[:,1] + self.a3v * self.a_root[:,1]
        self.a_root[:,0] = self.a1a * (self.x_root[:,0] - self.x_root[:,1]) + self.a2a * \
            self.v_root[:,1] + self.a3a * self.a_root[:,1]

    def _UpdateDisplacement(self, model_part_name, x):

        x_prev, v_prev, a_prev = self._GetKinematics(model_part_name, buffer=1)

        # Calculate the velocity and acceleration according to the gen-alpha method
        v = self.a1v * (x - x_prev) + self.a2v * v_prev + self.a3v * a_prev
        a = self.a1a * (x - x_prev) + self.a2a * v_prev + self.a3a * a_prev

        self._SetCompleteVector(model_part_name, KM.DISPLACEMENT, KM.ROTATION, x)
        self._SetCompleteVector(model_part_name, KM.VELOCITY, KM.ANGULAR_VELOCITY, v)
        self._SetCompleteVector(model_part_name, KM.ACCELERATION, KM.ANGULAR_ACCELERATION, a)

    
    def _GetKinematics(self, model_part_name, buffer=0):

        x = self._GetCompleteVector(model_part_name, KM.DISPLACEMENT, KM.ROTATION, buffer=buffer)
        v = self._GetCompleteVector(model_part_name, KM.VELOCITY, KM.ANGULAR_VELOCITY, buffer=buffer)
        a = self._GetCompleteVector(model_part_name, KM.ACCELERATION, KM.ANGULAR_ACCELERATION, buffer=buffer)

        return x, v, a


    def _CalculateEquivalentForceFromRootPointDisplacement(self):
        # Transform the movement of the root point into an equivalent force
        x_root, v_root, a_root = self._GetKinematics("root_point")
        equivalent_force = self.K.dot(x_root) + self.C.dot(v_root)
        return equivalent_force


    def _CalculateEffectiveLoad(self):

        # Calculate the total load
        self_weight = self.CalculateSelfWeight()
        external_load = self._GetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT)
        prescribed_load = self._GetCompleteVector("rigid_body", KMC.PRESCRIBED_FORCE, KMC.PRESCRIBED_MOMENT)
        self.total_load = external_load + prescribed_load + self_weight

        # Calculate the total root point displacement and the equivalent force it generates
        external_root_point_displ = self._GetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION)
        prescribed_root_point_displ = self._GetCompleteVector("root_point", KMC.PRESCRIBED_DISPLACEMENT, KMC.PRESCRIBED_ROTATION)
        self.total_root_point_displ[:,0] = external_root_point_displ + prescribed_root_point_displ
        self._UpdateDisplacement("root_point", self.total_root_point_displ[:,0])
        root_point_force = self._CalculateEquivalentForceFromRootPointDisplacement()
        
        # Sum up both loads
        effective_load = self.total_load + root_point_force

        return effective_load
    
    '''
    def _UpdateVelocityAndAcceleration(self, x):
        # Update the velocity and acceleration according to the gen-alpha method
        self.v[:,0] = self.a1v * (x - self.x[:,1]) + self.a2v * \
                        self.v[:,1] + self.a3v * self.a[:,1]
        self.a[:,0] = self.a1a * (x - self.x[:,1]) + self.a2a * \
                        self.v[:,1] + self.a3a * self.a[:,1]
    '''

    def CalculateReaction(self, buffer=0):
        x, v, a = self._GetKinematics("rigid_body", buffer=buffer)
        x_root, v_root, a_root = self._GetKinematics("root_point", buffer=buffer)
        # TODO: Check how does it work with constrained DOFs
        reaction = self.C.dot(v - v_root) + self.K.dot(x - x_root)
        return reaction


    def CalculateSelfWeight(self):
        self_weight = self.M.dot(self.modulus_self_weight)
        return self_weight


    # CAN BE ERASED
    def SetSolutionStepValue(self, identifier, values, buffer_idx=0):

        # Check that the selected buffer value is OK
        input_check._CheckBufferId(buffer_idx,identifier)
        
        # Check that the input's size is the expected one
        expected_size = self._ExpectedDataSize(identifier)
        if len(values) != expected_size:
            msg = 'The variable "' + identifier + '" does not have the '
            msg += 'right size. It has size ' + str(len(values))
            msg += ' but it should be of size ' + str(expected_size)
            raise Exception(msg)

        x = np.zeros(expected_size)

        # Loop through input active DOFs saving the values
        for index, value in enumerate(values):
            # Increase the index so it fits with the angular terms (if necessary)
            if identifier in self.angular_variables:
                index += self.linear_size
            # Save input variables in their corresponding spots
            if self.available_dofs[index] in self.active_dofs:
                if identifier in ["DISPLACEMENT", "ROTATION", "DISPLACEMENTS_ALL"]:
                    x[index] = value
                elif identifier == "VELOCITY":
                    self.v[index, buffer_idx] = value
                elif identifier == "ACCELERATION":
                    self.a[index, buffer_idx] = value
                elif identifier in ["FORCE", "MOMENT", "FORCE_ALL"]:
                    self.external_load[index] = value
                elif identifier in ["ROOT_POINT_DISPLACEMENT", "ROOT_POINT_ROTATION", "ROOT_POINT_DISPLACEMENT_ALL"]:
                    self.external_root_point_displ[index] = value
                else:
                    raise Exception("Identifier is unknown!")
        
        #self._SetCompleteVector("rigid_body", KM.DISPLACEMENT, KM.ROTATION, x)


    # CAN BE ERASED
    def GetSolutionStepValue(self, identifier, buffer_idx=0):

        # Check that the selected buffer value is OK
        input_check._CheckBufferId(buffer_idx, identifier)

        # Initialize output as a KM vector
        output = KM.Vector(self._ExpectedDataSize(identifier))

        x = self._GetCompleteVector("rigid_body", KM.DISPLACEMENT, KM.ROTATION)

        # Save all the values from the active DOFs
        for index in range(len(output)):
            out_index = index
            # Increase the index so it fits with the angular terms (if necessary)
            if identifier in self.angular_variables:
                index += self.linear_size
            # Fill the output with its corresponding values
            if identifier in ["DISPLACEMENT", "ROTATION", "DISPLACEMENTS_ALL"]:
                output[out_index] = x[index]
            elif identifier == "VELOCITY":
                output[out_index] = self.v[index, buffer_idx]
            elif identifier == "ACCELERATION":
                output[out_index] = self.a[index, buffer_idx]
            elif identifier in ["REACTION", "REACTION_MOMENT", "REACTION_ALL"]:
                output[out_index] = self.CalculateReaction(buffer=buffer_idx)[index]
            else:
                raise Exception("Identifier is unknown!")
        return output

    
    # CAN BE ERASED
    def _ExpectedDataSize(self, identifier):

        if identifier in self.linear_variables:
            expected_size = self.linear_size
        elif identifier in self.angular_variables:
            expected_size = self.angular_size
        else:
            expected_size = self.system_size

        return expected_size

    
    def _ResetExternalVariables(self):
        zero_vector = np.zeros(self.system_size)
        self._SetCompleteVector("rigid_body", KM.FORCE, KM.MOMENT, zero_vector)
        self._SetCompleteVector("rigid_body", KMC.PRESCRIBED_FORCE, KMC.PRESCRIBED_MOMENT, zero_vector)
        self._SetCompleteVector("root_point", KM.DISPLACEMENT, KM.ROTATION, zero_vector)
        self._SetCompleteVector("root_point", KMC.PRESCRIBED_DISPLACEMENT, KMC.PRESCRIBED_ROTATION, zero_vector)
    
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


