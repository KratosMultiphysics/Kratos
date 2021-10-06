# Kratos imports
import KratosMultiphysics

# RigidBody imports
from . import rigid_body_process
from . import input_check

# Other imports
import numpy as np
import json
import os

class RigidBodySolver(object):
    """
    This class implements a Rigid Body solver independent of Kratos.
    Several types of load applications are available, and they can be applyed to each degree of freedom.
    """

    def __init__(self, input_name):

        # Allow different formats for the input parameters
        if isinstance(input_name, dict):
            parameters = input_name
        elif isinstance(input_name, str):
            if not input_name.endswith(".json"):
                input_name += ".json"
            with open(input_name,'r') as parameter_file:
                parameters = json.load(parameter_file)
        else:
            raise Exception("The input has to be provided as a dict or a string")

        # Degrees of freedom that can be activated from the parameters
        self.available_dofs = ['displacement_x', 'displacement_y', 'displacement_z',
            'rotation_x', 'rotation_y', 'rotation_z']
        self.system_size = len(self.available_dofs)

        # How many dofs are linear (displacement/force) and how many are angular (rotation/moment)?
        # TODO: For future implementation of a 2D version. It will only be needed to change self.available_dofs
        self.linear_size = int(np.ceil(self.system_size/2))
        self.angular_size = int(np.floor(self.system_size/2))

        # Note which variable labels are linear/angular for later inpu/output use
        self.linear_variables = ["DISPLACEMENT", "ROOT_POINT_DISPLACEMENT", "FORCE", "REACTION"]
        self.angular_variables = ["ROTATION", "ROOT_POINT_ROTATION", "MOMENT", "REACTION_MOMENT"]

        # Check that the activated dofs are not repeated and are among the available ones
        self.active_dofs = input_check._CheckActiveDofs(parameters, self.available_dofs)

        # Fill with defaults and check that mandatory fields are given
        input_check._CheckMandatoryInputParameters(parameters)
        dof_params, sol_params = input_check._ValidateAndAssignRigidBodySolverDefaults(parameters, self.available_dofs)

        # Create all the processes stated in the project parameters
        if "processes" in parameters:
            self.process_list = []
            for process_settings in parameters["processes"]:
                process_settings = KratosMultiphysics.Parameters(json.dumps(process_settings))
                self.process_list.append(rigid_body_process.CreateRigidBodyProcess(self, process_settings))
        
        # Safe all the filled data in their respective class variables
        self._InitializeDofsVariables(dof_params)
        self._InitializeSolutionVariables(sol_params)

        # Prepare the parameters for the generalized-alpha method
        self._InitializeGeneralizedAlphaParameters()


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
        self.x = np.zeros((self.system_size, self.buffer_size))
        self.v = np.zeros((self.system_size, self.buffer_size))
        self.a = np.zeros((self.system_size, self.buffer_size))

        # Initialize the displacement, velocity and accelerration of the root point
        self.x_f = np.zeros((self.system_size, self.buffer_size))
        self.v_f = np.zeros((self.system_size, self.buffer_size))
        self.a_f = np.zeros((self.system_size, self.buffer_size))

        # Other variables that are used in SolveSolutionStep()
        self.total_root_point_displ = np.zeros((self.system_size, self.buffer_size))
        self.total_load = np.zeros((self.system_size, self.buffer_size))
        self.prescribed_load = np.zeros(self.system_size)
        self.prescribed_root_point_displ = np.zeros(self.system_size)
        self.external_load = np.zeros(self.system_size)
        self.external_root_point_displ = np.zeros(self.system_size)
        self.effective_load = np.zeros((self.system_size, self.buffer_size))

        # Apply initial conditions
        self.x[:,0] = self.initial_displacement
        self.v[:,0] = self.initial_velocity
        self.a[:,0] = self.initial_acceleration

        # Apply external load as an initial impulse
        self.total_load[:,0] = self.load_impulse
        self.effective_load[:,0] = self.total_load[:,0]

        # Create output file and wrrite the time=start_time step
        if self.write_output_file:
            self.InitializeOutput()


    def InitializeOutput(self):

        # Carry out this task only in the first rank (for parallel computing)
        data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
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
        data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:

            # Only write if the output is enabled
            if self.write_output_file:

                # Calculate the reaction for the output
                # TODO: this sould be calculated elsewhere or not outputed
                reaction = self.CalculateReaction()

                # Write the output only for the active DOFs
                for index, dof in enumerate(self.available_dofs):
                    if dof in self.active_dofs:
                        with open(self.output_file_name[dof], "a") as results_rigid_body:
                            results_rigid_body.write(str(np.around(self.time, 3)) + " " +
                                                    str(self.x[index,0]) + " " +
                                                    str(self.v[index,0]) + " " +
                                                    str(self.a[index,0]) + " " +
                                                    str(self.x_f[index,0]) + " " +
                                                    str(self.v_f[index,0]) + " " +
                                                    str(self.a_f[index,0]) + " " +
                                                    str(self.x[index,0] - self.x_f[index,0]) + " " +
                                                    str(self.v[index,0] - self.v_f[index,0]) + " " +
                                                    str(self.a[index,0] - self.a_f[index,0]) + " " +
                                                    str(reaction[index]) + "\n")


    def AdvanceInTime(self, current_time):
        # Similar to the Kratos CloneTimeStep function
        # Advances values along the buffer axis (so rolling columns) using numpy's roll
        # Column 0 is the current time step, column 1 the previous one...

        # Variables whith buffer. Column 0 will be overwriten later
        self.x = np.roll(self.x,1,axis=1)
        self.x_f = np.roll(self.x_f,1,axis=1)
        self.v = np.roll(self.v,1,axis=1)
        self.v_f = np.roll(self.v_f,1,axis=1)
        self.a = np.roll(self.a,1,axis=1)
        self.a_f = np.roll(self.a_f,1,axis=1)
        self.total_load = np.roll(self.total_load,1,axis=1)
        self.total_root_point_displ = np.roll(self.total_root_point_displ,1,axis=1)
        self.effective_load = np.roll(self.effective_load,1,axis=1)

        # Variables that need to be reseted. They might not be overwriten later so they
        # need to be zero to avoid values from previous time steps ar not continously used.
        self.prescribed_load = np.zeros(self.system_size)
        self.prescribed_root_point_displ = np.zeros(self.system_size)
        self.external_load = np.zeros(self.system_size)
        self.external_root_point_displ = np.zeros(self.system_size)

        # Update the time of the simulation
        self.time = current_time + self.delta_t
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

        # Creation of the RHS according to the gen-alpha method
        RHS = np.dot(self.M, (self.a1m * self.x[:,1] +
                            self.a2m * self.v[:,1] + self.a3m * self.a[:,1]))
        RHS += np.dot(self.C, (self.a1b * self.x[:,1] +
                            self.a2b * self.v[:,1] + self.a3b * self.a[:,1]))
        RHS += np.dot(self.a1k * self.K, self.x[:,1]) + F

        # Make zero the corrresponding values of RHS so the inactive dofs are not excited
        for index, dof in enumerate(self.available_dofs):
            if dof not in self.active_dofs:
                RHS[index] = 0

        # Solve the solution step and find the new displacements
        self.x[:,0] = np.linalg.solve(self.LHS, RHS)

        # constrained dofs will have the root_point_displacement as a total displacement
        for index, dof in enumerate(self.available_dofs):
            if self.is_constrained[dof]:
                self.x[index,0] = self.total_root_point_displ[index,0]

        # Update velocity and acceleration according to the gen-alpha method
        self.v[:,0] = self.a1v * (self.x[:,0] - self.x[:,1]) + self.a2v * \
                        self.v[:,1] + self.a3v * self.a[:,1]
        self.a[:,0] = self.a1a * (self.x[:,0] - self.x[:,1]) + self.a2a * \
                        self.v[:,1] + self.a3a * self.a[:,1]


    def _UpdateRootPointDisplacement(self, displ):
        # Save root point displacement
        self.x_f[:,0] = displ
        # Update the velocity and acceleration according to the gen-alpha method
        self.v_f[:,0] = self.a1v * (self.x_f[:,0] - self.x_f[:,1]) + self.a2v * \
            self.v_f[:,1] + self.a3v * self.a_f[:,1]
        self.a_f[:,0] = self.a1a * (self.x_f[:,0] - self.x_f[:,1]) + self.a2a * \
            self.v_f[:,1] + self.a3a * self.a_f[:,1]


    def _CalculateEquivalentForceFromRootPointDisplacement(self):
        # Transform the movement of the root point into an equivalent force
        equivalent_force = self.K.dot(self.x_f[:,0]) + self.C.dot(self.v_f[:,0])
        return equivalent_force


    def _CalculateEffectiveLoad(self):

        # Calculate the total load
        self_weight = self.CalculateSelfWeight()
        self.total_load[:,0] = self.external_load + self.prescribed_load + self_weight

        # Calculate the total root point displacement and the equivalent force it generates
        self.total_root_point_displ[:,0] = self.external_root_point_displ + self.prescribed_root_point_displ
        self._UpdateRootPointDisplacement(self.total_root_point_displ[:,0])
        root_point_force = self._CalculateEquivalentForceFromRootPointDisplacement()
        
        # Sum up both loads
        effective_load = self.total_load[:,0] + root_point_force

        return effective_load
    

    def _UpdateVelocityAndAcceleration(self, x):
        # Update the velocity and acceleration according to the gen-alpha method
        self.v[:,0] = self.a1v * (x - self.x[:,1]) + self.a2v * \
                        self.v[:,1] + self.a3v * self.a[:,1]
        self.a[:,0] = self.a1a * (x - self.x[:,1]) + self.a2a * \
                        self.v[:,1] + self.a3a * self.a[:,1]


    def CalculateReaction(self, buffer_idx=0):
        # TODO: Check how does it work with constrained DOFs
        reaction = self.C.dot(self.v[:,buffer_idx] - self.v_f[:,buffer_idx]) \
                + self.K.dot(self.x[:,buffer_idx] - self.x_f[:,buffer_idx])
        return reaction


    def CalculateSelfWeight(self):
        self_weight = self.M.dot(self.modulus_self_weight)
        return self_weight


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

        # Loop through input active DOFs saving the values
        for index, value in enumerate(values):
            # Increase the index so it fits with the angular terms (if necessary)
            if identifier in self.angular_variables:
                index += self.linear_size
            # Save input variables in their corresponding spots
            if self.available_dofs[index] in self.active_dofs:
                if identifier in ["DISPLACEMENT", "ROTATION", "DISPLACEMENTS_ALL"]:
                    self.x[index, buffer_idx] = value
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


    def GetSolutionStepValue(self, identifier, buffer_idx=0):

        # Check that the selected buffer value is OK
        input_check._CheckBufferId(buffer_idx, identifier)

        # Initialize output as a KratosMultiphysics vector
        output = KratosMultiphysics.Vector(self._ExpectedDataSize(identifier))

        # Save all the values from the active DOFs
        for index in range(len(output)):
            out_index = index
            # Increase the index so it fits with the angular terms (if necessary)
            if identifier in self.angular_variables:
                index += self.linear_size
            # Fill the output with its corresponding values
            if identifier in ["DISPLACEMENT", "ROTATION", "DISPLACEMENTS_ALL"]:
                output[out_index] = self.x[index, buffer_idx]
            elif identifier == "VELOCITY":
                output[out_index] = self.v[index, buffer_idx]
            elif identifier == "ACCELERATION":
                output[out_index] = self.a[index, buffer_idx]
            elif identifier in ["REACTION", "REACTION_MOMENT", "REACTION_ALL"]:
                output[out_index] = self.CalculateReaction(buffer_idx=buffer_idx)[index]
            else:
                raise Exception("Identifier is unknown!")
        return output

    
    def _ExpectedDataSize(self, identifier):

        if identifier in self.linear_variables:
            expected_size = self.linear_size
        elif identifier in self.angular_variables:
            expected_size = self.angular_size
        else:
            expected_size = self.system_size

        return expected_size

