# CoSimulation imports
from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction
import KratosMultiphysics

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
            with open(input_name,'r') as ProjectParameters:
                parameters = json.load(ProjectParameters)
        else:
            raise Exception("The input has to be provided as a dict or a string")

        # Degrees of freedom that can be activated from the parameters
        self.available_dofs = [
            'displacement_x',
            'displacement_y',
            'displacement_z',
            'rotation_x',
            'rotation_y',
            'rotation_z'
        ]
        self.system_size = len(self.available_dofs)

        # Check that the activated dofs are not repeated and are among the available ones
        self.active_dofs = self._CheckActiveDofs(parameters)

        # Fill with defaults and check that mandatory fields are given
        self._CheckMandatoryInputParameters(parameters)
        dof_params, sol_params = self._ValidateAndAssignRigidBodySolverDefaults(parameters)
        
        # Safe all the filled data in their respective class variables
        self._InitializeDofsVariables(dof_params)
        self._InitializeSolutionVariables(sol_params)

        # Create the stiffness, damping and mass matrix
        self._InitializeStructuralParameters(dof_params)

        # Prepare the parameters for the generalized-alpha method
        self._InitializeGeneralizedAlphaParameters()

    def _InitializeDofsVariables(self, dof_params):

        self.excitation_function_force = {}
        self.excitation_function_root_point_displ = {}
        self.load_impulse = {}
        self.omega_force = {}
        self.omega_root_point_displ = {}
        self.amplitude_root_point_displ = {}
        self.amplitude_force = {}

        for dof in self.available_dofs:

            self.excitation_function_force[dof] = dof_params[dof]["boundary_conditions"]["excitation_function_force"].GetString()
            self.excitation_function_root_point_displ[dof] = dof_params[dof]["boundary_conditions"]["excitation_function_root_point_displacement"].GetString()
            self.load_impulse[dof] = dof_params[dof]["boundary_conditions"]["load_impulse"].GetDouble()
            self.omega_force[dof] = dof_params[dof]["boundary_conditions"]["omega_force"].GetDouble()
            self.omega_root_point_displ[dof] = dof_params[dof]["boundary_conditions"]["omega_root_point_displacement"].GetDouble()
            self.amplitude_root_point_displ[dof] = dof_params[dof]["boundary_conditions"]["amplitude_root_point_displacement"].GetDouble()
            self.amplitude_force[dof] = dof_params[dof]["boundary_conditions"]["amplitude_force"].GetDouble()

    def _InitializeSolutionVariables(self, sol_params):
        
        self.rho_inf = sol_params["time_integration_parameters"]["rho_inf"].GetDouble()

        self.delta_t = sol_params["time_integration_parameters"]["time_step"].GetDouble()
        self.start_time = sol_params["time_integration_parameters"]["start_time"].GetDouble()
        self.buffer_size = sol_params["solver_parameters"]["buffer_size"].GetInt()
        self.output_file_path = sol_params["output_parameters"]["file_path"].GetString()
        self.write_output_file = sol_params["output_parameters"]["write_output_files"].GetBool()

    def _InitializeStructuralParameters(self, dof_params):

        self.M = np.zeros((self.system_size,self.system_size)) # Mass matrix
        self.C = np.zeros((self.system_size,self.system_size)) # Damping matrix
        self.K = np.zeros((self.system_size,self.system_size)) # Stiffness matrix
        self.modulus_self_weight = np.zeros(self.system_size)

        self.initial_displacement = np.zeros(self.system_size)
        self.initial_velocity = np.zeros(self.system_size)
        self.initial_acceleration = np.zeros(self.system_size)

        for index, dof in enumerate(self.available_dofs):

            self.M[index][index] = dof_params[dof]['system_parameters']['mass'].GetDouble()
            self.C[index][index] = dof_params[dof]['system_parameters']['damping'].GetDouble()
            self.K[index][index] = dof_params[dof]['system_parameters']['stiffness'].GetDouble()
            self.modulus_self_weight[index] = dof_params[dof]['system_parameters']['modulus_self_weight'].GetDouble()

            self.initial_displacement[index] = dof_params[dof]["initial_values"]["displacement"].GetDouble()
            self.initial_velocity[index] = dof_params[dof]["initial_values"]["velocity"].GetDouble()
            factor = self.load_impulse[dof] - self.K[index][index] * self.initial_displacement[index]
            self.initial_acceleration[index] = (1/self.M[index][index]) * factor

    def _InitializeGeneralizedAlphaParameters(self):

        self.alpha_f = self.rho_inf / (self.rho_inf + 1)
        self.alpha_m = (2*self.rho_inf - 1) / (self.rho_inf + 1)
        self.beta = 0.25 * (1- self.alpha_m + self.alpha_f)**2
        self.gamma =  0.50 - self.alpha_m + self.alpha_f
        
        # coefficients for LHS
        self.a1h = (1.0 - self.alpha_m) / (self.beta * self.delta_t**2)
        self.a2h = (1.0 - self.alpha_f) * self.gamma / (self.beta * self.delta_t)
        self.a3h = 1.0 - self.alpha_f

        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.delta_t
        self.a3m = (1.0 - self.alpha_m - 2.0 * self.beta) / (2.0 * self.beta)

        # coefficients for damping
        self.a1b = (1.0 - self.alpha_f) * self.gamma / (self.beta * self.delta_t)
        self.a2b = (1.0 - self.alpha_f) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alpha_f) * (0.5 * self.gamma / self.beta - 1.0) * self.delta_t

        # coefficient for stiffness
        self.a1k = -1.0 * self.alpha_f

        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.delta_t)
        self.a2v = 1.0 - self.gamma / self.beta
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.delta_t

        # coefficients for acceleration update
        self.a1a = self.a1v / (self.delta_t * self.gamma)
        self.a2a = -1.0 / (self.beta * self.delta_t)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)

        # In a linear case, the Left-Hand Side can be calculated now
        self.LHS = self.a1h * self.M + self.a2h * self.C + self.a3h * self.K
        
    def Initialize(self):

        self.time = self.start_time

        # Initialize total displacement, velocity and acceleration
        self.x = np.zeros((self.system_size, self.buffer_size))
        self.v = np.zeros((self.system_size, self.buffer_size))
        self.a = np.zeros((self.system_size, self.buffer_size))

        # Initialize the displacement, velocity and accelerration of the root point
        self.x_f = np.zeros((self.system_size, self.buffer_size))
        self.v_f = np.zeros((self.system_size, self.buffer_size))
        self.a_f = np.zeros((self.system_size, self.buffer_size))

        # Others
        self.total_root_point_displ = np.zeros((self.system_size, self.buffer_size))
        self.total_load = np.zeros((self.system_size, self.buffer_size))
        self.external_root_point_displ = np.zeros(self.system_size)
        self.external_load = np.zeros(self.system_size)
        self.effective_load = np.zeros((self.system_size, self.buffer_size))

        # Apply initial conditions
        self.x[:,0] = self.initial_displacement
        self.v[:,0] = self.initial_velocity
        self.a[:,0] = self.initial_acceleration

        #Apply external load as an initial impulse
        for index, dof in enumerate(self.available_dofs):
            self.total_load[index,0] = self.load_impulse[dof]
        self.effective_load[:,0] = self.total_load[index,0]

        # Create output file and wrrite the time=start_time step
        if self.write_output_file:
            self.InitializeOutput()

    def InitializeOutput(self):

        data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:
            self.output_file_name = {}
            for dof in self.available_dofs: #TODO: create only files for active dofs
                self.output_file_name[dof] = os.path.join(self.output_file_path, dof + '.dat')
                if os.path.isfile(self.output_file_name[dof]):
                    os.remove(self.output_file_name[dof])
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
            self.OutputSolutionStep()

    def OutputSolutionStep(self):
        data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:
            if self.write_output_file:
                reaction = self.CalculateReaction()
                for index, dof in enumerate(self.available_dofs):
                    # TODO: Ad a filter here so only active dofs are written
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
        # similar to the Kratos CloneTimeStep function
        # advances values along the buffer axis (so rolling columns) using numpy's roll
        self.x = np.roll(self.x,1,axis=1)
        self.x_f = np.roll(self.x_f,1,axis=1)
        self.v = np.roll(self.v,1,axis=1)
        self.v_f = np.roll(self.v_f,1,axis=1)
        self.a = np.roll(self.a,1,axis=1)
        self.a_f = np.roll(self.a_f,1,axis=1)
        self.total_load = np.roll(self.total_load,1,axis=1)
        self.total_root_point_displ = np.roll(self.total_root_point_displ,1,axis=1)
        self.effective_load = np.roll(self.effective_load,1,axis=1)

        self.time = current_time + self.delta_t
        return self.time

    def CalculateEquivalentForceFromRootPointExcitation(self, displ):
        self.x_f[:,0] = displ
        self.v_f[:,0] = self.v_f[:,1] + self.delta_t * (self.gamma * displ + (1-self.gamma) * self.a_f[:,1])
        self.a_f[:,0] = 1/(self.delta_t**2 * self.beta) * (displ - self.x_f[:,1])\
            - 1/(self.delta_t * self.beta) * self.v_f[:,1]\
            + (1-1/(2*self.beta)) * self.a_f[:,1]
        equivalent_force = self.K.dot(self.x_f[:,0]) + self.C.dot(self.v_f[:,0])
        return equivalent_force

    def ApplyRootPointExcitation(self):
        root_point_excitation = np.zeros(self.system_size)
        for index, dof in enumerate(self.available_dofs):
            scope_vars = {'t' : self.time, 'omega': self.omega_root_point_displ[dof], 'A': self.amplitude_root_point_displ[dof]}
            root_point_excitation[index] = GenericCallFunction(self.excitation_function_root_point_displ[dof], scope_vars, check=False)
        return root_point_excitation

    def ApplyForceExcitation(self):
        force_excitation = np.zeros(self.system_size)
        for index, dof in enumerate(self.available_dofs):
            scope_vars = {'t' : self.time, 'omega': self.omega_force[dof], 'A': self.amplitude_force[dof]}
            force_excitation[index] = GenericCallFunction(self.excitation_function_force[dof], scope_vars, check=False)
        return force_excitation

    def SolveSolutionStep(self):
        #external load
        prescribed_load = self.ApplyForceExcitation()
        self.total_load[:,0] = self.external_load + prescribed_load
        #root point displacement
        prescribed_root_point_displ = self.ApplyRootPointExcitation()
        self.total_root_point_displ[:,0] = self.external_root_point_displ + prescribed_root_point_displ
        #root point force
        root_point_force = self.CalculateEquivalentForceFromRootPointExcitation(self.total_root_point_displ[:,0])
        #equivalent force
        self.effective_load[:,0] = self.total_load[:,0] + root_point_force

        F = (1.0 - self.alpha_f) * self.effective_load[:,0] + self.alpha_f * self.effective_load[:,1]

        # system: in matrix form
        RHS = np.dot(self.M, (self.a1m * self.x[:,1] +
                            self.a2m * self.v[:,1] + self.a3m * self.a[:,1]))
        RHS += np.dot(self.C, (self.a1b * self.x[:,1] +
                            self.a2b * self.v[:,1] + self.a3b * self.a[:,1]))
        RHS += np.dot(self.a1k * self.K, self.x[:,1]) + F

        # Make zero the corrresponding values of RHS so the not active dofs are not excited
        for index, dof in enumerate(self.available_dofs):
            if dof not in self.active_dofs:
                RHS[index] = 0

        self.x[:,0] = np.linalg.solve(self.LHS, RHS)
        self.v[:,0] = self.UpdateVelocity(self.x[:,0])
        self.a[:,0] = self.UpdateAcceleration(self.x[:,0])
    
    def UpdateVelocity(self, x):
        v = self.a1v * (x - self.x[:,1]) + self.a2v * \
            self.v[:,1] + self.a3v * self.a[:,1]
        return v

    def UpdateAcceleration(self, x):
        a = self.a1a * (x - self.x[:,1]) + self.a2a * \
            self.v[:,1] + self.a3a * self.a[:,1]
        return a

    def CalculateReaction(self, buffer_idx=0):
        reaction = self.C.dot(self.v[:,buffer_idx] - self.v_f[:,buffer_idx]) \
                + self.K.dot(self.x[:,buffer_idx] - self.x_f[:,buffer_idx])
        return reaction

    def CalculateSelfWeight(self):
        self_weight = self.M.dot(self.modulus_self_weight)
        return self_weight

    def SetSolutionStepValue(self, identifier, values, buffer_idx=0):
        self._CheckBufferId(buffer_idx,identifier)
        for index, value in enumerate(values):
            if identifier == "MOMENT":
                index += 3
            # Maybe the following "if" is not necessary anymore
            if self.available_dofs[index] in self.active_dofs:
                if identifier == "DISPLACEMENT":
                    self.x[index, buffer_idx] = value
                elif identifier == "VELOCITY":
                    self.v[index, buffer_idx] = value
                elif identifier == "ACCELERATION":
                    self.a[index, buffer_idx] = value
                elif identifier in ["FORCE","MOMENT","RESULTANT"]:
                    self.external_load[index] = value
                elif identifier == "ROOT_POINT_DISPLACEMENT":
                    self.external_root_point_displ[index] = value
                else:
                    raise Exception("Identifier is unknown!")

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        self._CheckBufferId(buffer_idx, identifier)
        output = KratosMultiphysics.Vector(self.system_size)
        for index, dof in enumerate(self.available_dofs):
            if dof in self.active_dofs:
                if identifier == "DISPLACEMENT":
                    output[index] = self.x[index, buffer_idx]
                elif identifier == "VELOCITY":
                    output[index] = self.v[index, buffer_idx]
                elif identifier == "ACCELERATION":
                    output[index] = self.a[index, buffer_idx]
                elif identifier == "REACTION":
                    output[index] = self.CalculateReaction(buffer_idx=buffer_idx)[index]
                elif identifier == "VOLUME_WEIGHT":
                    output[index] = self.CalculateSelfWeight()[index]
                else:
                    raise Exception("Identifier is unknown!")
        return output

    def _CheckBufferId(self, buffer_idx, identifier):
        if identifier in ["VOLUME_WEIGHT", "ROOT_POINT_DISPLACEMENT","FORCE","MOMENT","RESULTANT"] and buffer_idx != 0:
            msg = 'The buffer_idx can only be 0 for the variable "' + identifier + '".'
            raise Exception(msg)
    
    def _CheckActiveDofs(self, parameters):

        active_dofs = []
        for dof in parameters["active_dofs"]:

            if dof not in self.available_dofs:
                msg = 'The degree of freedom "' + dof + '" is not among the available ones. '
                msg += 'Chose one of the following: ' + str(self.available_dofs)[1:-1] + '.'
                raise Exception(msg)

            if dof in active_dofs:
                msg = 'The active degree of freedom "' + dof + '" is repeated in the project parameters. '
                raise Exception(msg)
            active_dofs.append(dof)
        
        if len(active_dofs) == 0:
            msg = 'At least an active degree of freedom is needed to use the solver '
            msg += ' and none where provided in "active_dofs".'
            raise Exception(msg)
        
        return active_dofs

    def _CheckMandatoryInputParameters(self, parameters):

        for key in ["active_dofs", "solution_parameters"]:
            if key not in parameters:
                msg = 'The key "' + key + '" was not found in the project parameters '
                msg += 'and it is necessary to configure the RigidBodySolver.'
                raise Exception(msg)
        
        msg = '"time_step" should be given as par of "time_integration_parameters" '
        msg += 'in "solution_parameters".'
        if "time_integration_parameters" not in parameters["solution_parameters"]:
            raise Exception(msg)
        else:
            if "time_step" not in parameters["solution_parameters"]["time_integration_parameters"]:
                raise Exception(msg)
    
    def _ValidateAndAssignRigidBodySolverDefaults(self, parameters):

        default_single_dof_parameters = KratosMultiphysics.Parameters('''{
            "blocked": false,
            "system_parameters":{
                "mass"      : 1.0,
                "stiffness" : 1.0,
                "damping"   : 0.0,
                "modulus_self_weight": 0.0
            },
            "initial_values":{
                "displacement"  : 0.0,
                "velocity"      : 0.0
            },
            "boundary_conditions":{
                "load_impulse": 0.0,
                "omega_force": 0.0,
                "omega_root_point_displacement": 0.0,
                "excitation_function_force": "A * sin(omega * t)",
                "excitation_function_root_point_displacement": "A * sin(omega * t)",
                "amplitude_root_point_displacement": 0.0,
                "amplitude_force": 0.0
            }
        }''')
        default_solution_parameters = KratosMultiphysics.Parameters('''{
            "time_integration_parameters":{
                "rho_inf"   : 0.16,
                "start_time": 0.0,
                "time_step" : 0.05
            },
            "solver_parameters":{
                "buffer_size"   : 2
            },
            "output_parameters":{
                "write_output_files": true,
                "file_path" : "results/rigid_body"
            }
        }''')

        dof_parameters = {}
        for dof in self.available_dofs:

            if dof not in parameters["active_dofs"]:
                parameters["active_dofs"][dof] = {}

            single_dof_parameters = parameters["active_dofs"][dof]

            single_dof_parameters = KratosMultiphysics.Parameters(json.dumps(single_dof_parameters))
            single_dof_parameters.RecursivelyValidateAndAssignDefaults(default_single_dof_parameters)
            
            dof_parameters[dof] = single_dof_parameters
        
        solution_parameters = KratosMultiphysics.Parameters(json.dumps(parameters["solution_parameters"]))
        solution_parameters.RecursivelyValidateAndAssignDefaults(default_solution_parameters)

        return dof_parameters, solution_parameters