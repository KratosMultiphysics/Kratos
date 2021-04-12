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

        # Check that the activated dofs are not repeated and are among the available ones
        self.active_dofs = self._CheckActiveDofs(parameters)
        self.system_size = len(self.active_dofs)

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
        self.excitation_function_root_point_displacement = {}
        self.load_impulse = {}
        self.omega_force = {}
        self.omega_root_point_displacement = {}
        self.amplitude_root_point_displacement = {}
        self.amplitude_force = {}

        for dof in self.active_dofs:

            self.excitation_function_force[dof] = dof_params[dof]["boundary_conditions"]["excitation_function_force"].GetString()
            self.excitation_function_root_point_displacement[dof] = dof_params[dof]["boundary_conditions"]["excitation_function_root_point_displacement"].GetString()
            self.load_impulse[dof] = dof_params[dof]["boundary_conditions"]["load_impulse"].GetDouble()
            self.omega_force[dof] = dof_params[dof]["boundary_conditions"]["omega_force"].GetDouble()
            self.omega_root_point_displacement[dof] = dof_params[dof]["boundary_conditions"]["omega_root_point_displacement"].GetDouble()
            self.amplitude_root_point_displacement[dof] = dof_params[dof]["boundary_conditions"]["amplitude_root_point_displacement"].GetDouble()
            self.amplitude_force[dof] = dof_params[dof]["boundary_conditions"]["amplitude_force"].GetDouble()

    def _InitializeSolutionVariables(self, sol_params):
        
        rho_inf = sol_params["time_integration_parameters"]["rho_inf"].GetDouble()
        self.alpha_f = rho_inf / (rho_inf + 1)
        self.alpha_m = (2*rho_inf - 1) / (rho_inf + 1)
        self.beta = 0.25 * (1- self.alpha_m + self.alpha_f)**2
        self.gamma =  0.50 - self.alpha_m + self.alpha_f

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

        for index, dof in enumerate(self.active_dofs):

            self.M[index][index] = dof_params[dof]['system_parameters']['mass'].GetDouble()
            self.C[index][index] = dof_params[dof]['system_parameters']['damping'].GetDouble()
            self.K[index][index] = dof_params[dof]['system_parameters']['stiffness'].GetDouble()
            self.modulus_self_weight[index] = dof_params[dof]['system_parameters']['modulus_self_weight'].GetDouble()

            self.initial_displacement[index] = dof_params[dof]["initial_values"]["displacement"].GetDouble()
            self.initial_velocity[index] = dof_params[dof]["initial_values"]["velocity"].GetDouble()
            factor = self.load_impulse[dof] - self.K[index][index] * self.initial_displacement[index]
            self.initial_acceleration[index] = (1/self.M[index][index]) * factor

    def _InitializeGeneralizedAlphaParameters(self):
        
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

        # TODO: Check notation of the following variables to make it intuitive
        # Others
        self.root_point_displacement = np.zeros(self.system_size)
        self.load_vector = np.zeros(self.system_size)

        # Apply initial conditions
        self.x[:,0] = self.initial_displacement
        self.v[:,0] = self.initial_velocity
        self.a[:,0] = self.initial_acceleration

        #Apply external load as an initial impulse
        for index, dof in enumerate(self.active_dofs):
            self.load_vector[index] = self.load_impulse[dof]

        # Create output file and wrrite the time=start_time step
        if self.write_output_file:
            self.InitializeOutput()

    def InitializeOutput(self):
        
        data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:
            self.output_file_name = {}
            for dof in self.active_dofs:
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
                for index, dof in enumerate(self.active_dofs):
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
        for dof in self.active_dofs:
            # similar to the Kratos CloneTimeStep function
            # advances values along the buffer axis (so rolling columns) using numpy's roll
            self.x = np.roll(self.x,1,axis=1)
            self.x_f = np.roll(self.x_f,1,axis=1)
            self.v = np.roll(self.v,1,axis=1)
            self.v_f = np.roll(self.v_f,1,axis=1)
            self.a = np.roll(self.a,1,axis=1)
            self.a_f = np.roll(self.a_f,1,axis=1)

        self.time = current_time + self.delta_t
        return self.time

    def CalculateEquivalentForceFromRootPointExcitation(self, d_f_excitation):
        b_f = {}
        for dof in self.active_dofs:
            d_f = d_f_excitation[dof] + self.root_point_displacement[dof]
            v_f = self.x_f[dof][1,0] + self.delta_t * (self.gamma * d_f + (1-self.gamma) * self.x_f[dof][2,0])
            a_f = 1/(self.delta_t**2 * self.beta) * (d_f - self.x_f[dof][0,1])\
                - 1/(self.delta_t * self.beta) * self.x_f[dof][1,0]\
                + (1-1/(2*self.beta)) * self.x_f[dof][2,0]
            self.dx_f[dof] = np.array([d_f, v_f, a_f])
            b_f[dof] = np.array([0.0, 0.0, d_f * self.stiffness[dof] + v_f * self.damping[dof]])
        return b_f

    def ApplyRootPointExcitation(self):
        excitation = {}
        for dof in self.active_dofs:
            scope_vars = {'t' : self.time, 'omega': self.omega_root_point_displacement[dof], 'A': self.amplitude_root_point_displacement[dof]}
            excitation[dof] = GenericCallFunction(self.excitation_function_root_point_displacement[dof], scope_vars, check=False)
        return excitation

    def ApplyForceExcitation(self):
        excitation = {}
        for dof in self.active_dofs:
            scope_vars = {'t' : self.time, 'omega': self.omega_force[dof], 'A': self.amplitude_force[dof]}
            excitation[dof] = GenericCallFunction(self.excitation_function_force[dof], scope_vars, check=False)
        return excitation

    def SolveSolutionStep(self):
        #external load
        excitation_load = self.ApplyForceExcitation()
        #root point displacement
        d_f_excitation = self.ApplyRootPointExcitation()
        #root point force
        b_f = self.CalculateEquivalentForceFromRootPointExcitation(d_f_excitation)
        for dof in self.active_dofs:
            RHS = self.RHS_matrix[dof] @ self.x[dof][:,0]
            self.load_vector[dof][-1] += excitation_load[dof] # Why? Shouldn't it be just a "="?
            # print("external load= ", self.load_vector[-1])
            RHS += self.load_vector[dof]
            RHS += b_f[dof]
            self.dx[dof] = np.linalg.solve(self.LHS[dof], RHS)

    def CalculateReaction(self, buffer_idx=0):
        reaction = self.C.dot(self.v[:,buffer_idx][1] - self.v_f[:,buffer_idx][1]) \
                + self.K.dot(self.x[:,buffer_idx][0] - self.x_f[:,buffer_idx][0])
        return reaction

    def CalculateSelfWeight(self, dof):
        self_weight = self.M.dot(self.modulus_self_weight)
        return self_weight

    def SetSolutionStepValue(self, identifier, values, buffer_idx=0):
        self._CheckBufferId(buffer_idx,identifier)
        for val_index, value in enumerate(values):
            if identifier == "MOMENT":
                val_index += 3
            if self.available_dofs[val_index] in self.active_dofs:
                reduced_index = self._ReduceIndex(val_index)
                if identifier == "DISPLACEMENT":
                    self.x[reduced_index, buffer_idx] = value
                elif identifier == "VELOCITY":
                    self.v[reduced_index, buffer_idx] = value
                elif identifier == "ACCELERATION":
                    self.a[reduced_index, buffer_idx] = value
                elif identifier == "FORCE" or identifier == "MOMENT":
                    self.load_vector[reduced_index] = value
                elif identifier == "ROOT_POINT_DISPLACEMENT":
                    self.root_point_displacement[reduced_index] = value
                else:
                    raise Exception("Identifier is unknown!")

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        self._CheckBufferId(buffer_idx, identifier)
        output = np.zeros(len(self.available_dofs))
        for index, dof in enumerate(self.available_dofs):
            if dof in self.active_dofs:
                reduced_index = self._ReduceIndex(index)
                if identifier == "DISPLACEMENT":
                    output[index] = self.x[reduced_index, buffer_idx]
                elif identifier == "VELOCITY":
                    output[index] = self.v[reduced_index, buffer_idx]
                elif identifier == "ACCELERATION":
                    output[index] = self.a[reduced_index, buffer_idx]
                elif identifier == "REACTION":
                    output[index] = self.CalculateReaction(buffer_idx=buffer_idx)[reduced_index]
                elif identifier == "VOLUME_ACCELERATION":
                    output[index] = self.CalculateSelfWeight()[reduced_index]
                else:
                    raise Exception("Identifier is unknown!")
        return output

    def _CheckBufferId(self, buffer_idx, identifier):
        if identifier not in ["DISPLACEMENT","VELOCITY","ACCELERATION"] and buffer_idx != 0:
            msg = 'The buffer_idx can only be 0 for the variable "' + identifier + '".'
            raise Exception(msg)

    def _ReduceIndex(index):
        dof = self.available_dofs[index]
        if dof in self.active_dofs:
            return self.active_dofs.index(dof)
    '''
    def _ExpandIndex(index):
        dof = self.active_dofs[index]
        return self.available_dofs.index(dof)
    '''
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
        
        return active_dofs

    def _CheckMandatoryInputParameters(self, parameters):

        for key in ["active_dofs", "solution_parameters"]:
            if key not in parameters:
                msg = 'The key "' + key + '" was not found in the project parameters '
                msg += 'and it is necessary to configure the RigidBodySolver.'
                raise Exception(msg)
        
        if len(parameters["active_dofs"]) == 0:
            msg = 'At least an active degree of freedom is needed to use the solver '
            msg += ' and none where provided in "active_dofs".'
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
                "velocity"      : 0.0,
                "acceleration"  : 0.0
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
        for dof in parameters["active_dofs"]:

            single_dof_parameters = parameters["active_dofs"][dof]

            single_dof_parameters = KratosMultiphysics.Parameters(json.dumps(single_dof_parameters))
            single_dof_parameters.RecursivelyValidateAndAssignDefaults(default_single_dof_parameters)
            
            dof_parameters[dof] = single_dof_parameters
        
        solution_parameters = KratosMultiphysics.Parameters(json.dumps(parameters["solution_parameters"]))
        solution_parameters.RecursivelyValidateAndAssignDefaults(default_solution_parameters)

        return dof_parameters, solution_parameters