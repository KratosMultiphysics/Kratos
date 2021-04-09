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
        self.available_dofs = ['displacement_x', 'displacement_y', 'displacement_z']

        # Check that the activated dofs are not repeated and are among the available ones
        self.active_dofs = self._CheckActiveDofs(parameters)

        # Fill with defaults and check that mandatory fields are given
        self._CheckMandatoryInputParameters(parameters)
        dof_params, sol_params = self._ValidateAndAssignRigidBodySolverDefaults(parameters)
        
        # Safe all the filled data in their respective class variables
        self._InitializeDofsVariables(dof_params)
        self._InitializeSolutionVariables(sol_params)

        # Calculate parameters for the system of equations that performs the time integration.
        # The Left-hand-side (LHS) is always the same. The Right-hand-side (RHS) depends on the forces
        # and previous situation. Here, only the constant part (RHS_matrix) is calculated.
        self.LHS = {}
        self.RHS_matrix = {}
        for dof in self.active_dofs:
            self.LHS[dof] = np.array([[1.0, 0.0, -self.delta_t**2 * self.beta],
                                    [0.0, 1.0, -self.delta_t * self.gamma],
                                    [self.stiffness[dof],
                                    self.damping[dof],
                                    (1-self.alpha_m) * self.mass[dof]]])

            self.RHS_matrix[dof] = np.array([[1.0, self.delta_t, self.delta_t**2 * (0.5 - self.beta)],
                                            [0.0, 1.0, self.delta_t*(1-self.gamma)],
                                            [0.0,
                                            0.0,
                                            -self.alpha_m * self.mass[dof]]])

    def _InitializeDofsVariables(self, dof_params):
        
        self.mass = {}
        self.stiffness = {}
        self.damping = {}
        self.modulus_self_weight = {}

        self.excitation_function_force = {}
        self.excitation_function_root_point_displacement = {}
        self.load_impulse = {}
        self.omega_force = {}
        self.omega_root_point_displacement = {}
        self.amplitude_root_point_displacement = {}
        self.amplitude_force = {}

        self.initial_displacement = {}
        self.initial_velocity = {}
        self.initial_acceleration = {}

        for dof in self.active_dofs:
            self.mass[dof] = dof_params[dof]['system_parameters']['mass'].GetDouble()
            self.stiffness[dof] = dof_params[dof]['system_parameters']['stiffness'].GetDouble()
            self.damping[dof] = dof_params[dof]['system_parameters']['damping'].GetDouble()
            self.modulus_self_weight[dof] = dof_params[dof]['system_parameters']['modulus_self_weight'].GetDouble()

            self.excitation_function_force[dof] = dof_params[dof]["boundary_conditions"]["excitation_function_force"].GetString()
            self.excitation_function_root_point_displacement[dof] = dof_params[dof]["boundary_conditions"]["excitation_function_root_point_displacement"].GetString()
            self.load_impulse[dof] = dof_params[dof]["boundary_conditions"]["load_impulse"].GetDouble()
            self.omega_force[dof] = dof_params[dof]["boundary_conditions"]["omega_force"].GetDouble()
            self.omega_root_point_displacement[dof] = dof_params[dof]["boundary_conditions"]["omega_root_point_displacement"].GetDouble()
            self.amplitude_root_point_displacement[dof] = dof_params[dof]["boundary_conditions"]["amplitude_root_point_displacement"].GetDouble()
            self.amplitude_force[dof] = dof_params[dof]["boundary_conditions"]["amplitude_force"].GetDouble()

            self.initial_displacement[dof] = dof_params[dof]["initial_values"]["displacement"].GetDouble()
            self.initial_velocity[dof] = dof_params[dof]["initial_values"]["velocity"].GetDouble()
            factor = self.load_impulse[dof] - self.stiffness[dof] * self.initial_displacement[dof]
            self.initial_acceleration[dof] = (1/self.mass[dof]) * factor

    def _InitializeSolutionVariables(self, sol_params):
        
        self.alpha_m = sol_params["time_integration_parameters"]["alpha_m"].GetDouble()
        self.beta = 0.25 * (1- self.alpha_m)**2
        self.gamma =  0.50 - self.alpha_m

        self.delta_t = sol_params["time_integration_parameters"]["time_step"].GetDouble()
        self.start_time = sol_params["time_integration_parameters"]["start_time"].GetDouble()
        self.buffer_size = sol_params["solver_parameters"]["buffer_size"].GetInt()
        self.output_file_path = sol_params["output_parameters"]["file_path"].GetString()
        self.write_output_file = sol_params["output_parameters"]["write_output_files"].GetBool()
        
    def Initialize(self):

        self.time = self.start_time

        self.x = {}
        self.x_f = {}
        self.dx = {}
        self.dx_f = {}
        self.root_point_displacement = {}
        self.load_vector = {}

        for dof in self.active_dofs:
            #solution buffer
            self.x[dof] = np.zeros((3, self.buffer_size))
            #values at the root point buffer
            self.x_f[dof] = np.zeros((3, self.buffer_size))

            initial_values = np.array([self.initial_displacement[dof],
                                        self.initial_velocity[dof],
                                        self.initial_acceleration[dof]])
            self.dx[dof] = initial_values
            self.dx_f[dof] = np.zeros(3)

            self.root_point_displacement[dof] = 0.0

            #apply external load as an initial impulse
            self.load_vector[dof] = np.array([0,
                                            0,
                                            self.load_impulse[dof]])

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
                                            "root point displacement" + " " +
                                            "root point velocity" + " " +
                                            "root point acceleration" + " " +
                                            "relative displacement" + " " +
                                            "relative velocity" + " " +
                                            "relative accleration" + " " +
                                            "reaction" + "\n")
            self.OutputSolutionStep()

    def OutputSolutionStep(self):
        data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:
            if self.write_output_file:
                for dof in self.active_dofs:
                    reaction = self.CalculateReaction(dof)
                    with open(self.output_file_name[dof], "a") as results_rigid_body:
                        #outputs results
                        #x and dx contain: [displacement, velocity, acceleration]
                        results_rigid_body.write(str(np.around(self.time, 3)) + " " +
                                                str(self.dx[dof][0]) + " " +
                                                str(self.dx[dof][1]) + " " +
                                                str(self.dx[dof][2]) + " " +
                                                str(self.dx_f[dof][0]) + " " +
                                                str(self.dx_f[dof][1]) + " " +
                                                str(self.dx_f[dof][2]) + " " +
                                                str(self.dx[dof][0] - self.dx_f[dof][0]) + " " +
                                                str(self.dx[dof][1] - self.dx_f[dof][1]) + " " +
                                                str(self.dx[dof][2] - self.dx_f[dof][2]) + " " +
                                                str(reaction) + "\n")


    def AdvanceInTime(self, current_time):
        for dof in self.active_dofs:
            # similar to the Kratos CloneTimeStep function
            # advances values along the buffer axis (so rolling columns) using numpy's roll
            self.x[dof] = np.roll(self.x[dof],1,axis=1)
            self.x_f[dof] = np.roll(self.x_f[dof],1,axis=1)
            # overwriting at the buffer_idx=0 the newest values
            buffer_idx = 0
            self.x[dof][:,buffer_idx] = self.dx[dof]
            self.x_f[dof][:,buffer_idx] = self.dx_f[dof]

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

    def CalculateReaction(self, dof, buffer_idx=0):
        reaction = self.damping[dof] * (self.dx[dof][1] - self.dx_f[dof][1]) \
                + self.stiffness[dof] * (self.dx[dof][0] - self.dx_f[dof][0])
        return reaction

    def CalculateSelfWeight(self, dof):
        self_weight = self.mass[dof] * self.modulus_self_weight[dof]
        return self_weight

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        # Shouldn't it be df?
        output = []
        for dof in self.available_dofs:
            if dof in self.active_dofs:
                if identifier == "DISPLACEMENT":
                    output.append(self.x[dof][:,buffer_idx][0])
                elif identifier == "VELOCITY":
                    output.append(self.x[dof][:,buffer_idx][1])
                elif identifier == "ACCELERATION":
                    output.append(self.x[dof][:,buffer_idx][2])
                elif identifier == "REACTION":
                    output.append(self.CalculateReaction(dof))
                elif identifier == "VOLUME_ACCELERATION":
                    output.append(self.CalculateSelfWeight(dof))
                else:
                    raise Exception("Identifier is unknown!")
            else:
                output.append(0.0)
        return output

    def SetSolutionStepValue(self, identifier, values, buffer_idx=0):
        for val_id, value in enumerate(values):
            if self.available_dofs[val_id] in self.active_dofs:
                if identifier == "DISPLACEMENT":
                    self.x[self.available_dofs[val_id]][:,buffer_idx][0] = value
                elif identifier == "VELOCITY":
                    self.x[self.available_dofs[val_id]][:,buffer_idx][1] = value
                elif identifier == "ACCELERATION":
                    self.x[self.available_dofs[val_id]][:,buffer_idx][2] = value
                elif identifier == "LOAD":
                    self.load_vector[self.available_dofs[val_id]][-1] = 0.0
                    self.load_vector[self.available_dofs[val_id]][-1] = value
                elif identifier == "ROOT_POINT_DISPLACEMENT":
                    self.root_point_displacement[self.available_dofs[val_id]] = 0.0
                    self.root_point_displacement[self.available_dofs[val_id]] = value
                else:
                    raise Exception("Identifier is unknown!")

    def _CheckActiveDofs(self, parameters):

        active_dofs = []
        for single_dof_parameters in parameters["active_dof_list"]:

            if "dof" not in single_dof_parameters:
                msg = '"dof" must be specified in all the items in the list "active_dof_list". '
                msg += 'It represents the name of the degree of freedom that is activated.'
                raise Exception(msg)
            dof = single_dof_parameters["dof"]

            if dof not in self.available_dofs:
                msg = 'The degree of freedom "' + dof + '" is not among the available ones. '
                msg += 'Chose one of the following: ' + str(self.available_dofs)[1:-1] + '.'
                raise Exception(msg)

            if dof in active_dofs:
                msg = 'The degree of freedom "' + dof + '" is repeated in the project parameters. '
                raise Exception(msg)
            active_dofs.append(dof)
        
        return active_dofs

    def _ValidateAndAssignRigidBodySolverDefaults(self, parameters):

        default_single_dof_parameters = KratosMultiphysics.Parameters('''{
            "dof": "displacement_x",
            "system_parameters":{
                "mass"      : 100.0,
                "stiffness" : 4000.0,
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
                "alpha_m"   : -0.3,
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
        for single_dof_parameters in parameters["active_dof_list"]:

            dof = single_dof_parameters["dof"]

            single_dof_parameters = KratosMultiphysics.Parameters(json.dumps(single_dof_parameters))
            single_dof_parameters.RecursivelyValidateAndAssignDefaults(default_single_dof_parameters)
            
            dof_parameters[dof] = single_dof_parameters
        
        solution_parameters = KratosMultiphysics.Parameters(json.dumps(parameters["solution_parameters"]))
        solution_parameters.RecursivelyValidateAndAssignDefaults(default_solution_parameters)

        return dof_parameters, solution_parameters

    def _CheckMandatoryInputParameters(self, parameters):

        for key in ["active_dof_list", "solution_parameters"]:
            if key not in parameters:
                msg = 'The key "' + key + '" was not found in the project parameters '
                msg += 'and it is necessary to run configure the RigidBodySolver.'
                raise Exception(msg)
        
        if len(parameters["active_dof_list"]) == 0:
            msg = 'At least an active degree of freedom is needed to use the solver '
            msg += ' and none where provided in "active_dof_list".'
            raise Exception(msg)
        
        msg = '"time_step" should be given as par of "time_integration_parameters" '
        msg += 'in "solution_parameters".'
        if "time_integration_parameters" not in parameters["solution_parameters"]:
            raise Exception(msg)
        else:
            if "time_step" not in parameters["solution_parameters"]["time_integration_parameters"]:
                raise Exception(msg)