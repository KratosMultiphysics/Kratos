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

        self.available_dofs = ['displacement_x', 'displacement_y', 'displacement_z']

        self.dof_params, self.sol_params = self._ValidateAndAssignRigidBodySolverDefaults(parameters)
        
        self.mass = parameters["system_parameters"]["mass"]
        self.stiffness = parameters["system_parameters"]["stiffness"]
        self.damping = parameters["system_parameters"]["damping"]
        self.modulus_self_weight = parameters["system_parameters"]["modulus_self_weight"]

        self.alpha_m = parameters["time_integration_parameters"]["alpha_m"]
        self.delta_t = parameters["time_integration_parameters"]["time_step"]
        self.start_time = parameters["time_integration_parameters"]["start_time"]

        self.initial_displacement = parameters["initial_values"]["displacement"]
        self.initial_velocity = parameters["initial_values"]["velocity"]

        self.excitation_function_force = parameters["boundary_conditions"]["excitation_function_force"]
        self.excitation_function_root_point_displacement = parameters["boundary_conditions"]["excitation_function_root_point_displacement"]
        self.load_impulse = parameters["boundary_conditions"]["load_impulse"]
        self.omega_force = parameters["boundary_conditions"]["omega_force"]
        self.omega_root_point_displacement = parameters["boundary_conditions"]["omega_root_point_displacement"]
        self.amplitude_root_point_displacement = parameters["boundary_conditions"]["amplitude_root_point_displacement"]
        self.amplitude_force = parameters["boundary_conditions"]["amplitude_force"]

        #calculate initial acceleration
        factor = self.load_impulse - self.stiffness * self.initial_displacement
        self.initial_acceleration = (1/self.mass) * factor

        self.beta = 0.25 * (1- self.alpha_m)**2
        self.gamma =  0.50 - self.alpha_m

        self.LHS = np.array([[1.0, 0.0, -self.delta_t**2 * self.beta],
                             [0.0, 1.0, -self.delta_t * self.gamma],
                             [self.stiffness,
                              self.damping,
                               (1-self.alpha_m) * self.mass]])

        self.RHS_matrix = np.array([[1.0, self.delta_t, self.delta_t**2 * (0.5 - self.beta)],
                                    [0.0, 1.0, self.delta_t*(1-self.gamma)],
                                    [0.0,
                                     0.0,
                                     -self.alpha_m * self.mass]])

        self.buffer_size = parameters["solver_parameters"]["buffer_size"]
        self.output_file_name = parameters["output_parameters"]["file_name"]
        self.write_output_file = parameters["output_parameters"]["write_output_file"]
        
    def Initialize(self):
        #solution buffer
        self.x = np.zeros((3, self.buffer_size))
        #values at the root point buffer
        self.x_f = np.zeros((3, self.buffer_size))

        initial_values = np.array([self.initial_displacement,
                                   self.initial_velocity,
                                   self.initial_acceleration])
        self.dx = initial_values
        self.dx_f = np.zeros(3)
        self.time = self.start_time

        self.root_point_displacement = 0.0

        #x and dx contain: [displacement, velocity, acceleration]
        if self.write_output_file:
            if os.path.isfile(self.output_file_name):
                os.remove(self.output_file_name)
            self.InitializeOutput()

        #apply external load as an initial impulse
        self.load_vector = np.array([0,
                                     0,
                                     self.load_impulse])

    def InitializeOutput(self):
        data_comm = KratosMultiphysics.DataCommunicator.GetDefault()
        if data_comm.Rank()==0:
            with open(self.output_file_name, "w") as results_rigid_body:
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
            reaction = self.CalculateReaction()
            if self.write_output_file:
                with open(self.output_file_name, "a") as results_rigid_body:
                    #outputs results
                    results_rigid_body.write(str(np.around(self.time, 3)) + " " +
                                    str(self.dx[0]) + " " +
                                    str(self.dx[1]) + " " +
                                    str(self.dx[2]) + " " +
                                    str(self.dx_f[0]) + " " +
                                    str(self.dx_f[1]) + " " +
                                    str(self.dx_f[2]) + " " +
                                    str(self.dx[0] - self.dx_f[0]) + " " +
                                    str(self.dx[1] - self.dx_f[1]) + " " +
                                    str(self.dx[2] - self.dx_f[2]) + " " +
                                    str(reaction) + "\n")


    def AdvanceInTime(self, current_time):
        # similar to the Kratos CloneTimeStep function
        # advances values along the buffer axis (so rolling columns) using numpy's roll
        self.x = np.roll(self.x,1,axis=1)
        self.x_f = np.roll(self.x_f,1,axis=1)
        # overwriting at the buffer_idx=0 the newest values
        buffer_idx = 0
        self.x[:,buffer_idx] = self.dx
        self.x_f[:,buffer_idx] = self.dx_f

        self.time = current_time + self.delta_t
        return self.time

    def CalculateEquivalentForceFromRootPointExcitation(self, d_f):
        #d_f = self.root_point_displacement
        v_f = self.x_f[1,0] + self.delta_t * (self.gamma * d_f + (1-self.gamma) * self.x_f[2,0])
        a_f = 1/(self.delta_t**2 * self.beta) * (d_f - self.x_f[0,1])\
            - 1/(self.delta_t * self.beta) * self.x_f[1,0]\
            + (1-1/(2*self.beta)) * self.x_f[2,0]
        self.dx_f = np.array([d_f, v_f, a_f])
        b_f = np.array([0.0, 0.0, d_f * self.stiffness + v_f * self.damping])
        return b_f

    def ApplyRootPointExcitation(self):
        scope_vars = {'t' : self.time, 'omega': self.omega_root_point_displacement, 'A': self.amplitude_root_point_displacement}
        return GenericCallFunction(self.excitation_function_root_point_displacement, scope_vars, check=False)

    def ApplyForceExcitation(self):
        scope_vars = {'t' : self.time, 'omega': self.omega_force, 'A': self.amplitude_force}
        return GenericCallFunction(self.excitation_function_force, scope_vars, check=False)

    def SolveSolutionStep(self):
        b = self.RHS_matrix @ self.x[:,0]
        #external load
        excitation_load = self.ApplyForceExcitation()
        self.load_vector[-1] += excitation_load
        # print("external load= ", self.load_vector[-1])
        b += self.load_vector
        #root point displacement
        d_f_excitation = self.ApplyRootPointExcitation()
        d_f = d_f_excitation + self.root_point_displacement
        b_f = self.CalculateEquivalentForceFromRootPointExcitation(d_f)
        b += b_f
        self.dx = np.linalg.solve(self.LHS, b)

    def CalculateReaction(self, buffer_idx=0):
        reaction = self.damping * ( self.dx[1] - self.dx_f[1]) \
                 + self.stiffness * (self.dx[0] - self.dx_f[0])
        return reaction

    def CalculateSelfWeight(self):
        self_weight = self.mass * self.modulus_self_weight
        return self_weight

    def GetSolutionStepValue(self, identifier, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            return self.x[:,buffer_idx][0]
        elif identifier == "VELOCITY":
            return self.x[:,buffer_idx][1]
        elif identifier == "ACCELERATION":
            return self.x[:,buffer_idx][2]
        elif identifier == "REACTION":
            return self.CalculateReaction()
        elif identifier == "VOLUME_ACCELERATION":
            return self.CalculateSelfWeight()
        else:
            raise Exception("Identifier is unknown!")

    def SetSolutionStepValue(self, identifier, value, buffer_idx=0):
        if identifier == "DISPLACEMENT":
            self.x[:,buffer_idx][0] = value
        elif identifier == "VELOCITY":
            self.x[:,buffer_idx][1] = value
        elif identifier == "ACCELERATION":
            self.x[:,buffer_idx][2] = value
        elif identifier == "LOAD":
            self.load_vector[-1] = 0.0
            self.load_vector[-1] = value
        elif identifier == "ROOT_POINT_DISPLACEMENT":
            self.root_point_displacement = 0.0
            self.root_point_displacement = value
        else:
            raise Exception("Identifier is unknown!")

    def _ValidateAndAssignRigidBodySolverDefaults(self, parameters):
        
        self._CheckMandatoryInputParameters(parameters)

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
                "path" : "rigid_body_solver"
            }
        }''')

        dof_parameters = {}
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

            single_dof_parameters = KratosMultiphysics.Parameters(json.dumps(single_dof_parameters))
            single_dof_parameters.RecursivelyValidateAndAssignDefaults(default_single_dof_parameters)

            active_dofs.append(dof)
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