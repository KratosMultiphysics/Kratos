from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressFlow

# Importing the problem analysis stage class
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis

# Avoid printing of Kratos informations
# KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

# Import python libraries
from math import *

# TODO: use a r.v. for Mach and one for alpha

"""
SimulationScenario is inherited from the PotentialFlowAnalysis class and solves the potential flow
-lapl(u) = \varepsilon*f    u \in \Omega
                            u \in \partial(\Omega)
where
\varepsilon \sim
and
f =
the QoI is
lift coeficient
"""
class SimulationScenario(PotentialFlowAnalysis):
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(SimulationScenario,self).__init__(input_model,input_parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    """
    function introducing the stochasticity in the right hand side
    input:  self: an instance of the class
    """
    def ModifyInitialProperties(self):
        '''Introduce here the stochasticity in the Mach number and the angle of attack'''
        Mach = self.sample[0]
        a_infinity = 340 # [m/s] velocity of sound at infinity
        alpha =  self.sample[1]
        v_norm = Mach * a_infinity
        velocity = [v_norm*cos(alpha),v_norm*sin(alpha),0]
        boundary_processes = self.project_parameters["processes"]["boundary_conditions_process_list"]
        problem_name=self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
        for i in range(0,boundary_processes.size()):
            python_module = boundary_processes[i]["python_module"].GetString()
            if python_module == "apply_far_field_process":
                self.project_parameters["processes"]["boundary_conditions_process_list"][i]["Parameters"]["velocity_infinity"].SetVector(velocity)
        # self.project_parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(problem_name+'_M'+str(Mach)+'_A'+str(alpha))

    """
    function evaluating the QoI of the problem: lift coefficient
    input:  self: an instance of the class
    """
    def EvaluateQuantityOfInterest(self):
        Q = self._GetSolver().main_model_part.GetValue(KratosMultiphysics.FRICTION_COEFFICIENT)
        return Q
