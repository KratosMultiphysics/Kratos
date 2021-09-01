# Importing the Kratos Library
import KratosMultiphysics as KM

# debug
import pdb

# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Importing PfemFluidDynamics
if not CheckIfApplicationsAvailable("PfemFluidDynamicsApplication"):
    raise ImportError("The PfemFluidDynamicsApplication is not available!")
from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis import PfemFluidDynamicsAnalysis

def Create(settings, model, solver_name):
    return PfemFluidDynamicsWrapper(settings, model, solver_name)

class PfemFluidDynamicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the PfemFluidDynamicsApplication of Kratos"""

    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)
        self.repeat_time_step_flag = False #At initial timestep should be false

    def _CreateAnalysisStage(self):
        return PfemFluidDynamicsAnalysis(self.model, self.project_parameters)
    
    '''def Initialize(self):
        super().Initialize()
        self.export_config = {
            "type" : "repeat_time_step",
            "repeat_time_step" : False
        }
        # self.repeat_time_step = self.export_config["repeat_time_step"]
        # save nodes in model parts which need to be moved while simulating
        self.list_of_nodes_in_move_mesh_model_parts = [self.model[mp_name].Nodes for mp_name in self.settings["solver_wrapper_settings"]["move_mesh_model_part"].GetStringArray()]'''
    
    def SolveSolutionStep(self):
        # move the pfem mesh fluid back to original state?
        #TODO change to only do when repeating a timestep.
        # pdb.set_trace()
        # input_output = self.__io.repeat_time_step_flag
        #input_output = self._KratosBaseWrapper_CoSimulationSolverWrapper__io
        #input_output = self.__GetIO().repeat_time_step_flag
        #if self.__GetIO().repeat_time_step_flag:

        #repeat_time_step_flag = self.ImportData(dummy)
        
        # if self.repeat_time_step_flag:
            # pdb.set_trace()
            #TODO fix for the computing model part. If we have the SOLID and RIGID check this shouldn't make any difference
            # KM.PfemFluidDynamicsApplication.MoveMeshUtility().ResetPfemKinematicValues(self._analysis_stage._GetSolver().model["PfemFluidModelPart.Fluid"])

        # solve PFEM
        super().SolveSolutionStep()

    def _GetIOType(self):
        # only external solvers have to specify sth here / override this
        return "pfem_io"

    def __GetIO(self):
        super().__GetIO()
