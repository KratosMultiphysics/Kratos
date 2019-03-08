from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.EmpireApplication

# Importing the base class
from co_simulation_solvers.co_simulation_base_solver import CoSimulationBaseSolver

# Other imports
from co_simulation_tools import csprint, yellow

def CreateSolver(cosim_solver_settings, level):
    return KratosEmpireSolver(cosim_solver_settings, level)


class KratosEmpireSolver(CoSimulationBaseSolver):
    def __init__(self, cosim_solver_settings, level):
        super(KratosEmpireSolver, self).__init__(cosim_solver_settings, level)
        # TODO make client-model-part-name settable from outside?
        # TODO set Domain_size! (maybe from outside?) => this would also solve the quad/tetra problem in the empire-wrapper
        self.client_model_part = KratosMultiphysics.ModelPart("ClientModelPart")
        self.client_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)

        for var_name in self.cosim_solver_settings["nodal_variables"]:
            variable = KratosMultiphysics.KratosGlobals.GetVariable(var_name)
            self.client_model_part.AddNodalSolutionStepVariable(variable)
        self.client_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        self.client_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        self.model = KratosMultiphysics.Model()
        self.model.AddModelPart(self.client_model_part)

        self.xml_file_name = self.cosim_solver_settings["xml_file_name"]
        if not self.xml_file_name.endswith(".xml"):
            self.xml_file_name += ".xml"

        self.time_step = self.cosim_solver_settings["time_step"] # dummy, not really used


    def Initialize(self):
        csprint(self.lvl, yellow(self.__class__.__name__ + ":") + " Starting to initialize Empire")
        import empire_wrapper
        csprint(self.lvl, yellow(self.__class__.__name__ + ":") + " Wrapper-Import Successful")
        self.empire = empire_wrapper.EmpireWrapper()
        csprint(self.lvl, yellow(self.__class__.__name__ + ":") + " Wrapper Created")
        self.empire.Connect(self.xml_file_name)

        self.empire.ReceiveMesh("client_mesh", self.client_model_part)

    def Finalize(self):
        self.empire.Disconnect()

    def AdvanceInTime(self, current_time):
        new_time = current_time + self.time_step
        self.client_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.client_model_part.CloneTimeStep(new_time)
        self.step_is_repeated = False
        return new_time

    def ImportData(self, data_name, from_client):
        '''This function first calls the BaseClass, i.e. the IO to
        do sth with the data, e.g. mapping
        Then it calls Empire to send the data (to the connected client)
        '''
        # the following is needed for iterative coupling to tell the client
        # to repeat the same timestep
        if self.step_is_repeated:
            self.empire.SendConvergenceSignal(0) # 0 means that it is NOT converged
        self.step_is_repeated = True

        # this brings the data in the format of the client
        super(KratosEmpireSolver, self).ImportData(data_name, from_client)

        # send data to Client
        # TODO read params from input
        print("Before Sending")
        self.empire.SendDataField("client_mesh",
                                  "displacements",
                                  KratosMultiphysics.DISPLACEMENT)
        print("After Sending")


    def ExportData(self, data_name, to_client):
        '''This function first receives the data from Empire (from the connected client)
        Then it calls the BaseClass, i.e. the IO to do sth with
        the data, e.g. mapping
        '''
        # receive data from Client
        # TODO read params from input
        print("Before Receiving")
        self.empire.ReceiveDataField("client_mesh",
                                     "forces",
                                     KratosMultiphysics.DISPLACEMENT)
        print("After Receiving")

        # this brings the data in the format of the coupled-solver
        super(KratosEmpireSolver, self).ExportData(data_name, to_client)


    def FinalizeSolutionStep(self):
        '''If this function is called then it means that either convergence is achieved
        or the maximum number of iterations is achieved (in case of iterative coupling)
        In either way and also for explicit coupling, here the connected client is
        informed to advance in time
        '''
        self.empire.SendConvergenceSignal(1) # 1 means that it IS converged