# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.CoSimulationApplication as Kratos_CoSim

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities.data_communicator_utilities import GetRankZeroDataCommunicator
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

import sys
import shutil
import pysu2
from mpi4py import MPI
import numpy as np

from KratosMultiphysics.CoSimulationApplication import CoSimIO


def Create(settings, model, solver_name):
    print(dir(CoSimIO))
    raise RuntimeError(1)
    return SU2Wrapper(settings, model, solver_name)

class SU2Wrapper(CoSimulationSolverWrapper):
    def __init__(self, settings, model: Kratos.Model, solver_name):
        super().__init__(settings, model, solver_name)

        flow_filename = settings["solver_wrapper_settings"]["input_file"].GetString()

        self.MPI = MPI
        self.comm = MPI.COMM_WORLD			#MPI World communicator
        self.myid = self.comm.Get_rank()

        # Initialize the flow driver of SU2, this includes solver preprocessing

        self.FlowDriver = pysu2.CSinglezoneDriver(flow_filename, 1, self.comm)
        self.FlowMarkerID = None
        FlowMarkerName = settings["solver_wrapper_settings"]["interface_marker"].GetString()                                            # FSI marker (flow side)
        FlowMarkerList = self.FlowDriver.GetAllBoundaryMarkersTag()              # Get all the flow boundary tags
        FlowMarkerIDs = self.FlowDriver.GetAllBoundaryMarkers()                  # Get all the associated indices to the flow markers
        if FlowMarkerName in FlowMarkerList and FlowMarkerName in FlowMarkerIDs.keys():
            self.FlowMarkerID = FlowMarkerIDs[FlowMarkerName]                      # Check if the flow FSI marker exists
            self.nVertex_Marker_Flow = self.FlowDriver.GetNumberVertices(self.FlowMarkerID)    # Get the number of vertices of the flow FSI marker
        else:
            self.nVertex_Marker_Flow = 0

        # Getting mesh information from SU2 and store it as an Model Part

        self.comm.barrier()
        local_nodal_coords = np.zeros([self.nVertex_Marker_Flow, 3], dtype=np.float64)
        # local_nodal_coords = np.array(self.nVertex_Marker_Flow)
        for node_i in range(self.nVertex_Marker_Flow):
            coords = self.FlowDriver.GetInitialMeshCoord(self.FlowMarkerID, node_i)
            local_nodal_coords[node_i, 0] = coords[0]
            local_nodal_coords[node_i, 1] = coords[1]
            local_nodal_coords[node_i, 2] = coords[2]
        
        # print(f"myid is {self.myid}, local_nodal_coords = {local_nodal_coords}")


        self.model_part = self.model.CreateModelPart(FlowMarkerName)
        self.model_part.AddNodalSolutionStepVariable(Kratos.MESH_DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        for node_i in range(self.nVertex_Marker_Flow):
            self.model_part.CreateNewNode(node_i + 1, local_nodal_coords[node_i, 0], local_nodal_coords[node_i, 1], local_nodal_coords[node_i, 2])        

        # model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        # model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

    def SolveSolutionStep(self):

        # for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
        #     data_config = {
        #         "type" : "coupling_interface_data",
        #         "interface_data" : self.GetInterfaceData(data_name)
        #     }
        #     self.ImportData(data_config)

        # Get Kratos variable with displacement and set it as SetMeshDisplacement
        displacements = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(displacements, Kratos.MESH_DISPLACEMENT, True)
        displacement_numpy = displacements.Evaluate()

        for node_i in range(self.nVertex_Marker_Flow):
            if self.comm.Get_rank() == 0:
                print("Is Halo", self.FlowDriver.IsAHaloNode(self.FlowMarkerID, node_i))
            # if self.FlowDriver.IsAHaloNode(self.FlowMarkerID, node_i):
            #     print("Is Halo", self.FlowDriver.IsAHaloNode(self.FlowMarkerID, node_i))
            self.FlowDriver.SetMeshDisplacement(self.FlowMarkerID, node_i, displacement_numpy[node_i][0], displacement_numpy[node_i][1], displacement_numpy[node_i][2])
        
        raise RuntimeError(1)

        print(f"myid is {self.myid}, displacement = {displacement_numpy}")

        if self.myid == 0:
            print("\n------------------------------ Begin Solver -----------------------------\n")
        self.FlowDriver.ResetConvergence()
        self.FlowDriver.Preprocess(0)
        self.FlowDriver.Run()
        self.FlowDriver.Postprocess()
        stopCalc = self.FlowDriver.Monitor(0)



        # Recover the flow loads from SU2 and assign it to corresponding Kratos Variable   
        flow_loads = np.zeros([self.nVertex_Marker_Flow, 3], dtype=np.float64)
        for node_i in range(self.nVertex_Marker_Flow):
            vertex_load = self.FlowDriver.GetFlowLoad(self.FlowMarkerID, node_i)
            flow_loads[node_i, 0] = vertex_load[0]
            flow_loads[node_i, 1] = vertex_load[1]
            flow_loads[node_i, 2] = vertex_load[2]     

        forces = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.CArrayExpressionIO.Read(forces, flow_loads)
        Kratos.Expression.VariableExpressionIO.Write(forces, Kratos.REACTION, True)

        temp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(temp, Kratos.REACTION, True)

        print(f"myid is {self.myid}, local_force = {np.array(flow_loads)}")
        print(f"myid is {self.myid}, local_force = {temp.Evaluate()}")

        # for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
        #     data_config = {
        #         "type" : "coupling_interface_data",
        #         "interface_data" : self.GetInterfaceData(data_name)
        #     }
        #     self.ExportData(data_config)

    def AdvanceInTime(self, current_time):
        return 0.0 # TODO find a better solution here... maybe get time from solver through IO
    
    # def _GetIOType(self):
    #     return self.settings["io_settings"]["type"].GetString()

    def ImportCouplingInterface(self, interface_config):
        pass
        #for node_i in range(self.nVertex_Marker_Flow):
