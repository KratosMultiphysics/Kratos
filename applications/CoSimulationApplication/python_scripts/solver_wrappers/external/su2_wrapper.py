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



def Create(settings, model, solver_name):
    return SU2Wrapper(settings, model, solver_name)

class SU2Wrapper(CoSimulationSolverWrapper):
    def __init__(self, settings, model: Kratos.Model, solver_name):
        super().__init__(settings, model, solver_name)

        flow_filename = settings["solver_wrapper_settings"]["input_file"].GetString()

        self.MPI = MPI
        self.comm = MPI.COMM_WORLD
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

        local_nodal_coords = []
        global_id_is_local_node_list = []
        # local_nodal_coords = np.array(self.nVertex_Marker_Flow)
        for node_i in range(self.nVertex_Marker_Flow):
            coords = self.FlowDriver.GetInitialMeshCoord(self.FlowMarkerID, node_i)
            local_nodal_coords.append(coords)
            is_local = not self.FlowDriver.IsAHaloNode(self.FlowMarkerID, node_i)
            global_id = self.FlowDriver.GetVertexGlobalIndex(self.FlowMarkerID, node_i)
            global_id_is_local_node_list.append([global_id, is_local])

        self.model_part = self.model.CreateModelPart(FlowMarkerName)
        self.model_part.AddNodalSolutionStepVariable(Kratos.MESH_DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        Kratos_CoSim.ModelPartIO.FillCoSimIOModelPart(self.model_part, global_id_is_local_node_list, local_nodal_coords, Kratos.ParallelEnvironment.GetDefaultDataCommunicator())

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

        for node_i, node in enumerate(self.model_part.Nodes):
            mesh_displacement = node.GetSolutionStepValue(Kratos.MESH_DISPLACEMENT)
            self.FlowDriver.SetMeshDisplacement(self.FlowMarkerID, node_i, mesh_displacement[0], mesh_displacement[1], mesh_displacement[2])
        
        print(f"displacement = {displacement_numpy}")

        if self.myid == 0:
            print("\n------------------------------ Begin Solver -----------------------------\n")
        self.FlowDriver.ResetConvergence()
        self.FlowDriver.Preprocess(0)
        self.FlowDriver.Run()
        self.FlowDriver.Postprocess()
        stopCalc = self.FlowDriver.Monitor(0)



        # Recover the flow loads from SU2 and assign it to corresponding Kratos Variable   
        flow_loads = np.zeros([self.nVertex_Marker_Flow, 3], dtype=np.float64)
        for node_i, node in enumerate(self.model_part.Nodes):
            vertex_load = self.FlowDriver.GetFlowLoad(self.FlowMarkerID, node_i)
            node.SetSolutionStepValue(Kratos.REACTION, 0, Kratos.Array3(list(vertex_load)))
            
        forces = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(forces, Kratos.REACTION, True)

        print(f"local_force = {forces.Evaluate()}")

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
