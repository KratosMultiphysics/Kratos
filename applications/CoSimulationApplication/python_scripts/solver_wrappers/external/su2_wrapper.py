# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities.data_communicator_utilities import GetRankZeroDataCommunicator
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

import sys
import shutil
import pysu2
from mpi4py import MPI


def Create(settings, model, solver_name):
    return SU2Wrapper(settings, model, solver_name)

class SU2Wrapper(CoSimulationSolverWrapper):
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        flow_filename = settings["solver_wrapper_settings"]["input_file"].GetString()
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        # Initialize the flow driver of SU2, this includes solver preprocessing
        self.FlowDriver = pysu2.CSinglezoneDriver(flow_filename, 1, comm);
        FlowMarkerID = None
        FlowMarkerName = settings["solver_wrapper_settings"]["interface_marker"].GetString()                                            # FSI marker (flow side)
        FlowMarkerList = self.FlowDriver.GetAllBoundaryMarkersTag()              # Get all the flow boundary tags
        FlowMarkerIDs = self.FlowDriver.GetAllBoundaryMarkers()                  # Get all the associated indices to the flow markers
        if FlowMarkerName in FlowMarkerList and FlowMarkerName in FlowMarkerIDs.keys():
            FlowMarkerID = FlowMarkerIDs[FlowMarkerName]                      # Check if the flow FSI marker exists
        nVertex_Marker_Flow = self.FlowDriver.GetNumberVertices(FlowMarkerID)    # Get the number of vertices of the flow FSI marker

        model = KM.

    def SolveSolutionStep(self):
        print("\n------------------------------ Begin Solver -----------------------------\n")
        self.FlowDriver.ResetConvergence()
        self.FlowDriver.Preprocess(0)
        self.FlowDriver.Run()
        self.FlowDriver.Postprocess()
        stopCalc = self.FlowDriver.Monitor(0)