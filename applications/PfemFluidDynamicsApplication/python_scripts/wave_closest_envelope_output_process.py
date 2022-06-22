import KratosMultiphysics as KM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from KratosMultiphysics.PfemFluidDynamicsApplication.wave_envelope_output_process import WaveEnvelopeOutputProcess
from numpy import linspace

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WaveClosestEnvelopeOutputProcess(model, settings["Parameters"])

class WaveClosestEnvelopeOutputProcess(WaveEnvelopeOutputProcess):
    """WaveClosestEnvelopeOutputProcess

    This process records the wave height along a path.
    If a sampling point is not found, Nan will be printed.
    Possible specifications of the Parameters:
     - output_file_settings: a parameters encapsulating the 'file_name', 'output_path' and
                             other settings according to 'TimeBasedAsciiFileWritterUtility'.
     - wave_calculation_settings: a parameters according to 'CalculateWaveClosestHeightOutputProcess'
    """

    def ExecuteBeforeSolutionLoop(self):
        """Initialize the utility to calculate the water height"""
        wave_settings = self.settings["wave_calculation_settings"]
        ## this gives the average height of the closest nodes to the gauge over a user-provided tolerance distance
        self.wave_height_utility = PFEM.CalculateWaveClosestHeightUtility(self.model_part, wave_settings)
        self.max_values = [-1e6] * len(self.positions)
