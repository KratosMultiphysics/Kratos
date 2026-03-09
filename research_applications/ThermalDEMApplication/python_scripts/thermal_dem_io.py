# Imports
import os
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.ThermalDEMApplication import *

# Import correct DEM_procedures file (as it is done in DEM_analysis_stage for importing procedures)
if IsDistributedRun():
    if "DO_NOT_PARTITION_DOMAIN" in os.environ:
        from KratosMultiphysics.mpi import *
        import KratosMultiphysics.DEMApplication.DEM_procedures_mpi_no_partitions as DEM_procedures
    else:
        from KratosMultiphysics.MetisApplication import *
        from KratosMultiphysics.MPISearchApplication import *
        from KratosMultiphysics.mpi import *
        import KratosMultiphysics.DEMApplication.DEM_procedures_mpi as DEM_procedures
else:
    import KratosMultiphysics.DEMApplication.DEM_procedures as DEM_procedures

# Set base class
BaseIO = DEM_procedures.DEMIo

# Auxiliary functions
def GetBoolParameterIfItExists(parameters, key):
    if key in parameters.keys():
        return parameters[key].GetBool()
    else:
        return False

# IO class
class ThermalDEMIo(BaseIO):

    def __init__(self, model, DEM_parameters, post_path, all_model_parts):
        # Set standard pos options
        BaseIO.__init__(self, model, DEM_parameters, post_path, all_model_parts)

        # Set thermal pos options
        self.PostTemperature = GetBoolParameterIfItExists(self.DEM_parameters, "PostTemperature")
        self.PostHeatFlux    = GetBoolParameterIfItExists(self.DEM_parameters, "PostHeatFlux")

    def Initialize(self, DEM_parameters):
        # Add standard pos variables
        BaseIO.Initialize(self, DEM_parameters)

        # Add thermal pos variables
        self.PushPrintVar(self.PostTemperature, TEMPERATURE, self.global_variables)
        self.PushPrintVar(self.PostHeatFlux,    HEATFLUX,    self.global_variables)