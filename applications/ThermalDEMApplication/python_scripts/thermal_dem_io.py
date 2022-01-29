import os
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.ThermalDEMApplication import *

# Import correct DEM_procedures file
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

class ThermalDEMIo(BaseIO):

    def __init__(self, model, DEM_parameters, post_path, all_model_parts):
        BaseIO.__init__(self, model, DEM_parameters, post_path, all_model_parts)
