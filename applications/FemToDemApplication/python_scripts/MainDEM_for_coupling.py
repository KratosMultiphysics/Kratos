from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics.DEMApplication
import KratosMultiphysics.DEMApplication.main_script as MainDEM
import KratosMultiphysics.FemToDemApplication as FEMDEM

if IsDistributedRun():
    if "DO_NOT_PARTITION_DOMAIN" in os.environ:
        Logger.PrintInfo("DEM", "Running under MPI........")
        from KratosMultiphysics.mpi import *
        import KratosMultiphysics.DEMApplication.DEM_procedures_mpi_no_partitions as DEM_procedures
    else:
        Logger.PrintInfo("DEM", "Running under OpenMP........")
        from KratosMultiphysics.MetisApplication import *
        from KratosMultiphysics.MPISearchApplication import *
        from KratosMultiphysics.mpi import *
        import KratosMultiphysics.DEMApplication.DEM_procedures_mpi as DEM_procedures
else:
    Logger.PrintInfo("DEM", "Running under OpenMP........")
    import KratosMultiphysics.FemToDemApplication.DEM_procedures_for_coupling as DEM_procedures

class DEM_for_coupling_Solution(MainDEM.Solution):

    def Info(self):
        print("DEM part of the FEM-DEM application")

    def SetAnalyticParticleWatcher(self):
        pass
