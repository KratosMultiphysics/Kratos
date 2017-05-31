from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import main_script

sys.path.insert(0,'')
import DEM_explicit_solver_var as DEM_parameters

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ or "I_MPI_INFO_NUMA_NODE_NUM" in os.environ:
    print("Running under MPI...........")
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
    def model_part_reader(modelpart, nodeid=0, elemid=0, condid=0):
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid)

else:
    print("Running under OpenMP........")
    import DEM_procedures
    import DEM_material_test_script
    #def model_part_reader(modelpart, a=0, b=0, c=0):
    #    return ModelPartIO(modelpart)
    def model_part_reader(modelpart, nodeid=0, elemid=0, condid=0):
        #return ModelPartIO(modelpart)
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid)


class Solution(main_script.Solution):        
    
    def SetGraphicalOutput(self):
        pass
        
    def GraphicalOutputInitialize(self):
        pass
        
    def PrintResultsForGid(self, time):
        pass
        
    def GraphicalOutputFinalize(self):
        pass