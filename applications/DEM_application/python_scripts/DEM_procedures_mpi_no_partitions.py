from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.mpi import *

import DEM_procedures
import DEM_procedures_mpi

class PostUtils(DEM_procedures.PostUtils):
    pass

class Procedures(DEM_procedures_mpi.Procedures):
    def AddMpiVariables(self, model_part):
        pass

    def FindMaxNodeIdInModelPart(self, model_part):
        DEM_procedures.Procedures.FindMaxNodeIdInModelPart(model_part)

class MultifileList(object):
    def __init__(self,name,step):
        if (mpi.rank == 0):
            DEM_procedures.MultifileList.__init__(name,step)

class DEMFEMProcedures(DEM_procedures_mpi.DEMFEMProcedures):
    pass

class Report(DEM_procedures.Report):
    pass
    '''
    def __init__(self):
        super(Report,self).__init__()
        '''

class MaterialTest(DEM_procedures.MaterialTest):
    pass
    '''
    def __init__(self):
        super(MaterialTest,self).__init__()
        '''

class DEMIo(DEM_procedures.DEMIo):
    def PrintResults(self, all_model_parts, creator_destructor, dem_fem_search, time, bounding_box_time_limits):
        if (mpi.rank == 0):
            super(DEMIo,self).PrintResults(all_model_parts, creator_destructor, dem_fem_search, time, bounding_box_time_limits)

class ParallelUtils(DEM_procedures.ParallelUtils):
    pass
    '''
    def __init__(self):
        super(ParallelUtils,self).__init__()
        '''
class SetOfModelParts(DEM_procedures.SetOfModelParts):
    pass

class DEMEnergyCalculator(DEM_procedures.DEMEnergyCalculator):
    pass