from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.mpi import *

import DEM_procedures
import DEM_procedures_mpi

PostUtils = DEM_procedures.PostUtils

class Procedures(DEM_procedures_mpi.Procedures):
    def AddMpiVariables(self, model_part):
        pass

    def FindMaxNodeIdInModelPart(self, model_part):
        DEM_procedures.Procedures.FindMaxNodeIdInModelPart(model_part)

class MultifileList(DEM_procedures.MultifileList):
    def __init__(self, name, step):
        if (mpi.rank == 0):
            DEM_procedures.MultifileList.__init__(name,step)

DEMFEMProcedures = DEM_procedures_mpi.DEMFEMProcedures
Report = DEM_procedures.Report
MaterialTest = DEM_procedures.MaterialTest

class DEMIo(DEM_procedures.DEMIo):
    def PrintResults(self, all_model_parts, creator_destructor, dem_fem_search, time, bounding_box_time_limits):
        if (mpi.rank == 0):
            super(DEMIo,self).PrintResults(all_model_parts, creator_destructor, dem_fem_search, time, bounding_box_time_limits)

ParallelUtils = DEM_procedures.ParallelUtils
SetOfModelParts = DEM_procedures.SetOfModelParts
DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator