import sys
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.DEMApplication as Dem
import KratosMultiphysics.StructuralMechanicsApplication as Structural
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

from KratosMultiphysics.DemStructuresCouplingApplication.dem_fem_coupling_algorithm import Algorithm

Algorithm().Run()
