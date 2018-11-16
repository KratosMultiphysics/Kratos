from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as Kratos
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.MeshingApplication as Meshing
import KratosMultiphysics.ShallowWaterApplication as Shallow

from shallow_water_analysis import ShallowWaterAnalysis

class ShallowWaterMultigridAnalysis(ShallowWaterAnalysis):
    ''' Main script for shallow water simulations '''
    
