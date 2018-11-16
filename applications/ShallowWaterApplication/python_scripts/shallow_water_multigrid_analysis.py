from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as Kratos
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.MeshingApplication as Meshing
import KratosMultiphysics.ShallowWaterApplication as Shallow

from shallow_water_analysis import ShallowWaterAnalysis
# from multiscale_refining_process import MultiscaleRefiningProcess

class ShallowWaterMultigridAnalysis(ShallowWaterAnalysis):
    ''' Main script for shallow water simulations '''

    def _GetOrderOfProcessesInitialization(self):
        return ["multigrid_process_list",
                "bathymetry_process_list",
                "initial_conditions_process_list",
                "boundary_conditions_process_list"]
