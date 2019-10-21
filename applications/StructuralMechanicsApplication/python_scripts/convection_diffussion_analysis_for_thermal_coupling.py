from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion as solver_wrapper

# Importing the base class
import KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis as convection_diffusion_analysis

class ConvectionDiffusionAnalysisForThermalCoupling(convection_diffusion_analysis.ConvectionDiffusionAnalysis):
    """
    This class is the main-script of the ConvectionDiffusionApplication put in a class
    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        # Making sure that older cases still work by properly initalizing the parameters
        super(ConvectionDiffusionAnalysisForThermalCoupling, self).__init__(model, project_parameters)