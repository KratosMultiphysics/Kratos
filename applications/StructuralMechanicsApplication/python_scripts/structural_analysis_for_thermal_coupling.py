from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis as structural_mechanics_analysis

class StructuralMechanicsAnalysisForThermalCoupling(structural_mechanics_analysis.StructuralMechanicsAnalysis):
    """
    This class is the main-script of the StructuralMechanicsApplication put in a class

    It can be imported and used as "black-box"
    """
    def __init__(self, model, project_parameters):
        super(StructuralMechanicsAnalysisForThermalCoupling, self).__init__(model, project_parameters)