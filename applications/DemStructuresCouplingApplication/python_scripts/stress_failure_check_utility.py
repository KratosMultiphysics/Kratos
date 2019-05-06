import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

class StressFailureCheckUtility(object):
    def __init__(self, model_part):

        self.model_part = model_part
        self.settings = KratosMultiphysics.Parameters( """
        {
            "cylinder_center": [0.0,0.0,0.0],
            "min_radius": 0.00381,
            "max_radius": 0.00481
        }  """ )

        self.utility = DemFem.StressFailureCheckUtilities(self.model_part,self.settings)

    def ExecuteFinalizeSolutionStep(self):

        self.utility.ExecuteFinalizeSolutionStep()