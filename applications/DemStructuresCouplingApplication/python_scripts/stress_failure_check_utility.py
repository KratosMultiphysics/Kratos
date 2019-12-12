import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

class StressFailureCheckUtility(object):
    def __init__(self, model_part, test_number):

        if not test_number:
            return

        self.model_part = model_part

        if test_number == 1: # CTW16
            minimum_radius = 0.00381
            maximum_radius = 0.01381 #0.00501
        elif test_number == 2: # CTW10
            minimum_radius = 0.0127
            maximum_radius = 0.0139
        else: # Blind test
            minimum_radius = 0.0381
            maximum_radius = 0.0393

        self.settings = KratosMultiphysics.Parameters( """
        {
            "cylinder_center": [0.0,0.0,0.0]
        }  """ )

        self.settings.AddEmptyValue("min_radius")
        self.settings["min_radius"].SetDouble(minimum_radius)
        self.settings.AddEmptyValue("max_radius")
        self.settings["max_radius"].SetDouble(maximum_radius)

        self.utility = DemFem.StressFailureCheckUtilities(self.model_part, self.settings)

    def ExecuteFinalizeSolutionStep(self):

        self.utility.ExecuteFinalizeSolutionStep()