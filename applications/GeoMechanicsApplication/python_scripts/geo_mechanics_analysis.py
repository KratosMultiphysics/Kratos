from KratosMultiphysics.GeoMechanicsApplication import geomechanics_analysis

class GeoMechanicsAnalysis(geomechanics_analysis.GeoMechanicsAnalysis):

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)
