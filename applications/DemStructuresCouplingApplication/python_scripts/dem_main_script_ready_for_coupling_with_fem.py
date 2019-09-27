from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as Dem

from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage

class StructuresCoupledDEMAnalysisStage(DEMAnalysisStage):

    def __init__(self, model,parameters):
        super(StructuresCoupledDEMAnalysisStage,self).__init__(model, parameters)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        self.coupling_analysis.ReadDemModelParts()

    def BaseReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        super(StructuresCoupledDEMAnalysisStage, self).ReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def PrintResultsForGid(self, time):
        self.coupling_analysis.gid_output.Writeresults(time)

if __name__ == "__main__":
    parameter_file_name = "ProjectParametersDEM.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = StructuresCoupledDEMAnalysisStage(model, parameters)
    simulation.Run()
