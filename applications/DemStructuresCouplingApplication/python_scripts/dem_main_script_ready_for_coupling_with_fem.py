from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import KratosMultiphysics as kratos
import KratosMultiphysics.DEMApplication as Dem

from KratosMultiphysics.DEMApplication import DEM_analysis_stage

BaseAlgorithm = DEM_analysis_stage.Solution

class Solution(BaseAlgorithm):

    def __init__(self, model):
        with open("ProjectParametersDEM.json",'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        model = KratosMultiphysics.Model()
        super(Solution,self).__init__(model, project_parameters)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        self.coupling_algorithm.ReadDemModelParts()

    def BaseReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        super(Solution, self).ReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)

    def PrintResultsForGid(self, time):
        self.coupling_algorithm.gid_output.Writeresults(time)

if __name__ == "__main__":
    with open("ProjectParametersDEM.json",'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    Solution(model, project_parameters).Run()
