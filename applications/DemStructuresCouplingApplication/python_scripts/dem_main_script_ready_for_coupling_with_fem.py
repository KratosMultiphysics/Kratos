from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import os
import KratosMultiphysics as kratos
import KratosMultiphysics.DEMApplication as Dem

import main_script

BaseAlgorithm = main_script.Solution

class Solution(BaseAlgorithm):

    def __init__(self, model):
        super(Solution, self).__init__(model)

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        self.coupling_algorithm.ReadDemModelParts()

    def BaseReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        #super(Solution, self).ReadModelParts(max_node_Id, max_elem_Id, max_cond_Id)
        os.chdir(self.main_path)

        # Reading the model_part
        spheres_mp_filename = self.GetMpFilename()
        model_part_io_spheres = self.model_part_reader(spheres_mp_filename, max_node_Id, max_elem_Id, max_cond_Id)

        if "do_not_perform_initial_partition" in self.DEM_parameters.keys() and self.DEM_parameters["do_not_perform_initial_partition"].GetBool():
            pass
        else:
            self.parallelutils.PerformInitialPartition(model_part_io_spheres)

        os.chdir(self.main_path)
        [model_part_io_spheres, self.spheres_model_part, MPICommSetup] = self.parallelutils.SetCommunicator(self.spheres_model_part, model_part_io_spheres, spheres_mp_filename)
        model_part_io_spheres.ReadModelPart(self.spheres_model_part)

        max_node_Id += self.creator_destructor.FindMaxNodeIdInModelPart(self.spheres_model_part)
        max_elem_Id += self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part)
        old_max_elem_Id_spheres = max_elem_Id
        max_cond_Id += self.creator_destructor.FindMaxConditionIdInModelPart(self.spheres_model_part)

        #removed part reading the Rigid Faces

        clusters_mp_filename = self.GetClusterFilename()
        model_part_io_clusters = self.model_part_reader(clusters_mp_filename, max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)
        model_part_io_clusters.ReadModelPart(self.cluster_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part)
        if max_elem_Id != old_max_elem_Id_spheres:
            self.creator_destructor.RenumberElementIdsFromGivenValue(self.cluster_model_part, max_elem_Id)

        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(self.cluster_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.cluster_model_part)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(self.cluster_model_part)
        DEM_Inlet_filename = self.GetInletFilename()
        model_part_io_demInlet = self.model_part_reader(DEM_Inlet_filename, max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)
        model_part_io_demInlet.ReadModelPart(self.DEM_inlet_model_part)

        self.model_parts_have_been_read = True
        self.all_model_parts.ComputeMaxIds()

    def PrintResultsForGid(self, time):
        self.coupling_algorithm.gid_output.Writeresults(time)

if __name__ == "__main__":
    model = Model()
    Solution(model).Run()
