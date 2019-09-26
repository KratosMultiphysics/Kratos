from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.DemStructuresCouplingApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import gid_output


class DemStructuresCouplingGiDOutput(gid_output.GiDOutput):

    def __init__(self,
                 file_name,
                 vol_output,
                 post_mode,
                 multifile,
                 deformed_mesh,
                 write_conditions,
                 structures_model_part,
                 balls_model_part,
                 clusters_model_part,
                 rigid_faces_model_part,
                 contact_model_part,
                 mixed_model_part
                 ):
        gid_output.GiDOutput.__init__(self,
                                        file_name,
                                        vol_output,
                                        post_mode,
                                        multifile,
                                        deformed_mesh,
                                        write_conditions)

        self.GiDMultiFileFlag = multifile

        self.outerlistfilename = os.path.join("..", self.listfilename)

        self.structures_model_part = structures_model_part
        self.balls_model_part = balls_model_part
        self.clusters_model_part = clusters_model_part
        self.rigid_faces_model_part = rigid_faces_model_part
        self.contact_model_part = contact_model_part
        self.mixed_model_part = mixed_model_part

        self.structures_nodal_results = []
        self.dem_nodal_results = []
        self.clusters_nodal_results = []
        self.rigid_faces_nodal_results = []
        self.contact_gauss_points_results = []
        self.mixed_nodal_results = []
        self.structures_gauss_points_results = []

    def initialize_dem_fem_results(self,
                                    structures_nodal_results,
                                    dem_nodal_results,
                                    clusters_nodal_results,
                                    rigid_faces_nodal_results,
                                    contact_gauss_points_results,
                                    mixed_nodal_results,
                                    structures_gauss_points_results):
        self.structures_nodal_results = structures_nodal_results
        self.dem_nodal_results = dem_nodal_results
        self.clusters_nodal_results = clusters_nodal_results
        self.rigid_faces_nodal_results = rigid_faces_nodal_results
        self.contact_gauss_points_results = contact_gauss_points_results
        self.mixed_nodal_results = mixed_nodal_results
        self.structures_gauss_points_results = structures_gauss_points_results

        if self.multi_file == MultiFileFlag.SingleFile:
            print("Singlefile option is not available for the DEM-Structures Coupling application!")
            mesh_name = 0.0
            self.io.InitializeMesh(mesh_name)
            self.io.WriteSphereMesh(self.balls_model_part.GetMesh())
            self.io.WriteMesh(self.clusters_model_part.GetMesh())
            self.io.WriteMesh(self.rigid_faces_model_part.GetMesh())
            self.io.WriteMesh(self.contact_model_part.GetMesh())
            self.io.WriteMesh(self.mixed_model_part.GetMesh())
            self.io.FinalizeMesh()
            self.io.InitializeResults(mesh_name, self.mixed_model_part.GetMesh())

        # Initialize list file
        with open(self.listfilename, "w") as listfile:

            if self.multi_file == MultiFileFlag.MultipleFiles:
                listfile.write("Multiple\n")

            elif self.multi_file == MultiFileFlag.SingleFile:
                listfile.write("Single\n")

        if self.multi_file == MultiFileFlag.SingleFile:

            if self.post_mode == GiDPostMode.GiD_PostBinary:
                self.write_step_to_list()

            else:
                self.write_step_to_list(0)


    def write_step_to_list(self, step_label):

        if self.post_mode   == GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"

        with open(self.listfilename, "a") as listfile:
            listfile.write("Multiple\n")
            listfile.write(self.filename+"_"+"%.12g"%step_label+ext+"\n")

    def write_step_to_outer_list(self, step_label):

        if self.post_mode   == GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"

        with open(self.outerlistfilename, "a") as listfile:
            listfile.write("Multiple\n")
            folder_name = self.filename + "_Post_Files"
            full_string_to_write = os.path.join(folder_name,self.filename+"_"+"%.12g"%step_label+ext)
            listfile.write(full_string_to_write+"\n")


    def Writeresults(self, time):

        Logger.PrintInfo("DEM-Struct","")
        Logger.PrintInfo("DEM-Struct","*******************  PRINTING RESULTS FOR GID  ***************************")
        Logger.Flush()

        if self.GiDMultiFileFlag == "Multiples":
            self.mixed_model_part.Elements.clear()
            self.mixed_model_part.Nodes.clear()
            # here order is important!
            PostUtilities().AddModelPartToModelPart(self.mixed_model_part, self.balls_model_part)
            PostUtilities().AddModelPartToModelPart(self.mixed_model_part, self.rigid_faces_model_part)
            PostUtilities().AddModelPartToModelPart(self.mixed_model_part, self.contact_model_part)
            PostUtilities().AddModelPartToModelPart(self.mixed_model_part, self.structures_model_part)

        self.write_dem_fem_results(time)


    def write_dem_fem_results(self, label):
        # label = str(label) #it should be a C double
        # update cut data if necessary
        out_model_part = self.get_out_model_part(self.structures_model_part)

        # update cut data if necessary
        if not self.volume_output:
            self.cut_app.UpdateCutData(out_model_part, self.structures_model_part)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self.io.InitializeMesh(label)
            self.io.WriteSphereMesh(self.balls_model_part.GetMesh())
            self.io.WriteMesh(self.mixed_model_part.GetMesh())
            self.io.WriteMesh(self.rigid_faces_model_part.GetMesh())
            self.io.WriteMesh(self.contact_model_part.GetMesh())
            self.io.FinalizeMesh()
            self.io.InitializeResults(label, self.mixed_model_part.GetMesh())

        for var in  self.structures_nodal_results:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, self.structures_model_part, kratos_variable)

        for var in self.dem_nodal_results:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, self.balls_model_part, kratos_variable)

        for var in self.clusters_nodal_results:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, self.clusters_model_part, kratos_variable)

        for var in self.rigid_faces_nodal_results:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, self.rigid_faces_model_part, kratos_variable)

        for var in self.contact_gauss_points_results:
            kratos_variable = globals()[var]
            self._write_gp_results(label, self.contact_model_part, kratos_variable)

        for var in self.mixed_nodal_results:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, self.mixed_model_part, kratos_variable)

        for var in self.structures_gauss_points_results:
            kratos_variable = globals()[var]
            self._write_gp_results(label, self.structures_model_part, kratos_variable)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self._finalize_results()

            with open(self.listfilename, "a") as listfile:
                self.write_step_to_list(label)

            self.write_step_to_outer_list(label)

