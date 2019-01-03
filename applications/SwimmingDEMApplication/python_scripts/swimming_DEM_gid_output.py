from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import gid_output


class SwimmingDEMGiDOutput(gid_output.GiDOutput):

    def __init__(self,
                 file_name,
                 vol_output=True,
                 post_mode="Binary",
                 multifile="Single",
                 deformed_mesh=False,
                 write_conditions=True
                 ):
        gid_output.GiDOutput.__init__(self,
                                      file_name,
                 vol_output,
                 post_mode,
                 multifile,
                 deformed_mesh,
                 write_conditions)

        self.outerlistfilename = os.path.join("..", self.listfilename)


    def initialize_swimming_DEM_results(self, DEM_model_part, clusters_model_part, rigid_faces_model_part, mixed_model_part):

        if self.multi_file == MultiFileFlag.SingleFile:
            print("Singlefile option is not available for the swimming DEM application!")
            mesh_name = 0.0
            self.io.InitializeMesh(mesh_name)
            self.io.WriteSphereMesh(DEM_model_part.GetMesh())
            self.io.WriteMesh(clusters_model_part.GetMesh())
            self.io.WriteMesh(rigid_faces_model_part.GetMesh())
            self.io.WriteMesh(mixed_model_part.GetMesh())
            self.io.FinalizeMesh()
            self.io.InitializeResults(mesh_name, mixed_model_part.GetMesh())

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


    def write_swimming_DEM_results(self, label,
                                  fluid_model_part,
                                  DEM_model_part,
                                  clusters_model_part,
                                  rigid_faces_model_part,
                                  mixed_model_part,
                                  fluid_nodal_variables,
                                  DEM_nodal_variables,
                                  cluster_variables,
                                  rigid_faces_variables,
                                  mixed_nodal_variables,
                                  fluid_gp_variables):
        # label = str(label) #it should be a C double
        # update cut data if necessary
        out_model_part = self.get_out_model_part(fluid_model_part)

        # update cut data if necessary
        if not self.volume_output:
            self.cut_app.UpdateCutData(out_model_part, fluid_model_part)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self.io.InitializeMesh(label)
            self.io.WriteSphereMesh(DEM_model_part.GetMesh())
            self.io.WriteMesh(mixed_model_part.GetMesh())
            self.io.WriteMesh(rigid_faces_model_part.GetMesh())
            self.io.FinalizeMesh()
            self.io.InitializeResults(label, mixed_model_part.GetMesh())

        for var in fluid_nodal_variables:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, fluid_model_part, kratos_variable)

        for var in DEM_nodal_variables:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, DEM_model_part, kratos_variable)

        for var in cluster_variables:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, clusters_model_part, kratos_variable)

        for var in rigid_faces_variables:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, rigid_faces_model_part, kratos_variable)

        for var in mixed_nodal_variables:
            kratos_variable = globals()[var]
            self._write_nodal_results(label, mixed_model_part, kratos_variable)

        for var in fluid_gp_variables:
            kratos_variable = globals()[var]
            self._write_gp_results(label, fluid_model_part, kratos_variable)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self._finalize_results()

            with open(self.listfilename, "a") as listfile:
                self.write_step_to_list(label)

            self.write_step_to_outer_list(label)

