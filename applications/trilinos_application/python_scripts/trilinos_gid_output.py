from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.TrilinosApplication import *
CheckForPreviousImport()
import gid_output


class TrilinosGiDOutput(gid_output.GiDOutput):

    def write_step_to_list(self, step_label=None):
        if mpi.rank == 0:
            if self.post_mode == GiDPostMode.GiD_PostBinary:
                ext = ".post.bin"
            elif self.post_mode == GiDPostMode.GiD_PostAscii:
                ext = ".post.res"
            elif self.post_mode == GiDPostMode.GiD_PostAsciiZipped:
                ext = ".post.res"  # ??? CHECK!

            with open(self.listfilename, "a") as listfile:
                for r in range(0, mpi.size):
                    if step_label is not None:
                        listfile.write("{0}_{1}_{2}{3}\n".format(self.filename, r, step_label, ext))
                    else:
                        listfile.write("{0}_{1}{2}\n".format(self.filename, r, ext))

        mpi.world.barrier()

    def __init__(self,
                 file_name,
                 vol_output=True,
                 post_mode="Binary",
                 multifile="Single",
                 deformed_mesh=False,
                 write_conditions=True
                 ):
        # the GidIO C++ class needs the "partitioned" filename,
        # ending in the rank id. Passing it through the constructor
        # of the base class
        partitioned_file_name = "{0}_{1}".format(file_name, mpi.rank)
        gid_output.GiDOutput.__init__(self,
                                      partitioned_file_name,
                                      vol_output,
                                      post_mode,
                                      multifile,
                                      deformed_mesh,
                                      write_conditions)
        self.filename = file_name

    def define_cuts(self, model_part, cut_list):
        if self.volume_output:
            raise Exception("GiD IO Error: define_cuts was called for a problem set up to print volume output.")
        elif self.cut_number > 0:
            raise Exception("GiD IO Error: define_cuts was called twice. Please provide all cuts in a single list.")

        if self.cut_model_part is None:
            self.cut_model_part = ModelPart("CutPart")
            self.comm = CreateCommunicator()
            self.cut_app = TrilinosCuttingApplication(self.comm)
            self.cut_app.FindSmallestEdge(model_part)
            self.cut_app.AddVariablesToCutModelPart(model_part,self.cut_model_part)

        for cut_data in cut_list:
            self.cut_number += 1
            n = Vector(cut_data[0])
            p = Vector(cut_data[1])
            self.cut_app.GenerateCut(model_part, self.cut_model_part, n, p, self.cut_number, 0.01)
            mpi.world.barrier()

        self.cut_number += 1
        self.cut_app.AddSkinConditions(model_part, self.cut_model_part, self.cut_number)
        mpi.world.barrier()

    def _write_nodal_results(self, label, model_part, variable):
        self.io.WriteNodalResults(variable, model_part.GetCommunicator().LocalMesh().Nodes, label, 0)

    def initialize_results(self, model_part):

        if self.multi_file == MultiFileFlag.SingleFile:
            out_model_part = self.get_out_model_part(model_part)
            mesh_name = 0.0

            self._write_mesh(mesh_name, out_model_part)
            self._initialize_results(mesh_name, out_model_part)

        # Initialize list file
        if mpi.rank == 0:
            with open(self.listfilename, "w") as listfile:
                listfile.write("Merge\n")
                if self.multi_file == MultiFileFlag.MultipleFiles:
                    msg = """WARNING: In MPI, printing results in Multiple files
                             is not supported by the .post.lst list file.
                             Results will have to be open manually from GiD.\n"""
                    print (msg)

        if self.multi_file == MultiFileFlag.SingleFile:
            if self.post_mode == GiDPostMode.GiD_PostBinary:
                self.write_step_to_list()
            else:
                self.write_step_to_list(0)
