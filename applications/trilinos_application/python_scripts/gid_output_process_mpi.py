from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
CheckForPreviousImport()

import gid_output_process

class GiDOutputProcessMPI(gid_output_process.GiDOutputProcess):

    def __init__(self,model_part,file_name,param = None):

        if param is None:
            param = self.defaults
        else:
            # Note: this only validadtes the first level of the JSON tree.
            # I'm not going for recursive validation because some branches may
            # not exist and I don't want the validator assinging defaults there.
            param.ValidateAndAssignDefaults(self.defaults)

        self.param = param
        self.base_file_name = file_name + "_" + str(mpi.rank)

        self.model_part = model_part
        self.body_io = None
        self.volume_list_files = []

        # The following are only used if we asked to print results on surfaces
        self.cut_model_part = None
        self.cut_io = None
        self.output_surface_index = 0
        self.cut_list_files = []

        point_data_configuration = self.param["point_data_configuration"]
        if point_data_configuration.size() > 0:
            import point_output_process
            self.point_output_process = point_output_process.PointOutputProcess(self.model_part,point_data_configuration)
        else:
            self.point_output_process = None

        self.step_count = 0
        self.printed_step_count = 0
        self.next_output = 0.0

    #~ def __initialize_list_files(self,additional_frequencies):
        #~ '''Set up .post.lst files for global and cut results.
        #~ If we have only one tipe of output (volume or cut), the
        #~ list file is called <gid_model_name>.post.lst. When we have
        #~ both types, call the volume one <gid_model_name>.post.lst and
        #~ the cut one <gid_model_name>_cuts.post.lst.
        #~ additional_frequencies should contain an array of ints N representing
        #~ the frequencies of additional list files. If it is not empty and GidIO
        #~ is configured to print multiple files, extra list files are written,
        #~ listing one in every N output files.'''
        #~ # Get a name for the GiD list file
        #~ # if the model folder is model.gid, the list file should be called
        #~ # model.post.lst
        #~ path, folder_name = os.path.split(os.getcwd())
        #~ model_name, ext = os.path.splitext(folder_name)
        #~ name_base = model_name
        #~ name_ext = ".post.lst"

        #~ # Remove 1 from extra frequencies (and remove duplicates)
        #~ used_frequencies = [1,]
        #~ extra_frequencies = []
        #~ for f in additional_frequencies:
            #~ if f not in used_frequencies:
                #~ used_frequencies.append(f)
                #~ extra_frequencies.append(f)

        #~ if mpi.rank == 0:
            #~ if self.body_io is not None:
                #~ for rank in range(0, mpi.size):
                    #~ print(rank)
                    #~ list_file = open(name_base + "_" + str(rank) + name_ext,"w")

                    #~ if self.multifile_flag == MultiFileFlag.MultipleFiles:
                        #~ list_file.write("Multiple\n")
                    #~ elif self.multifile_flag == MultiFileFlag.SingleFile:
                        #~ list_file.write("Single\n")

                    #~ self.volume_list_files.append( [1,list_file] )

                    #~ if self.multifile_flag == MultiFileFlag.MultipleFiles:
                        #~ for freq in extra_frequencies:
                            #~ list_file_name = "{0}_list_{1}{2}".format(name_base,freq,name_ext)
                            #~ list_file = open(list_file_name,"w")
                            #~ list_file.write("Multiple\n")

                            #~ self.volume_list_files.append( [freq,list_file] )

            #~ name_base = "{0}_cuts".format(model_name)

            #~ if self.cut_io is not None:
                #~ for rank in range(0, mpi.size):
                    #~ list_file = open(name_base + "_" + str(rank) + name_ext,"w")

                    #~ if self.multifile_flag == MultiFileFlag.MultipleFiles:
                        #~ list_file.write("Multiple\n")
                    #~ elif self.multifile_flag == MultiFileFlag.SingleFile:
                        #~ list_file.write("Single\n")

                    #~ self.cut_list_files.append( [1,list_file] )

                    #~ if self.multifile_flag == MultiFileFlag.MultipleFiles:
                        #~ for freq in extra_frequencies:
                            #~ list_file_name = "{0}_list_{1}{2}".format(name_base,freq,name_ext)
                            #~ list_file = open(list_file_name,"w")
                            #~ list_file.write("Multiple\n")

                            #~ self.cut_list_files.append( [freq,list_file] )
        #~ mpi.world.barrier()
