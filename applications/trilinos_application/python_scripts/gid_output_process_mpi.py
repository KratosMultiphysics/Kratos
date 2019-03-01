from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.TrilinosApplication import *
CheckForPreviousImport()

import gid_output_process

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
    output_name = settings["Parameters"]["output_name"].GetString()
    postprocess_parameters = settings["Parameters"]["postprocess_parameters"]
    return GiDOutputProcessMPI(model_part, output_name, postprocess_parameters)

class GiDOutputProcessMPI(gid_output_process.GiDOutputProcess):

    def __init__(self, model_part, file_name, param=None):
        super(GiDOutputProcessMPI, self).__init__(model_part, file_name, param)
        self.serial_file_name = file_name
        self.base_file_name += "_" + str(mpi.rank) # overwriting the one from the baseclass

    # Note: I need to reimplement this just to be able to call the correct function to write the listfile
    def ExecuteInitialize(self):
        result_file_configuration = self.param["result_file_configuration"]
        result_file_configuration.ValidateAndAssignDefaults(self.defaults["result_file_configuration"])

        # If either of these is True, we will have a volume output file
        self.body_output = result_file_configuration["body_output"].GetBool()
        self.node_output = result_file_configuration["node_output"].GetBool()

        # If skin_output is True or we have output planes, we will have a cut output file
        self.skin_output = result_file_configuration["skin_output"].GetBool()
        plane_output_configuration = result_file_configuration["plane_output"] # should be of array type
        self.num_planes = plane_output_configuration.size()

        # Generate the cuts and store them in self.cut_model_part
        if self.skin_output or self.num_planes > 0:
            self.cut_model_part = ModelPart("CutPart")
            self.epetra_comm = CreateCommunicator()
            self.cut_manager = TrilinosCuttingApplication(self.epetra_comm)
            self._initialize_cut_output(plane_output_configuration)

        # Retrieve gidpost flags and setup GiD output tool
        gidpost_flags = result_file_configuration["gidpost_flags"]
        gidpost_flags.ValidateAndAssignDefaults(self.defaults["result_file_configuration"]["gidpost_flags"])

        self._InitializeGiDIO(gidpost_flags,gidpost_flags)

        # Process nodal and gauss point output
        self.nodal_variables = self._GenerateVariableListFromInput(result_file_configuration["nodal_results"])
        self.gauss_point_variables = self._GenerateVariableListFromInput(result_file_configuration["gauss_point_results"])
        self.nodal_nonhistorical_variables = self._GenerateVariableListFromInput(result_file_configuration["nodal_nonhistorical_results"])
        self.nodal_flags = self._GenerateFlagsListFromInput(result_file_configuration["nodal_flags_results"])
        self.nodal_flags_names =[]
        for i in range(result_file_configuration["nodal_flags_results"].size()):
            self.nodal_flags_names.append(result_file_configuration["nodal_flags_results"][i].GetString())
        self.elemental_conditional_flags = self._GenerateFlagsListFromInput(result_file_configuration["elemental_conditional_flags_results"])
        self.elemental_conditional_flags_names =[]
        for i in range(result_file_configuration["elemental_conditional_flags_results"].size()):
            self.elemental_conditional_flags_names.append(result_file_configuration["elemental_conditional_flags_results"][i].GetString())

        # Set up output frequency and format
        output_file_label = result_file_configuration["file_label"].GetString()
        if output_file_label == "time":
            self.output_label_is_time = True
        elif output_file_label == "step":
            self.output_label_is_time = False
        else:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,output_file_label,"file_label")
            raise Exception(msg)

        output_control_type = result_file_configuration["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_is_time = True
        elif output_control_type == "step":
            self.output_control_is_time = False
        else:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,output_file_label,"file_label")
            raise Exception(msg)

        self.output_frequency = result_file_configuration["output_frequency"].GetDouble()

        self.__initialize_list_files()

        # Process point recording data
        if self.point_output_process is not None:
            self.point_output_process.ExecuteInitialize()

    def __initialize_list_files(self):
        '''Set up .post.lst files for global and cut results.
        If we have only one tipe of output (volume or cut), the
        list file is called <gid_model_name>.post.lst. When we have
        both types, call the volume one <gid_model_name>.post.lst and
        the cut one <gid_model_name>_cuts.post.lst.
        Additional frequencies not supported in MPI (it is recommended
        you use single file mode anyways, for more convenient
        visualization.'''
        # Get a name for the GiD list file
        # if the model folder is model.gid, the list file should be called
        # model.post.lst
        path, folder_name = os.path.split(os.getcwd())
        model_name, ext = os.path.splitext(folder_name)
        name_base = model_name
        name_ext = ".post.lst"

        if self.post_mode == GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"  # ??? CHECK!
        else:
            return # No support for list_files in this format

        if mpi.rank == 0:
            if self.body_io is not None:
                with open(name_base + name_ext,"w") as list_file:

                    if self.multifile_flag == MultiFileFlag.MultipleFiles:
                        msg = """WARNING: In MPI, printing results in Multiple files
                             is not supported by the .post.lst list file.
                             Results will have to be open manually from GiD.\n"""
                        print (msg)

                    elif self.multifile_flag == MultiFileFlag.SingleFile:
                        list_file.write("Merge\n")
                        for rank in range(0, mpi.size):
                            list_file.write(self.serial_file_name+"_"+str(rank)+ext+'\n')

        mpi.world.barrier()
