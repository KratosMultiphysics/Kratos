from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("TrilinosApplication"):
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

from KratosMultiphysics.gid_output_process import GiDOutputProcess

import os

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
    output_name = settings["Parameters"]["output_name"].GetString()
    postprocess_parameters = settings["Parameters"]["postprocess_parameters"]
    return DistributedGiDOutputProcess(model_part, output_name, postprocess_parameters)

class DistributedGiDOutputProcess(GiDOutputProcess):
    def __init__(self, model_part, file_name, param=None):
        super(DistributedGiDOutputProcess, self).__init__(model_part, file_name, param)
        self.serial_file_name = file_name
        self.base_file_name += "_" + str(KM.ParallelEnvironment.GetDefaultDataCommunicator().Rank()) # overwriting the one from the baseclass

    def _InitializeListFiles(self, additional_frequencies):
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

        if self.post_mode == KM.GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == KM.GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == KM.GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"  # ??? CHECK!
        else:
            return # No support for list_files in this format

        if KM.ParallelEnvironment.GetDefaultDataCommunicator().Rank() == 0:
            if self.body_io is not None:
                with open(name_base + name_ext,"w") as list_file:

                    if self.multifile_flag == KM.MultiFileFlag.MultipleFiles:
                        msg = """WARNING: In MPI, printing results in Multiple files
                             is not supported by the .post.lst list file.
                             Results will have to be open manually from GiD.\n"""
                        KM.Logger.PrintWarning("GiDOutputProcessMPI",msg)

                    elif self.multifile_flag == KM.MultiFileFlag.SingleFile:
                        list_file.write("Merge\n")
                        for rank in range(0, KM.ParallelEnvironment.GetDefaultDataCommunicator().Size()):
                            list_file.write(self.serial_file_name+"_"+str(rank)+ext+'\n')

        KM.ParallelEnvironment.GetDefaultDataCommunicator().Barrier()

    def _CreateCuttingUtility(self):
        if not CheckIfApplicationsAvailable("TrilinosApplication"):
            raise Exception("Defining output cuts requires the TrilinosApplication, which appears to be unavailable.")
        self.epetra_comm = KratosTrilinos.CreateCommunicator()
        return KratosTrilinos.TrilinosCuttingApplication(self.epetra_comm)

    def _SetCurrentTimeParameters(self, additional_list_files):
        ''' doing nothing here in MPI'''
        pass
