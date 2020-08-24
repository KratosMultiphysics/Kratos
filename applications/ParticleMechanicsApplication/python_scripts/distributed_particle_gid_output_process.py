from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if CheckIfApplicationsAvailable("TrilinosApplication"):
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

from KratosMultiphysics.ParticleMechanicsApplication.particle_gid_output_process import ParticleGiDOutputProcess

import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
import os

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
    output_name = settings["Parameters"]["output_name"].GetString()
    postprocess_parameters = settings["Parameters"]["postprocess_parameters"]
    return DistributedParticleGiDOutputProcess(model_part, output_name, postprocess_parameters)

class DistributedParticleGiDOutputProcess(ParticleGiDOutputProcess):
    def __init__(self, model_part, file_name, param=None):
        super(DistributedParticleGiDOutputProcess, self).__init__(model_part, file_name, param)
        self.serial_file_name = file_name
        self.base_file_name += "_" + str(KM.ParallelEnvironment.GetDefaultDataCommunicator().Rank()) # overwriting the one from the baseclass

    def ExecuteBeforeSolutionLoop(self):
        # Initiate Output Mesh
        self.mesh_file = open(self.base_file_name + ".post.msh",'w')
        self.mesh_file.write("MESH \"")
        self.mesh_file.write("outmesh")
        self.mesh_file.write("\" dimension 3 ElemType Point Nnode 1\n")
        self.mesh_file.write("Coordinates\n")
        rank = str(KM.ParallelEnvironment.GetDefaultDataCommunicator().Rank())
        for mpm in self.model_part.Elements:
            coord = mpm.CalculateOnIntegrationPoints(KratosParticle.MP_COORD,self.model_part.ProcessInfo)[0]
            # TODO: find better way than this:
            Id = str(mpm.Id) + str(rank) + str(rank) + str(rank)
            self.mesh_file.write("{} {} {} {}\n".format( Id, coord[0], coord[1], coord[2]))
        self.mesh_file.write("End Coordinates\n")
        self.mesh_file.write("Elements\n")
        for mpm in self.model_part.Elements:
            Id = str(mpm.Id) + str(rank) + str(rank) + str(rank)
            self.mesh_file.write("{} {}\n".format(Id, Id))
        self.mesh_file.write("End Elements\n")
        self.mesh_file.flush()

        # Initiate Output File
        self.result_file = open(self.base_file_name + ".post.res",'w')
        self.result_file.write("GiD Post Results File 1.0\n")

    def _write_mp_results(self, step_label=None):
        clock_time = self._start_time_measure()
        rank = str(KM.ParallelEnvironment.GetDefaultDataCommunicator().Rank())
        for i in range(self.variable_name_list.size()):
            var_name = self.variable_name_list[i].GetString()
            variable = self.variable_list[i]

            is_scalar = self._is_scalar(variable)

            # Write in result file
            self.result_file.write("Result \"")
            self.result_file.write(var_name)

            if is_scalar:
                self.result_file.write('" "Kratos" {} Scalar OnNodes\n'.format(step_label))
            else:
                self.result_file.write('" "Kratos" {} Vector OnNodes\n'.format(step_label))

            self.result_file.write("Values\n")
            for mpm in self.model_part.Elements:
                print_variable = mpm.CalculateOnIntegrationPoints(variable,self.model_part.ProcessInfo)[0]
                # Check whether variable is a scalar or vector
                if isinstance(print_variable, float) or isinstance(print_variable, int):
                    print_size = 1
                else:
                    print_size = print_variable.Size()
                # TODO: find better way than this:
                Id = str(mpm.Id) + str(rank) + str(rank) + str(rank)
                # Write variable as formated
                if print_size == 1:
                    self.result_file.write("{} {}\n".format(Id, print_variable))
                elif print_size == 3:
                    self.result_file.write("{} {} {} {}\n".format(Id, print_variable[0], print_variable[1], print_variable[2]))
                elif print_size == 6:
                    self.result_file.write("{} {} {} {} {} {} {}\n".format(Id, print_variable[0], print_variable[1], print_variable[2], print_variable[3], print_variable[4], print_variable[5]))
                else:
                    KratosMultiphysics.Logger.PrintInfo("Warning in mpm gid output", "Printing format is not defined for variable: ", var_name, "with size: ", print_size)

            self.result_file.write("End Values\n")

        self._stop_time_measure(clock_time)

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
        name_ext = ".post.rest"

        # if self.post_mode == KM.GiDPostMode.GiD_PostBinary:
        #     ext = ".post.bin"
        # elif self.post_mode == KM.GiDPostMode.GiD_PostAscii:
        #     ext = ".post.res"
        # elif self.post_mode == KM.GiDPostMode.GiD_PostAsciiZipped:
        #     ext = ".post.res"  # ??? CHECK!
        # else:
        #     return # No support for list_files in this format

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
