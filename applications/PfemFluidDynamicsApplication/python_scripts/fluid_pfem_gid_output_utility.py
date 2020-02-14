from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PfemFluidDynamicsApplication import *

class GidOutputUtility(object):
    _post_mode = {"Binary": GiDPostMode.GiD_PostBinary,
                 "Ascii": GiDPostMode.GiD_PostAscii,
                 "AsciiZipped": GiDPostMode.GiD_PostAsciiZipped,
                  }

    _write_deformed_mesh = {True: WriteDeformedMeshFlag.WriteDeformed,
                           False: WriteDeformedMeshFlag.WriteUndeformed,
                            }

    _write_conditions = {True: WriteConditionsFlag.WriteConditions,
                        False: WriteConditionsFlag.WriteElementsOnly,
                         }

    _multi_file_flag = {"Multiples": MultiFileFlag.MultipleFiles,
                       "Single": MultiFileFlag.SingleFile,
                        }

    def _set_flag(self, name, value, flag_dict):
        # print "set_flag",name,value,flag_dict[value]
        try:
            flag = flag_dict[value]
        except KeyError:
            msg = """Trying to set GiD IO flag {0} to unkonwn value {1}\n
                     Acceptable values of {0} are:\n""".format(name, value)
            for key in list(flag_dict.keys()):
                msg.append("  {0}\n".format(str(key)))
            raise KeyError(msg)
        self.__dict__[name] = flag

    def __setattr__(self, name, value):
        # print "setattr",name,value
        if name == "post_mode":
            self._set_flag(name, value, self._post_mode)
        elif name == "write_deformed_mesh":
            self._set_flag(name, value, self._write_deformed_mesh)
        elif name == "write_conditions":
            self._set_flag(name, value, self._write_conditions)
        elif name == "multi_file":
            self._set_flag(name, value, self._multi_file_flag)
        else:
            self.__dict__[name] = value

    #
    def __init__(self, file_name, gid_output_settings):

        # set problem name
        self.file_name = file_name

        # set gid options:
        self.post_mode = gid_output_settings.GiDPostMode   # ascii or binary
        self.multi_file = gid_output_settings.GiDMultiFileFlag  # single or multiple files

        # set gid print flags
        self.write_deformed_mesh = gid_output_settings.GiDWriteMeshFlag
        self.write_conditions = gid_output_settings.GiDWriteConditionsFlag

        # set gid print options:
        self.write_particles = gid_output_settings.GiDWriteParticlesFlag

        # set gid input-output class
        self.io = GidIO(self.file_name, self.post_mode, self.multi_file, self.write_deformed_mesh, self.write_conditions)

    #
    #
    def _write_mesh(self, label, model_part):
        self.io.InitializeMesh(label)
        if(self.write_particles):
            self.io.WriteNodeMesh(model_part.GetMesh());
        self.io.WriteMesh(model_part.GetMesh())
        self.io.FinalizeMesh()

    #
    def _write_gp_results(self, label, model_part, variable):
        self.io.PrintOnGaussPoints(variable, model_part, label)

    #
    def _write_nodal_results(self, label, model_part, variable):
        self.io.WriteNodalResults(variable, model_part.Nodes, label, 0)

    #
    #
    def _initialize_results(self, label, model_part):
        self.io.InitializeResults(label, model_part.GetMesh())

    #
    def _finalize_results(self):
        self.io.FinalizeResults()

    #
    def initialize_results(self, model_part):

        if self.multi_file == MultiFileFlag.SingleFile:
            mesh_name = 0.0
            self._write_mesh(mesh_name, model_part)
            self._initialize_results(mesh_name, model_part)

    #
    def write_results(self, model_part, nodal_variables, gp_variables, current_time, current_step, current_id):

        # set multi-file label
        label = current_id

        print("::[GID Output Utility]:: WRITING RESULTS: [ID: ", label, "] [STEP: ", current_step, "] [TIME: %.8e" % current_time, "]")

        # update cut data if necessary
        if self.multi_file == MultiFileFlag.MultipleFiles:
            self._write_mesh(label, model_part)
            self._initialize_results(label, model_part)

        for var in nodal_variables:
            kratos_variable = globals()[var]
            self._write_nodal_results(current_time, model_part, kratos_variable)

        for var in gp_variables:
            kratos_variable = globals()[var]
            self._write_gp_results(current_time, model_part, kratos_variable)

        # flush gid writing
        self.io.Flush()

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self._finalize_results()

        # print("::GID Outuput Utility:: -Run_GID_for_viewing_the_results_of_the_analysis-")

    #
    def finalize_results(self):
        if self.multi_file == MultiFileFlag.SingleFile:
            self.io.FinalizeResults()
