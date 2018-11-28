from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
import KratosMultiphysics.PfemFluidDynamicsApplication as PfemFluid
CheckForPreviousImport()
import gid_output_process

class GiDOutputProcess(gid_output_process.GiDOutputProcess):

    def _InitializeGiDIO(self,gidpost_flags,param):
        '''Initialize GidIO objects (for volume and cut outputs) and related data.'''
        self.volume_file_name = self.base_file_name
        self.cut_file_name = self.volume_file_name+"_cuts"
        self.post_mode = self.__get_gidpost_flag(param, "GiDPostMode", self.__post_mode)
        self.write_deformed_mesh = self.__get_gidpost_flag(param, "WriteDeformedMeshFlag", self.__write_deformed_mesh)
        self.write_conditions = self.__get_gidpost_flag(param,"WriteConditionsFlag",self.__write_conditions)
        self.multifile_flag = self.__get_gidpost_flag(param,"MultiFileFlag", self.__multi_file_flag)

        if self.body_output or self.node_output:
            self.body_io = PfemFluid.PfemFluidGidIO( self.volume_file_name,
                                    self.post_mode,
                                    self.multifile_flag,
                                    self.write_deformed_mesh,
                                    self.write_conditions)

        if self.skin_output or self.num_planes > 0:
            self.cut_io = PfemFluid.PfemFluidGidIO(self.cut_file_name,
                                self.post_mode,
                                self.multifile_flag,
                                self.write_deformed_mesh,
                                WriteConditionsFlag.WriteConditionsOnly) # Cuts are conditions, so we always print conditions in the cut ModelPart


