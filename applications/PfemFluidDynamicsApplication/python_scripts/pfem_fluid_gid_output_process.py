from KratosMultiphysics import *
import KratosMultiphysics.PfemFluidDynamicsApplication as PfemFluid
from KratosMultiphysics import gid_output_process

class GiDOutputProcess(gid_output_process.GiDOutputProcess):

    def __init__(self, model_part, file_name, param=None):
        super().__init__(model_part, file_name, param)
        self.do_first_output = True

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

    def IsOutputStep(self):

        if self.do_first_output:
            self.do_first_output = False
            return True
        else:
            return self.controller.Evaluate()

    def PrintOutput(self):
        if self.point_output_process is not None:
            self.point_output_process.ExecuteBeforeOutputStep()

        # Print the output
        time = self.__get_pretty_time(self.model_part.ProcessInfo[TIME])
        # self.printed_step_count += 1
        self.model_part.ProcessInfo[PRINTED_STEP] = self.printed_step_count

        tolerance=0.0000000000001
        if (time<tolerance):
            time=0

        if self.output_label_is_time:
            label = time
        else:
            label = self.printed_step_count

        self.printed_step_count += 1 # this is after because result0 is the input mesh with no results

        if self.multifile_flag == MultiFileFlag.MultipleFiles:
            self.__write_mesh(label)
            self.__initialize_results(label)

        self.__write_nodal_results(time)
        self.__write_gp_results(time)
        self.__write_nonhistorical_nodal_results(time)
        self.__write_nodal_flags(time)
        self.__write_elemental_conditional_flags(time)

        if self.flush_after_output and self.multifile_flag != MultiFileFlag.MultipleFiles:
            if self.body_io is not None:
                self.body_io.Flush()
            if self.cut_io is not None:
                self.cut_io.Flush()

        if self.multifile_flag == MultiFileFlag.MultipleFiles:
            self.__finalize_results()
            self.__write_step_to_list(label)

        # Schedule next output
        self.controller.Update()

        if self.point_output_process is not None:
            self.point_output_process.ExecuteAfterOutputStep()

