from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from KratosMultiphysics import *

from gid_output_process import GiDOutputProcess


## LOCAL DERIVATE CLASS IN ORDER TO BE CONSISTENT WITH RESTART (I.E. THE FREQUENCY WORKS THE SAME IN BOTH )
class PfemGiDOutputProcess(GiDOutputProcess):

    def __init__(self,model_part,file_name,param = None):

        GiDOutputProcess.__init__(self,model_part,file_name,param)



    #Local modification in order to not write the first step
    def ExecuteInitialize(self):


        GiDOutputProcess.ExecuteInitialize(self)

        # Set current time parameters
        if(  self.model_part.ProcessInfo[IS_RESTARTED] == True ):
            self.step_count = self.model_part.ProcessInfo[STEP]
            self.printed_step_count = self.model_part.ProcessInfo[PRINTED_STEP]
            if self.output_control_is_time:
                self.next_output = self.model_part.ProcessInfo[TIME]
            else:
                self.next_output = self.model_part.ProcessInfo[STEP]
        else:
            self.next_output = 0;


        self.next_output += self.output_frequency;

    #Local modification to plot GP data at the first step
    def ExecuteBeforeSolutionLoop(self):
        '''Initialize output meshes.'''

        if self.multifile_flag == MultiFileFlag.SingleFile:
            mesh_name = 0.0
            self._GiDOutputProcess__write_mesh(mesh_name)
            self._GiDOutputProcess__initialize_results(mesh_name)

            if self.post_mode == GiDPostMode.GiD_PostBinary:
                self._GiDOutputProcess__write_step_to_list()
            else:
                self._GiDOutputProcess__write_step_to_list(0)

        if self.multifile_flag == MultiFileFlag.MultipleFiles:
            label = 0.0
            self._GiDOutputProcess__write_mesh(label)
            self._GiDOutputProcess__initialize_results(label)

            self._GiDOutputProcess__write_nodal_results(label)
            self._GiDOutputProcess__write_gp_results(label)
            self._GiDOutputProcess__write_nonhistorical_nodal_results(label)
            self._GiDOutputProcess__write_nodal_flags(label)
            self._GiDOutputProcess__finalize_results()

        if self.point_output_process is not None:
            self.point_output_process.ExecuteBeforeSolutionLoop()


