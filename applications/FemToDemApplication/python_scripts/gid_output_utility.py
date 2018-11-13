from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FemToDemApplication import *
CheckForPreviousImport()


class GidOutputUtility(object):
    
    _post_mode = {"GiD_PostBinary": GiDPostMode.GiD_PostBinary,
                  "GiD_PostAscii": GiDPostMode.GiD_PostAscii,
                 }

    _write_deformed_mesh = {"WriteDeformed": WriteDeformedMeshFlag.WriteDeformed,
                            "WriteUndeformed": WriteDeformedMeshFlag.WriteUndeformed,
                           }

    _write_conditions = {True: WriteConditionsFlag.WriteConditions,
                         False: WriteConditionsFlag.WriteElementsOnly,
                        }

    _multi_file_flag = {"MultipleFiles": MultiFileFlag.MultipleFiles,
                        "SingleFile": MultiFileFlag.SingleFile,
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
        elif name == "multi_file_flag":
            self._set_flag(name, value, self._multi_file_flag)
        else:
            self.__dict__[name] = value
            

    def __init__(self, ProjectParameters, file_name, starting_time, ending_time, delta_time):
        ## GiD output initialization
        # set gid options:
        self.post_mode = ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].GetString() # (Binary/Ascii)
        self.multi_file_flag = ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["MultiFileFlag"].GetString() # (Single/Multiples)

        # set gid print flags
        self.write_deformed_mesh = ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["WriteDeformedMeshFlag"].GetString()
        self.write_conditions = False #HARD CODED TODO

        # set gid input-output class
        self.io = GidIO(file_name, self.post_mode, self.multi_file_flag, self.write_deformed_mesh, self.write_conditions)
        
        ## Time operations initialization
        # set ending time
        self.ending_time = ending_time

        # set time frequency
        self.write_frequency = ProjectParameters["AMR_data"]["refinement_frequency"].GetDouble()
        if(self.write_frequency < delta_time):
            self.write_frequency = delta_time

        # set time counters
        self.time_counter = starting_time + self.write_frequency

        # set time operation tolerance
        self.tolerance = delta_time * 1e-10;


    def _write_mesh(self, label, model_part):
        self.io.InitializeMesh(label)
        self.io.WriteMesh(model_part.GetMesh())
        self.io.FinalizeMesh()


    def initialize_results(self, model_part, current_id):
        if self.multi_file_flag == MultiFileFlag.SingleFile:
            self._write_mesh(current_id, model_part)
            self.io.InitializeResults(current_id, model_part.GetMesh())


    def CheckWriteResults(self, current_time):
        write = False

        if(current_time + self.tolerance >= self.time_counter):
            self.time_counter = self.time_counter + self.write_frequency
            write = True
        elif(current_time + self.tolerance >= self.ending_time):
            write = True

        return write
        

    def write_results(self, model_part, nodal_variables, gp_variables, current_time, current_step, current_id):
        print("WRITING RESULTS: ID", current_id, " - STEP", current_step, " - TIME", "%.5f" % current_time)

        # update cut data if necessary
        if self.multi_file_flag == MultiFileFlag.MultipleFiles:
            self._write_mesh(current_id, model_part)
            self.io.InitializeResults(current_id, model_part.GetMesh())

        for var in nodal_variables:
          if var!="NO_RESULT":
            kratos_variable = globals()[var]
            self.io.WriteNodalResults(kratos_variable, model_part.Nodes, current_time, 0)

        for var in gp_variables:
          if var!="NO_RESULT":
            kratos_variable = globals()[var]
            self.io.PrintOnGaussPoints(kratos_variable, model_part, current_time)

        # flush gid writing
        self.io.Flush()

        if self.multi_file_flag == MultiFileFlag.MultipleFiles:
            self.io.FinalizeResults()


    def finalize_results(self):
        if self.multi_file_flag == MultiFileFlag.SingleFile:
            self.io.FinalizeResults()
