import os
from KratosMultiphysics import *

class GiDOutput(object):
    _post_mode = {"Binary": GiDPostMode.GiD_PostBinary,
                  "Ascii": GiDPostMode.GiD_PostAscii,
                  "AsciiZipped": GiDPostMode.GiD_PostAsciiZipped,
                  }

    _write_deformed_mesh = {True: WriteDeformedMeshFlag.WriteDeformed,
                            False: WriteDeformedMeshFlag.WriteUndeformed,
                            }

    _write_conditions = {True: WriteConditionsFlag.WriteConditions,
                         False: WriteConditionsFlag.WriteElementsOnly,
                         #"CondOnly" : WriteConditionsFlag.WriteConditionsOnly,
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

    def write_step_to_list(self, step_label=None):
        if self.post_mode == GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"  # ??? CHECK!

        with open(self.listfilename, "a") as listfile:
            if step_label is not None:
                listfile.write(
                    "{0}_{1}{2}\n".format(
                        self.filename,
                        step_label,
                        ext))
            else:
                listfile.write("{0}{1}\n".format(self.filename, ext))

    def __init__(self,
                 file_name,
                 vol_output=True,
                 post_mode="Binary",
                 multifile="Single",
                 deformed_mesh=False,
                 write_conditions=True
                 ):

        depr_msg  = '"kratos/python_scripts/gid_output.py" is deprecated and will be removed\n'
        depr_msg += 'Please use "kratos/python_scripts/gid_output_process.py" instead!'
        Logger.PrintWarning('DEPRECATION-WARNING', depr_msg)

        self.filename = file_name
        self.volume_output = vol_output
        self.post_mode = post_mode
        self.write_deformed_mesh = deformed_mesh
        self.write_conditions = write_conditions
        if vol_output is False and write_conditions is False:
            self.write_conditions = True
            msg = "**\n** WARNING -- GiD output:Changing write_conditions from False to True because\n**            volume output is set to False (Cuts are printed as conditions).\n**"
            print(msg)

        self.multi_file = multifile

        # Initialize GiD io
        self.io = GidIO(
            self.filename,
            self.post_mode,
            self.multi_file,
            self.write_deformed_mesh,
            self.write_conditions)

        # Initialize cut data
        self.cut_model_part = None
        self.cut_number = 0

        # Get a name for the GiD list file
        # if the model folder is model.gid, the list file should be called
        # model.post.lst
        path, folder_name = os.path.split(os.getcwd())
        model_name, ext = os.path.splitext(folder_name)
        self.listfilename = "{0}.post.lst".format(model_name)

    def define_cuts(self, model_part, cut_list):
        if self.volume_output:
            raise Exception(
                "GiD IO Error: define_cuts was called for a problem set up to print volume output.")
        elif self.cut_number > 0:
            raise Exception(
                "GiD IO Error: define_cuts was called twice. Please provide all cuts in a single list.")

        if self.cut_model_part is None:
            self.cut_model_part = ModelPart("CutPart")
            self.cut_app = CuttingUtility()
            self.cut_app.FindSmallestEdge(model_part)
            self.cut_app.AddVariablesToCutModelPart(model_part,self.cut_model_part)

        for cut_data in cut_list:
            self.cut_number += 1
            n = Vector(cut_data[0])
            p = Vector(cut_data[1])
            self.cut_app.GenerateCut(
                model_part,
                self.cut_model_part,
                n,
                p,
                self.cut_number,
                0.01)

        self.cut_number += 1
        self.cut_app.AddSkinConditions(
            model_part,
            self.cut_model_part,
            self.cut_number)

    def get_out_model_part(self, model_part):
        if self.volume_output:
            return model_part
        else:
            if self.cut_number == 0:
                raise Exception(
                    "GiD IO Error: writing mesh in cut mode, but no cuts were defined. Please call define_cuts first")
            return self.cut_model_part

    def _write_mesh(self, label, model_part):
        self.io.InitializeMesh(label)
        self.io.WriteMesh(model_part.GetMesh())
        self.io.FinalizeMesh()

    def _initialize_results(self, label, model_part):
        self.io.InitializeResults(label, model_part.GetMesh())

    def _write_nodal_results(self, label, model_part, variable):
        self.io.WriteNodalResults(variable, model_part.Nodes, label, 0)

    def _write_gp_results(self, label, model_part, variable):
        self.io.PrintOnGaussPoints(variable, model_part, label)

    def _finalize_results(self):
        self.io.FinalizeResults()

    def initialize_results(self, model_part):

        if self.multi_file == MultiFileFlag.SingleFile:
            out_model_part = self.get_out_model_part(model_part)
            mesh_name = 0.0

            self._write_mesh(mesh_name, out_model_part)
            self._initialize_results(mesh_name, out_model_part)

        # Initialize list file
        with open(self.listfilename, "w") as listfile:
            if self.multi_file == MultiFileFlag.MultipleFiles:
                listfile.write("Multiple\n")
            elif self.multi_file == MultiFileFlag.SingleFile:
                listfile.write("Single\n")

        if self.multi_file == MultiFileFlag.SingleFile:
            if self.post_mode == GiDPostMode.GiD_PostBinary:
                self.write_step_to_list()
            else:
                self.write_step_to_list(0)

    def write_results(self, label, model_part, nodal_variables, gp_variables):
        # label = str(label) #it should be a C double
        out_model_part = self.get_out_model_part(model_part)

        # update cut data if necessary
        if not self.volume_output:
            self.cut_app.UpdateCutData(out_model_part, model_part)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self._write_mesh(label, out_model_part)
            self._initialize_results(label, out_model_part)

        for var in nodal_variables:
            kratos_variable = KratosGlobals.GetVariable(var)
            self._write_nodal_results(label, out_model_part, kratos_variable)

        for var in gp_variables:
            kratos_variable = KratosGlobals.GetVariable(var)
            self._write_gp_results(label, out_model_part, kratos_variable)

        if self.multi_file == MultiFileFlag.MultipleFiles:
            self._finalize_results()

            with open(self.listfilename, "w") as listfile:
                self.write_step_to_list(label)

    def finalize_results(self):
        if self.multi_file == MultiFileFlag.SingleFile:
            self.io.FinalizeResults()
