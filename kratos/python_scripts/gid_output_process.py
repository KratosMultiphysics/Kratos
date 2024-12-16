import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities
from KratosMultiphysics.deprecation_management import DeprecationManager
import os
from pathlib import Path

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
    output_name = settings["Parameters"]["output_name"].GetString()
    postprocess_parameters = settings["Parameters"]["postprocess_parameters"]

    if model_part.IsDistributed():
        from KratosMultiphysics.mpi.distributed_gid_output_process import DistributedGiDOutputProcess
        return DistributedGiDOutputProcess(model_part, output_name, postprocess_parameters)
    else:
        return GiDOutputProcess(model_part, output_name, postprocess_parameters)

class GiDOutputProcess(KM.OutputProcess):

    defaults = KM.Parameters('''{
        "result_file_configuration": {
            "gidpost_flags": {
                "GiDPostMode": "GiD_PostBinary",
                "WriteDeformedMeshFlag": "WriteUndeformed",
                "WriteConditionsFlag": "WriteElementsOnly",
                "MultiFileFlag": "SingleFile"
            },
            "file_label": "time",
            "time_label_format": "{:.12f}",
            "output_control_type": "step",
            "output_interval": 1.0,
            "flush_after_output": false,
            "body_output": true,
            "node_output": false,
            "skin_output": false,
            "plane_output": [],
            "nodal_results": [],
            "nodal_nonhistorical_results": [],
            "nodal_flags_results": [],
            "elemental_conditional_flags_results": [],
            "gauss_point_results": [],
            "additional_list_files": []
        },
        "point_data_configuration": []
    }''')

    default_plane_output_data = KM.Parameters('''{
        "normal": [0.0, 0.0, 0.0],
        "point" : [0.0, 0.0, 0.0]
    }''')

    __post_mode = {
                    # JSON input
                    "GiD_PostAscii":        KM.GiDPostMode.GiD_PostAscii,
                    "GiD_PostAsciiZipped":  KM.GiDPostMode.GiD_PostAsciiZipped,
                    "GiD_PostBinary":       KM.GiDPostMode.GiD_PostBinary,
                    "GiD_PostHDF5":         KM.GiDPostMode.GiD_PostHDF5,
                    # Legacy
                    "Binary":               KM.GiDPostMode.GiD_PostBinary,
                    "Ascii":                KM.GiDPostMode.GiD_PostAscii,
                    "AsciiZipped":          KM.GiDPostMode.GiD_PostAsciiZipped,
                    }

    __write_deformed_mesh = {
                    # JSON input
                    "WriteDeformed":        KM.WriteDeformedMeshFlag.WriteDeformed,
                    "WriteUndeformed":      KM.WriteDeformedMeshFlag.WriteUndeformed,
                    # Legacy
                    True:                   KM.WriteDeformedMeshFlag.WriteDeformed,
                    False:                  KM.WriteDeformedMeshFlag.WriteUndeformed,
                    }

    __write_conditions = {
                    # JSON input
                    "WriteConditions":      KM.WriteConditionsFlag.WriteConditions,
                    "WriteElementsOnly":    KM.WriteConditionsFlag.WriteElementsOnly,
                    "WriteConditionsOnly":  KM.WriteConditionsFlag.WriteConditionsOnly,
                    # Legacy
                    True:                   KM.WriteConditionsFlag.WriteConditions,
                    False:                  KM.WriteConditionsFlag.WriteElementsOnly,
                    }

    __multi_file_flag = {
                    # JSON input
                    "SingleFile":           KM.MultiFileFlag.SingleFile,
                    "MultipleFiles":        KM.MultiFileFlag.MultipleFiles,
                    # Legacy
                    "Multiples":            KM.MultiFileFlag.MultipleFiles,
                    "Single":               KM.MultiFileFlag.SingleFile,
                    }

    def __init__(self,model_part,file_name,param = None):
        super().__init__()

        if param is None:
            param = self.defaults
        else:
            # Warning: we may be changing the parameters object here:
            self.TranslateLegacyVariablesAccordingToCurrentStandard(param)
            # Note: this only validates the first level of the JSON tree.
            # I'm not going for recursive validation because some branches may
            # not exist and I don't want the validator assigning defaults there.
            param.ValidateAndAssignDefaults(self.defaults)

        self.param = param
        self.base_file_name = self._ValidateFileName(file_name)

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
            from KratosMultiphysics.point_output_process import PointOutputProcess
            self.point_output_process = PointOutputProcess(self.model_part,point_data_configuration)
        else:
            self.point_output_process = None

        self.printed_step_count = 0

        controller_settings = KM.Parameters("""{}""")
        controller_settings.AddString("model_part_name", model_part.FullName())
        if param["result_file_configuration"].Has("output_control_type"): controller_settings.AddValue("output_control_type", param["result_file_configuration"]["output_control_type"])
        if param["result_file_configuration"].Has("output_interval"): controller_settings.AddValue("output_interval", param["result_file_configuration"]["output_interval"])
        self.controller = KM.OutputController(model_part.GetModel(), controller_settings)

    # This function can be extended with new deprecated variables as they are generated
    def TranslateLegacyVariablesAccordingToCurrentStandard(self, settings):
        # Defining a string to help the user understand where the warnings come from (in case any is thrown)
        context_string = type(self).__name__

        if settings.Has('result_file_configuration'):
            sub_settings_where_var_is = settings['result_file_configuration']
            old_name = 'output_frequency'
            new_name = 'output_interval'

            if DeprecationManager.HasDeprecatedVariable(context_string, sub_settings_where_var_is, old_name, new_name):
                DeprecationManager.ReplaceDeprecatedVariableName(sub_settings_where_var_is, old_name, new_name)

        if settings.Has('result_file_configuration'):
            sub_settings_where_var_is = settings['result_file_configuration']
            old_name = 'write_properties_id'
            new_name = 'write_ids'

            if DeprecationManager.HasDeprecatedVariable(context_string, sub_settings_where_var_is, old_name, new_name):
                DeprecationManager.ReplaceDeprecatedVariableName(sub_settings_where_var_is, old_name, new_name)

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
            model = self.model_part.GetModel()
            self.cut_model_part = model.CreateModelPart("CutPart")
            self.cut_manager = self._CreateCuttingUtility()
            self._initialize_cut_output(plane_output_configuration)

        # Retrieve gidpost flags and setup GiD output tool
        gidpost_flags = result_file_configuration["gidpost_flags"]
        gidpost_flags.ValidateAndAssignDefaults(self.defaults["result_file_configuration"]["gidpost_flags"])

        self._InitializeGiDIO(gidpost_flags,gidpost_flags)

        # Process nodal and gauss point output
        self.nodal_variables = kratos_utilities.GenerateVariableListFromInput(result_file_configuration["nodal_results"])
        self.gauss_point_variables = kratos_utilities.GenerateVariableListFromInput(result_file_configuration["gauss_point_results"])
        self.nodal_nonhistorical_variables = kratos_utilities.GenerateVariableListFromInput(result_file_configuration["nodal_nonhistorical_results"])
        self.nodal_flags = kratos_utilities.GenerateFlagsListFromInput(result_file_configuration["nodal_flags_results"])
        self.nodal_flags_names =[]
        for i in range(result_file_configuration["nodal_flags_results"].size()):
            self.nodal_flags_names.append(result_file_configuration["nodal_flags_results"][i].GetString())
        self.elemental_conditional_flags = kratos_utilities.GenerateFlagsListFromInput(result_file_configuration["elemental_conditional_flags_results"])
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

        self.time_label_format = result_file_configuration["time_label_format"].GetString()

        output_control_type = result_file_configuration["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_is_time = True
        elif output_control_type == "step":
            self.output_control_is_time = False
        else:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,output_file_label,"file_label")
            raise Exception(msg)

        self.output_interval = result_file_configuration["output_interval"].GetDouble()

        self.flush_after_output = result_file_configuration["flush_after_output"].GetBool()

        # get .post.lst files
        additional_list_file_data = result_file_configuration["additional_list_files"]
        additional_list_files = [ additional_list_file_data[i].GetInt() for i in range(0,additional_list_file_data.size()) ]

        # Set current time parameters
        if self.model_part.ProcessInfo[KM.IS_RESTARTED]:
            self._SetCurrentTimeParameters(additional_list_files)
        else:
            # Create .post.lst files
            self._InitializeListFiles(additional_list_files)

        # Process point recording data
        if self.point_output_process is not None:
            self.point_output_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        '''Initialize output meshes.'''

        if self.multifile_flag == KM.MultiFileFlag.SingleFile:
            mesh_name = 0.0
            self.__write_mesh(mesh_name)
            self.__initialize_results(mesh_name)
            self.__write_step_to_list()

        if self.multifile_flag == KM.MultiFileFlag.MultipleFiles:
            label = 0.0
            self.__write_mesh(label)
            self.__initialize_results(label)
            self.__write_nodal_results(label)
            self.__write_gp_results(label)
            self.__write_nonhistorical_nodal_results(label)
            self.__write_nodal_flags(label)
            self.__write_elemental_conditional_flags(label)
            self.__finalize_results()

            if self.output_label_is_time:
                file_label = 0.0
            else:
                file_label = 0
            self.__write_step_to_list(file_label)

        if self.point_output_process is not None:
            self.point_output_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        if self.point_output_process is not None:
            self.point_output_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        if self.point_output_process is not None:
            self.point_output_process.ExecuteFinalizeSolutionStep()

    def IsOutputStep(self):
        return self.controller.Evaluate()

    def PrintOutput(self):
        if self.point_output_process is not None:
            self.point_output_process.ExecuteBeforeOutputStep()
            if self.point_output_process.IsOutputStep():
                self.point_output_process.PrintOutput()

        # Print the output
        time = self.__get_pretty_time(self.model_part.ProcessInfo[KM.TIME])
        self.printed_step_count += 1
        self.model_part.ProcessInfo[KM.PRINTED_STEP] = self.printed_step_count
        if self.output_label_is_time:
            label = time
        else:
            label = self.printed_step_count

        if self.multifile_flag == KM.MultiFileFlag.MultipleFiles:
            self.__write_mesh(label)
            self.__initialize_results(label)

        self.__write_nodal_results(time)
        self.__write_gp_results(time)
        self.__write_nonhistorical_nodal_results(time)
        self.__write_nodal_flags(time)
        self.__write_elemental_conditional_flags(time)

        if self.flush_after_output and self.multifile_flag != KM.MultiFileFlag.MultipleFiles:
            if self.body_io is not None:
                self.body_io.Flush()
            if self.cut_io is not None:
                self.cut_io.Flush()

        if self.multifile_flag == KM.MultiFileFlag.MultipleFiles:
            self.__finalize_results()
            self.__write_step_to_list(label)

        # Schedule next output
        self.controller.Update()

        if self.point_output_process is not None:
            self.point_output_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        '''Finalize files and free resources.'''

        if self.multifile_flag == KM.MultiFileFlag.SingleFile:
            self.__finalize_results()

        if self.point_output_process is not None:
            self.point_output_process.ExecuteFinalize()

        for freq,f in self.volume_list_files:
            f.close()
        for freq,f in self.cut_list_files:
            f.close()

        # Note: it is important to call the GidIO destructor, since it closes output files
        # Since Python's garbage collection DOES NOT ensure that the destructor will be called,
        # I'm deallocating the GidIO instances explicitly. This is VERY BAD PRACTICE
        # and effectively breaks the class if called after this point, but we haven't found
        # a better solution yet (jcotela 12/V/2016)
        del self.body_io
        del self.cut_io

    def _InitializeGiDIO(self,gidpost_flags,param):
        '''Initialize GidIO objects (for volume and cut outputs) and related data.'''
        self.volume_file_name = self.base_file_name
        self.cut_file_name = self.volume_file_name+"_cuts"
        self.post_mode = self.__get_gidpost_flag(param, "GiDPostMode", self.__post_mode)
        self.write_deformed_mesh = self.__get_gidpost_flag(param, "WriteDeformedMeshFlag", self.__write_deformed_mesh)
        self.write_conditions = self.__get_gidpost_flag(param,"WriteConditionsFlag",self.__write_conditions)
        self.multifile_flag = self.__get_gidpost_flag(param,"MultiFileFlag", self.__multi_file_flag)

        if self.body_output or self.node_output:
            self.body_io = KM.GidIO( self.volume_file_name,
                                    self.post_mode,
                                    self.multifile_flag,
                                    self.write_deformed_mesh,
                                    self.write_conditions,
                                    self.param["result_file_configuration"]["gauss_point_results"].size()>0)

        if self.skin_output or self.num_planes > 0:
            self.cut_io = KM.GidIO(self.cut_file_name,
                                self.post_mode,
                                self.multifile_flag,
                                self.write_deformed_mesh,
                                KM.WriteConditionsFlag.WriteConditionsOnly,
                                self.param["result_file_configuration"]["gauss_point_results"].size()>0) # Cuts are conditions, so we always print conditions in the cut ModelPart

    def __get_pretty_time(self,time):
        pretty_time = self.time_label_format.format(time)
        pretty_time = float(pretty_time)
        return pretty_time

    def __get_gidpost_flag(self, param, label, dictionary):
        '''Parse gidpost settings using an auxiliary dictionary of acceptable values.'''

        keystring = param[label].GetString()
        try:
            value = dictionary[keystring]
        except KeyError:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,keystring,label)
            raise Exception(msg)

        return value

    def _initialize_cut_output(self,plane_output_configuration):
        '''Set up tools used to produce output in skin and cut planes.'''

        self.cut_manager.FindSmallestEdge(self.model_part)
        self.cut_manager.AddVariablesToCutModelPart(self.model_part,self.cut_model_part)
        if self.skin_output:
            self.cut_manager.AddSkinConditions(self.model_part,self.cut_model_part,self.output_surface_index)
            self.output_surface_index += 1

        for plane_id in range(0,plane_output_configuration.size()):

            cut_data = plane_output_configuration[plane_id]
            # This part of the input data is validated term by term
            cut_data.ValidateAndAssignDefaults(self.default_plane_output_data)

            self.__define_output_plane(cut_data)

    def _InitializeListFiles(self,additional_frequencies):
        '''Set up .post.lst files for global and cut results.
        If we have only one type of output (volume or cut), the
        list file is called <gid_model_name>.post.lst. When we have
        both types, call the volume one <gid_model_name>.post.lst and
        the cut one <gid_model_name>_cuts.post.lst.
        additional_frequencies should contain an array of ints N representing
        the frequencies of additional list files. If it is not empty and GidIO
        is configured to print multiple files, extra list files are written,
        listing one in every N output files.'''
        # Get a name for the GiD list file
        # if the model folder is model.gid, the list file should be called
        # model.post.lst
        path, folder_name = os.path.split(os.getcwd())
        model_name, ext = os.path.splitext(folder_name)
        name_base = model_name
        name_ext = ".post.lst"

        # Remove 1 from extra frequencies (and remove duplicates)
        used_frequencies = [1,]
        extra_frequencies = []
        for f in additional_frequencies:
            if f not in used_frequencies:
                used_frequencies.append(f)
                extra_frequencies.append(f)

        if self.body_io is not None:
            list_file = open(name_base+name_ext,"w")

            if self.multifile_flag == KM.MultiFileFlag.MultipleFiles:
                list_file.write("Multiple\n")
            elif self.multifile_flag == KM.MultiFileFlag.SingleFile:
                list_file.write("Single\n")

            self.volume_list_files.append( [1,list_file] )

            if self.multifile_flag == KM.MultiFileFlag.MultipleFiles:
                for freq in extra_frequencies:
                    list_file_name = "{0}_list_{1}{2}".format(name_base,freq,name_ext)
                    list_file = open(list_file_name,"w")
                    list_file.write("Multiple\n")

                    self.volume_list_files.append( [freq,list_file] )

            name_base = "{0}_cuts".format(model_name)

        if self.cut_io is not None:

            list_file = open(name_base+name_ext,"w")

            if self.multifile_flag == KM.MultiFileFlag.MultipleFiles:
                list_file.write("Multiple\n")
            elif self.multifile_flag == KM.MultiFileFlag.SingleFile:
                list_file.write("Single\n")

            self.cut_list_files.append( [1,list_file] )

            if self.multifile_flag == KM.MultiFileFlag.MultipleFiles:
                for freq in extra_frequencies:
                    list_file_name = "{0}_list_{1}{2}".format(name_base,freq,name_ext)
                    list_file = open(list_file_name,"w")
                    list_file.write("Multiple\n")

                    self.cut_list_files.append( [freq,list_file] )

    def __define_output_plane(self,cut_data):
        '''Add a plane to the output plane list.'''

        normal_data = cut_data["normal"]
        n = KM.Vector(3)
        n[0] = normal_data[0].GetDouble()
        n[1] = normal_data[1].GetDouble()
        n[2] = normal_data[2].GetDouble()
        # Sanity check on normal
        norm2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2]
        if norm2 < 1e-12:
            raise Exception("{0} Error: something went wrong in output plane definition: plane {1} has a normal with module 0.".format(self.__class__.__name__,self.output_surface_index))

        point_data = cut_data["point"]
        p = KM.Vector(3)
        p[0] = point_data[0].GetDouble()
        p[1] = point_data[1].GetDouble()
        p[2] = point_data[2].GetDouble()

        # And finally, define the plane
        self.output_surface_index += 1
        self.cut_manager.GenerateCut(   self.model_part,
                                        self.cut_model_part,
                                        n,
                                        p,
                                        self.output_surface_index,
                                        0.01)

    def __write_mesh(self, label):
        if self.body_io is not None:
            self.body_io.InitializeMesh(label)
            if self.body_output:
                self.body_io.WriteMesh(self.model_part.GetMesh())
            if self.node_output:
                self.body_io.WriteNodeMesh(self.model_part.GetMesh())
            self.body_io.FinalizeMesh()

        if self.cut_io is not None:
            self.cut_io.InitializeMesh(label)
            self.cut_io.WriteMesh(self.cut_model_part.GetMesh())
            self.cut_io.FinalizeMesh()

    def __initialize_results(self, label):
        if self.body_io is not None:
            self.body_io.InitializeResults(label, self.model_part.GetMesh())

        if self.cut_io is not None:
            self.cut_io.InitializeResults(label,self.cut_model_part.GetMesh())

    def __write_nodal_results(self, label):
        if self.body_io is not None:
            for variable in self.nodal_variables:
                self.body_io.WriteNodalResults(variable, self.model_part.GetCommunicator().LocalMesh().Nodes, label, 0)

        if self.cut_io is not None:
            self.cut_manager.UpdateCutData(self.cut_model_part, self.model_part)
            for variable in self.nodal_variables:
                self.cut_io.WriteNodalResults(variable, self.cut_model_part.GetCommunicator().LocalMesh().Nodes, label, 0)

    def __write_gp_results(self, label):
        #if self.body_io is not None:
        if self.body_output: # Note: if we only print nodes, there are no GaussPoints
            for variable in self.gauss_point_variables:
                self.body_io.PrintOnGaussPoints(variable, self.model_part, label)

        # Gauss point results depend on the type of element!
        # they are not implemented for cuts (which are generic Condition3D)

    def __write_nonhistorical_nodal_results(self, label):
        if self.body_io is not None:
            for variable in self.nodal_nonhistorical_variables:
                self.body_io.WriteNodalResultsNonHistorical(variable, self.model_part.Nodes, label)

        if self.cut_io is not None:
            self.cut_manager.UpdateCutData(self.cut_model_part, self.model_part)
            for variable in self.nodal_nonhistorical_variables:
                self.cut_io.WriteNodalResultsNonHistorical(variable, self.cut_model_part.Nodes, label)

    def __write_nodal_flags(self, label):
        if self.body_io is not None:
            for flag in range(len(self.nodal_flags)):
                self.body_io.WriteNodalFlags(self.nodal_flags[flag], self.nodal_flags_names[flag], self.model_part.Nodes, label)

        if self.cut_io is not None:
            self.cut_manager.UpdateCutData(self.cut_model_part, self.model_part)
            for flag in range(len(self.nodal_flags)):
                self.cut_io.WriteNodalFlags(self.nodal_flags[flag], self.nodal_flags_names[flag], self.cut_model_part.Nodes, label)

    def __write_elemental_conditional_flags(self, label):
        if self.body_io is not None:
            for flag in range(len(self.elemental_conditional_flags)):
                self.body_io.PrintFlagsOnGaussPoints(self.elemental_conditional_flags[flag], self.elemental_conditional_flags_names[flag], self.model_part, label)

        if self.cut_io is not None:
            self.cut_manager.UpdateCutData(self.cut_model_part, self.model_part)
            for flag in range(len(self.elemental_conditional_flags)):
                self.cut_io.PrintFlagsOnGaussPoints(self.elemental_conditional_flags[flag], self.elemental_conditional_flags_names[flag], self.cut_model_part, label)

    def __finalize_results(self):
        if self.body_io is not None:
            self.body_io.FinalizeResults()

        if self.cut_io is not None:
            self.cut_io.FinalizeResults()

    def __write_step_to_list(self,step_label=None):
        if self.post_mode == KM.GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == KM.GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == KM.GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"  # ??? CHECK!
        else:
            return # No support for list_files in this format

        if step_label is None:
            pretty_label = ""
        elif self.output_label_is_time:
            pretty_label = "_{0:.12g}".format(step_label) # floating point format
        else:
            pretty_label = "_{0}".format(step_label) # int format

        if self.body_io is not None:
            for freq,f in self.volume_list_files:
                if (self.printed_step_count % freq) == 0:
                    f.write("{0}{1}{2}\n".format(self.volume_file_name,pretty_label,ext))
                    f.flush()

        if self.cut_io is not None:
            for freq,f in self.cut_list_files:
                if (self.printed_step_count % freq) == 0:
                    f.write("{0}{1}{2}\n".format(self.cut_file_name,pretty_label,ext))
                    f.flush()

    def __restart_list_files(self,additional_frequencies):
        self.__remove_list_files()

        self._InitializeListFiles(additional_frequencies)

        if self.post_mode == KM.GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == KM.GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == KM.GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"
        else:
            return # No support for list_files in this format


        # Rebuild List Files from existing problem build
        file_id   = []

        path = os.getcwd()

        for f in os.listdir(path):

            if(f.endswith(ext)):

                #if f.name = problem_tested_145.post.bin
                file_parts = f.split('_')  # you get ["problem","tested","145.post.bin"]
                num_parts  = len(file_parts)

                end_parts  = file_parts[num_parts-1].split(".") # you get ["145","post","bin"]
                print_id   = end_parts[0] # you get "145"

                try:
                    label = int(print_id)
                    if label != 0:
                        file_id.append(label)
                except ValueError:
                    # Whatever we got was not convertible to int, probably the input file has a different format
                    pass

            file_id.sort()

        for step_label in file_id:

            if step_label is None:
                pretty_label = ""
            elif self.output_label_is_time:
                pretty_label = "_{0:.12g}".format(step_label) # floating point format
            else:
                pretty_label = "_{0}".format(step_label) # int format

            if self.body_io is not None:
                for freq,f in self.volume_list_files:
                    list_file_name = "{0}{1}{2}\n".format(self.volume_file_name,pretty_label,ext)
                    if (step_label % freq) == 0:
                        f.write(list_file_name)
                        f.flush()

            if self.cut_io is not None:
                for freq,f in self.cut_list_files:
                    list_file_name = "{0}{1}{2}\n".format(self.cut_file_name,pretty_label,ext)
                    if (step_label % freq) == 0:
                        f.write(list_file_name)
                        f.flush()

    def __remove_list_files(self):
        path = os.getcwd()

        # remove previous list files:
        if(os.path.exists(path) == False):
            print(" Problem Path do not exists , check the Problem Path selected ")
        else:
            filelist = [f for f in os.listdir(os.getcwd()) if f.endswith(".lst")]
            for f in filelist:
                if(os.path.exists(f)):
                    try:
                        os.remove(f)
                    except OSError:
                        pass

    def __remove_post_results_files(self, step_label):
        path = os.getcwd()

        if self.post_mode == KM.GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == KM.GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == KM.GiDPostMode.GiD_PostAsciiZipped:
            ext = ".post.res"

        # remove post result files:
        if(os.path.exists(path) == False):
            print(" Problem Path do not exists , check the Problem Path selected ")
        else:
            filelist = []
            for f in os.listdir(os.getcwd()):
                if(f.endswith(ext)):

                    #if f.name = problem_tested_145.post.bin
                    file_parts = f.split('_')
                    # you get the parts ["problem","tested","145.post.bin"]
                    num_parts  = len(file_parts)
                    # take the last part
                    end_parts  = file_parts[num_parts-1].split(".")
                    # you get the parts ["145","post","bin"]
                    print_id   = end_parts[0]
                    # you get "145"

                    try:
                        float(print_id)
                    except ValueError:
                        break

                    if( float(print_id) >  float(step_label) ):
                        filelist.append(f)

            for f in filelist:
                if(os.path.exists(f)):
                    try:
                        os.remove(f)
                    except OSError:
                        pass

    def _CreateCuttingUtility(self):
        return KM.CuttingUtility()

    def _SetCurrentTimeParameters(self, additional_list_files):
        self.printed_step_count = self.model_part.ProcessInfo[KM.PRINTED_STEP]

        # Remove post results
        if self.output_label_is_time:
            label = self.model_part.ProcessInfo[KM.TIME]
        else:
            label = self.printed_step_count

        self.__remove_post_results_files(label)

        # Restart .post.lst files
        self.__restart_list_files(additional_list_files) # FIXME

    @staticmethod
    def _ValidateFileName(file_name):
        # A file name must be specified
        if not file_name:
            raise Exception('No "file_name" was specified!')

        # Make sure that the path to the desired output folder exists
        output_path = Path(file_name).parent
        KM.FilesystemExtensions.MPISafeCreateDirectories(output_path)

        return file_name
