import os
import KratosMultiphysics
from KratosMultiphysics.deprecation_management import DeprecationManager

class GiDDamOutputProcess(KratosMultiphysics.Process):

    defaults = KratosMultiphysics.Parameters('''{
        "result_file_configuration": {
            "gidpost_flags": {
                "GiDPostMode": "GiD_PostBinary",
                "WriteDeformedMeshFlag": "WriteUndeformed",
                "WriteConditionsFlag": "WriteElementsOnly",
                "MultiFileFlag": "SingleFile"
            },
            "output_control_type": "time_s",
            "output_interval": 1.0,
            "start_output_results": 0,
            "body_output": true,
            "node_output": false,
            "skin_output": false,
            "plane_output": [],
            "nodal_results": [],
            "nodal_nonhistorical_results": [],
            "nodal_flags_results": [],
            "gauss_point_results": [],
            "additional_list_files": []
        },
        "point_data_configuration": []
    }''')

    default_plane_output_data = KratosMultiphysics.Parameters('''{
        "normal": [0.0, 0.0, 0.0],
        "point" : [0.0, 0.0, 0.0]
    }''')

    __post_mode = {
                    # JSON input
                    "GiD_PostAscii":        KratosMultiphysics.GiDPostMode.GiD_PostAscii,
                    "GiD_PostAsciiZipped":  KratosMultiphysics.GiDPostMode.GiD_PostAsciiZipped,
                    "GiD_PostBinary":       KratosMultiphysics.GiDPostMode.GiD_PostBinary,
                    "GiD_PostHDF5":         KratosMultiphysics.GiDPostMode.GiD_PostHDF5,
                    # Legacy
                    "Binary":               KratosMultiphysics.GiDPostMode.GiD_PostBinary,
                    "Ascii":                KratosMultiphysics.GiDPostMode.GiD_PostAscii,
                    "AsciiZipped":          KratosMultiphysics.GiDPostMode.GiD_PostAsciiZipped,
                    }

    __write_deformed_mesh = {
                    # JSON input
                    "WriteDeformed":        KratosMultiphysics.WriteDeformedMeshFlag.WriteDeformed,
                    "WriteUndeformed":      KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed,
                    # Legacy
                    True:                   KratosMultiphysics.WriteDeformedMeshFlag.WriteDeformed,
                    False:                  KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed,
                    }

    __write_conditions = {
                    # JSON input
                    "WriteConditions":      KratosMultiphysics.WriteConditionsFlag.WriteConditions,
                    "WriteElementsOnly":    KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly,
                    "WriteConditionsOnly":  KratosMultiphysics.WriteConditionsFlag.WriteConditionsOnly,
                    # Legacy
                    True:                   KratosMultiphysics.WriteConditionsFlag.WriteConditions,
                    False:                  KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly,
                    }

    __multi_file_flag = {
                    # JSON input
                    "SingleFile":           KratosMultiphysics.MultiFileFlag.SingleFile,
                    "MultipleFiles":        KratosMultiphysics.MultiFileFlag.MultipleFiles,
                    # Legacy
                    "Multiples":            KratosMultiphysics.MultiFileFlag.MultipleFiles,
                    "Single":               KratosMultiphysics.MultiFileFlag.SingleFile,
                    }

    def __init__(self,model_part,file_name,start_time,param = None):
        KratosMultiphysics.Process.__init__(self)
        if param is None:
            param = self.defaults
        else:
            # Warning: we may be changing the parameters object here:
            self.TranslateLegacyVariablesAccordingToCurrentStandard(param)
            # Note: this only validates the first level of the JSON tree.
            # I'm not going for recursive validation because some branches may
            # not exist and I don't want the validator assinging defaults there.
            param.ValidateAndAssignDefaults(self.defaults)

        self.param = param
        self.base_file_name = file_name

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
            from KratosMultiphysics import point_output_process
            self.point_output_process = point_output_process.PointOutputProcess(self.model_part,point_data_configuration)
        else:
            self.point_output_process = None

        self.start_time = start_time
        self.step_count = 0
        self.printed_step_count = 0
        self.next_output = 0.0

    @classmethod
    def Flush(cls,a):
        a.flush()

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
            self.__initialize_cut_output(plane_output_configuration)

        # Retrieve gidpost flags and setup GiD output tool
        gidpost_flags = result_file_configuration["gidpost_flags"]
        gidpost_flags.ValidateAndAssignDefaults(self.defaults["result_file_configuration"]["gidpost_flags"])

        self.__initialize_gidio(gidpost_flags,gidpost_flags)

        # Process nodal and gauss point output
        self.nodal_variables = self.__generate_variable_list_from_input(result_file_configuration["nodal_results"])
        self.gauss_point_variables = self.__generate_variable_list_from_input(result_file_configuration["gauss_point_results"])
        self.nodal_nonhistorical_variables = self.__generate_variable_list_from_input(result_file_configuration["nodal_nonhistorical_results"])
        self.nodal_flags = self.__generate_flags_list_from_input(result_file_configuration["nodal_flags_results"])
        self.nodal_flags_names =[]
        for i in range(result_file_configuration["nodal_flags_results"].size()):
            self.nodal_flags_names.append(result_file_configuration["nodal_flags_results"][i].GetString())

        self.output_frequency = result_file_configuration["output_interval"].GetDouble()
        self.start_output_results = result_file_configuration["start_output_results"].GetDouble()

        if self.start_time >= self.start_output_results:
            self.start_output_results = self.start_time
            self.next_output += self.start_output_results + self.output_frequency
        else:
            self.next_output += self.start_output_results

        ## get .post.lst files
        #additional_list_file_data = result_file_configuration["additional_list_files"]
        #additional_list_files = [ additional_list_file_data[i].GetInt() for i in range(0,additional_list_file_data.size()) ]


        # Set current time parameters
        if(  self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.step_count = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
            self.printed_step_count = self.model_part.ProcessInfo[KratosMultiphysics.PRINTED_STEP]

            self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

            label = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

            self.__remove_post_results_files(label)

            ## Restart .post.lst files
            #self.__restart_list_files(additional_list_files)
        #else:
            ## Create .post.lst files
            #self.__initialize_list_files(additional_list_files)

        # Process point recording data
        if self.point_output_process is not None:
            self.point_output_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        '''Initialize output meshes.'''
        label = max(self.start_output_results, self.start_time)
        if self.multifile_flag == KratosMultiphysics.MultiFileFlag.SingleFile:
            mesh_name = 0
            self.__write_mesh(mesh_name)
            self.__initialize_results(mesh_name)
            if (self.start_output_results == 0) or (not self.start_time == 0):
                self.__write_nodal_results(label)
                self.__write_gp_results(label)
                self.__write_nonhistorical_nodal_results(label)
                self.__write_nodal_flags(label)

        if self.multifile_flag == KratosMultiphysics.MultiFileFlag.MultipleFiles:
            if (self.start_output_results == 0) or (not self.start_time == 0):
                self.__write_mesh(label)
                self.__initialize_results(label)
                self.__write_nodal_results(label)
                self.__write_nonhistorical_nodal_results(label)
                self.__write_nodal_flags(label)
                self.__finalize_results()

            result_file_configuration = self.param["result_file_configuration"]
            output_control_type = result_file_configuration["output_control_type"].GetString()
            if output_control_type == "time_s":
                self.multifiles = (
                    MultifileList(self.base_file_name, 1),
                    MultifileList(self.base_file_name, 60),
                    MultifileList(self.base_file_name, 3600),
                    MultifileList(self.base_file_name, 43200),
                    MultifileList(self.base_file_name, 86400),
                )
            elif output_control_type == "time_h":
                self.multifiles = (
                    MultifileList(self.base_file_name, 1),
                    MultifileList(self.base_file_name, 12),
                    MultifileList(self.base_file_name, 24),
                    MultifileList(self.base_file_name, 168),
                )
            elif output_control_type == "time_d":
                self.multifiles = (
                    MultifileList(self.base_file_name, 1),
                    MultifileList(self.base_file_name, 7),
                )
            elif output_control_type == "time_w":
                self.multifiles = (
                    MultifileList(self.base_file_name, 1),
                )
            elif output_control_type == "step":
                self.multifiles = (
                    MultifileList(self.base_file_name, 1),
                )
            else:
                msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,output_control_type,"output_control_type")
                raise Exception(msg)

            self.multifilelists = []
            self.__set_multifile_lists(self.multifiles)

            if self.start_output_results == 0:
                self.__write_inital_step_in_multifile_lists(label)

        if self.point_output_process is not None:
            self.point_output_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1

        if self.point_output_process is not None:
            self.point_output_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        if self.point_output_process is not None:
            self.point_output_process.ExecuteFinalizeSolutionStep()

    def IsOutputStep(self):

        result_file_configuration = self.param["result_file_configuration"]
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if result_file_configuration["output_control_type"].GetString() == "time_s":
            time = time
        elif result_file_configuration["output_control_type"].GetString() == "time_h":
            time = time/3600.0
        elif result_file_configuration["output_control_type"].GetString() == "time_d":
            time = time/86400.0
        elif result_file_configuration["output_control_type"].GetString() == "time_w":
            time = time/604800.0

        return ( time >= self.next_output )

    def PrintOutput(self):

        result_file_configuration = self.param["result_file_configuration"]
        self.output_frequency = result_file_configuration["output_interval"].GetDouble()

        if self.point_output_process is not None:
            self.point_output_process.ExecuteBeforeOutputStep()

        # Print the output
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.printed_step_count += 1
        self.model_part.ProcessInfo[KratosMultiphysics.PRINTED_STEP] = self.printed_step_count

        if result_file_configuration["output_control_type"].GetString() == "time_s":
            label = time
            time = time
        elif result_file_configuration["output_control_type"].GetString() == "time_h":
            label = time/3600
            time = time/3600.0
        elif result_file_configuration["output_control_type"].GetString() == "time_d":
            label = time/86400
            time = time/86400.0
        elif result_file_configuration["output_control_type"].GetString() == "time_w":
            label = time/604800
            time = time/604800.0

        if self.multifile_flag == KratosMultiphysics.MultiFileFlag.MultipleFiles:
            self.__write_mesh(label)
            self.__initialize_results(label)

        self.__write_nodal_results(time)
        self.__write_gp_results(time)
        self.__write_nonhistorical_nodal_results(time)
        self.__write_nodal_flags(time)

        if self.multifile_flag == KratosMultiphysics.MultiFileFlag.MultipleFiles:
            self.__finalize_results()
            self.__write_multifile_lists(label)

        # Schedule next output
        if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
            while self.next_output <= time:
                self.next_output += self.output_frequency

        if self.point_output_process is not None:
            self.point_output_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        '''Finalize files and free resources.'''

        if self.multifile_flag == KratosMultiphysics.MultiFileFlag.SingleFile:
            self.__finalize_results()

        if self.multifile_flag == KratosMultiphysics.MultiFileFlag.MultipleFiles:
            self.__close_multifiles()

        if self.point_output_process is not None:
            self.point_output_process.ExecuteFinalize()

        # NOTE (PR): The followin lines were commented due to a warnign from 'codacy' while doing a PR. In the original 'gid_outpu_process.py these lines are uncommented.'
        #~ for freq,f in self.volume_list_files:
            #~ f.close()
        #~ for freq,f in self.cut_list_files:
            #~ f.close()

        # Note: it is important to call the GidIO destructor, since it closes output files
        # Since Python's garbage colletion DOES NOT ensure that the destructor will be called,
        # I'm deallocating the GidIO instances explicitly. This is VERY BAD PRACTICE
        # and effectively breaks the class if called after this point, but we haven't found
        # a better solution yet (jcotela 12/V/2016)
        del self.body_io
        del self.cut_io

    def __initialize_gidio(self,gidpost_flags,param):
        '''Initialize GidIO objects (for volume and cut outputs) and related data.'''
        self.volume_file_name = self.base_file_name
        self.cut_file_name = self.volume_file_name+"_cuts"
        self.post_mode = self.__get_gidpost_flag(param, "GiDPostMode", self.__post_mode)
        self.write_deformed_mesh = self.__get_gidpost_flag(param, "WriteDeformedMeshFlag", self.__write_deformed_mesh)
        self.write_conditions = self.__get_gidpost_flag(param,"WriteConditionsFlag",self.__write_conditions)
        self.multifile_flag = self.__get_gidpost_flag(param,"MultiFileFlag", self.__multi_file_flag)

        if self.body_output or self.node_output:
            self.body_io = KratosMultiphysics.GidIO(self.volume_file_name,
                                                    self.post_mode,
                                                    self.multifile_flag,
                                                    self.write_deformed_mesh,
                                                    self.write_conditions)

        if self.skin_output or self.num_planes > 0:
            self.cut_io = KratosMultiphysics.GidIO(self.cut_file_name,
                                                   self.post_mode,
                                                   self.multifile_flag,
                                                   self.write_deformed_mesh,
                                                   KratosMultiphysics.WriteConditionsFlag.WriteConditionsOnly) # Cuts are conditions, so we always print conditions in the cut ModelPart

    def __get_gidpost_flag(self, param, label, dictionary):
        '''Parse gidpost settings using an auxiliary dictionary of acceptable values.'''

        keystring = param[label].GetString()
        try:
            value = dictionary[keystring]
        except KeyError:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,value,label)
            raise Exception(msg)

        return value

    def __initialize_cut_output(self,plane_output_configuration):
        '''Set up tools used to produce output in skin and cut planes.'''

        self.cut_model_part = ModelPart("CutPart")
        self.cut_manager = CuttingUtility()
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

    def __initialize_list_files(self,additional_frequencies):
        '''Set up .post.lst files for global and cut results.
        If we have only one tipe of output (volume or cut), the
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
        # NOTE (PR): In the followin lines, two variables were removed (path and ext). See 'gid_output_process.py'. They were put only because the functions have two output arguments. We can call just the one we need.
        folder_name = os.path.split(os.getcwd())[1]
        model_name = os.path.splitext(folder_name)[0]
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

            if self.multifile_flag == KratosMultiphysics.MultiFileFlag.MultipleFiles:
                list_file.write("Multiple\n")
            elif self.multifile_flag == KratosMultiphysics.MultiFileFlag.SingleFile:
                list_file.write("Single\n")

            self.volume_list_files.append( [1,list_file] )

            if self.multifile_flag == KratosMultiphysics.MultiFileFlag.MultipleFiles:
                for freq in extra_frequencies:
                    list_file_name = "{0}_list_{1}{2}".format(name_base,freq,name_ext)
                    list_file = open(list_file_name,"w")
                    list_file.write("Multiple\n")

                    self.volume_list_files.append( [freq,list_file] )

            name_base = "{0}_cuts".format(model_name)

        if self.cut_io is not None:

            list_file = open(name_base+name_ext,"w")

            if self.multifile_flag == KratosMultiphysics.MultiFileFlag.MultipleFiles:
                list_file.write("Multiple\n")
            elif self.multifile_flag == KratosMultiphysics.MultiFileFlag.SingleFile:
                list_file.write("Single\n")

            self.cut_list_files.append( [1,list_file] )

            if self.multifile_flag == KratosMultiphysics.MultiFileFlag.MultipleFiles:
                for freq in extra_frequencies:
                    list_file_name = "{0}_list_{1}{2}".format(name_base,freq,name_ext)
                    list_file = open(list_file_name,"w")
                    list_file.write("Multiple\n")

                    self.cut_list_files.append( [freq,list_file] )

    def __define_output_plane(self,cut_data):
        '''Add a plane to the output plane list.'''

        normal_data = cut_data["normal"]
        n = Vector(3)
        n[0] = normal_data[0].GetDouble()
        n[1] = normal_data[1].GetDouble()
        n[2] = normal_data[2].GetDouble()
        # Sanity check on normal
        norm2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2]
        if norm2 < 1e-12:
            raise Exception("{0} Error: something went wrong in output plane definition: plane {1} has a normal with module 0.".format(self.__class__.__name__,self.output_surface_index))

        point_data = cut_data["point"]
        p = Vector(3)
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

    def __generate_variable_list_from_input(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]

    def __generate_flags_list_from_input(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ globals()[ param[i].GetString() ] for i in range( 0,param.size() ) ]

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

    def __finalize_results(self):

        if self.body_io is not None:
            self.body_io.FinalizeResults()

        if self.cut_io is not None:
            self.cut_io.FinalizeResults()

    def __restart_list_files(self,additional_frequencies):

        self.__remove_list_files()

        self.__initialize_list_files(additional_frequencies)

        if self.post_mode == KratosMultiphysics.GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == KratosMultiphysics.GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == KratosMultiphysics.GiDPostMode.GiD_PostAsciiZipped:
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

                if( print_id != "0" ):
                    file_id.append(int(print_id))

            file_id.sort()

    def __set_multifile_lists(self,multifile_list):
        for mfilelist in multifile_list:
            self.multifilelists.append(mfilelist)

        for mfilelist in self.multifilelists:
            mfilelist.file.write("Multiple\n")
            mfilelist.index = 1

    def __write_inital_step_in_multifile_lists(self, label):
        for mfilelist in self.multifilelists:

            if (self.post_mode == KratosMultiphysics.GiDPostMode.GiD_PostBinary):
                text_to_print = self.__get_multifile_list_name(mfilelist.name)+"_"+"%.12g"%label+".post.bin\n"
                mfilelist.file.write(text_to_print)
            else:
                text_to_print1 = self.__get_multifile_list_name(mfilelist.name)+"_"+"%.12g"%label+".post.msh\n"
                text_to_print2 = self.__get_multifile_list_name(mfilelist.name)+"_"+"%.12g"%label+".post.res\n"
                mfilelist.file.write(text_to_print1)
                mfilelist.file.write(text_to_print2)
            self.Flush(mfilelist.file)

    def __write_multifile_lists(self, label):

        for mfilelist in self.multifilelists:
            if (label % mfilelist.step) == 0:

                if (self.post_mode == KratosMultiphysics.GiDPostMode.GiD_PostBinary):
                    text_to_print = self.__get_multifile_list_name(mfilelist.name)+"_"+"%.12g"%label+".post.bin\n"
                    mfilelist.file.write(text_to_print)
                else:
                    text_to_print1 = self.__get_multifile_list_name(mfilelist.name)+"_"+"%.12g"%label+".post.msh\n"
                    text_to_print2 = self.__get_multifile_list_name(mfilelist.name)+"_"+"%.12g"%label+".post.res\n"
                    mfilelist.file.write(text_to_print1)
                    mfilelist.file.write(text_to_print2)
                self.Flush(mfilelist.file)

    @classmethod
    def __get_multifile_list_name(cls, name):
        return name

    def __close_multifiles(self):
        for mfilelist in self.multifilelists:
            mfilelist.file.close()

    # NOTE (PR): 'Codacy' suggest to change the following method, from a standard method to a 'classmethod' or 'staticmethod' as it does not refer to any of the class attributes
    @classmethod
    def __remove_list_files(cls):

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
                    except WindowsError:
                        pass

    def __remove_post_results_files(self, step_label):

        path = os.getcwd()

        if self.post_mode == KratosMultiphysics.GiDPostMode.GiD_PostBinary:
            ext = ".post.bin"
        elif self.post_mode == KratosMultiphysics.GiDPostMode.GiD_PostAscii:
            ext = ".post.res"
        elif self.post_mode == KratosMultiphysics.GiDPostMode.GiD_PostAsciiZipped:
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
                    except WindowsError:
                        pass


class MultifileList(object):

    def __init__(self, name, step):
        self.index = 0
        self.step = step
        self.name = name
        absolute_path_to_file = os.path.join("_list_" + self.name + "_" + str(step) + ".post.lst")
        self.file = open(absolute_path_to_file,"w")
