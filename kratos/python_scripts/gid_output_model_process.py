from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
from KratosMultiphysics import *
CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return GiDOutputModelProcess(Model, settings["Parameters"])

import gid_output_process
class GiDOutputModelProcess(gid_output_process.GiDOutputProcess):

    defaults = Parameters('''{
        "filename" : "Kratos_GiD_output",
        "model_part_name" : "",
        "result_file_configuration": {
            "gidpost_flags": {
                "GiDPostMode": "GiD_PostBinary",
                "WriteDeformedMeshFlag": "WriteUndeformed",
                "WriteConditionsFlag": "WriteElementsOnly",
                "MultiFileFlag": "SingleFile"
            },
            "file_label": "time",
            "output_control_type": "step",
            "output_frequency": 1.0,
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

    default_plane_output_data = Parameters('''{
        "normal": [0.0, 0.0, 0.0],
        "point" : [0.0, 0.0, 0.0]
    }''')

    __post_mode = {
        # JSON input
        "GiD_PostAscii":        GiDPostMode.GiD_PostAscii,
        "GiD_PostAsciiZipped":  GiDPostMode.GiD_PostAsciiZipped,
        "GiD_PostBinary":       GiDPostMode.GiD_PostBinary,
        "GiD_PostHDF5":         GiDPostMode.GiD_PostHDF5,
        # Legacy
        "Binary":               GiDPostMode.GiD_PostBinary,
        "Ascii":                GiDPostMode.GiD_PostAscii,
        "AsciiZipped":          GiDPostMode.GiD_PostAsciiZipped,
    }

    __write_deformed_mesh = {
        # JSON input
        "WriteDeformed":        WriteDeformedMeshFlag.WriteDeformed,
        "WriteUndeformed":      WriteDeformedMeshFlag.WriteUndeformed,
        # Legacy
        True:                   WriteDeformedMeshFlag.WriteDeformed,
        False:                  WriteDeformedMeshFlag.WriteUndeformed,
    }

    __write_conditions = {
        # JSON input
        "WriteConditions":      WriteConditionsFlag.WriteConditions,
        "WriteElementsOnly":    WriteConditionsFlag.WriteElementsOnly,
        "WriteConditionsOnly":  WriteConditionsFlag.WriteConditionsOnly,
        True:                   WriteConditionsFlag.WriteConditions,
        False:                  WriteConditionsFlag.WriteElementsOnly,
    }
        # Legacy

    __multi_file_flag = {
        # JSON input
        "SingleFile":           MultiFileFlag.SingleFile,
        "MultipleFiles":        MultiFileFlag.MultipleFiles,
        # Legacy
        "Multiples":            MultiFileFlag.MultipleFiles,
        "Single":               MultiFileFlag.SingleFile,
    }

    def __init__(self, Model, param = None):
        if param is None:
            param = self.defaults
        else:
            # Note: this only validates the first level of the JSON tree.
            # I'm not going for recursive validation because some branches may
            # not exist and I don't want the validator assinging defaults there.
            param.ValidateAndAssignDefaults(self.defaults)

        self.param = param
        self.base_file_name = self.param["filename"].GetString()

        self.model_part = Model[param["model_part_name"].GetString()]
        self.body_io = None
        self.volume_list_files = []

        # The following are only used if we asked to print results on surfaces
        self.cut_model_part = None
        self.cut_io = None
        self.output_surface_index = 0
        self.cut_list_files = []

        point_data_configuration = self.param["point_data_configuration"]
        if point_data_configuration.size() > 0:
            import point_output_process
            self.point_output_process = point_output_process.PointOutputProcess(self.model_part, point_data_configuration)
        else:
            self.point_output_process = None

        self.step_count = 0
        self.printed_step_count = 0
        self.next_output = 0.0

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
            self.cut_model_part = ModelPart("CutPart")
            self.cut_manager = CuttingUtility()
            self._initialize_cut_output(plane_output_configuration)

        # Retrieve gidpost flags and setup GiD output tool
        gidpost_flags = result_file_configuration["gidpost_flags"]
        gidpost_flags.ValidateAndAssignDefaults(self.defaults["result_file_configuration"]["gidpost_flags"])

        self._InitializeGiDIO(gidpost_flags,gidpost_flags)

        # Process nodal and gauss point output
        self.nodal_variables = self._GenerateVariableListFromInput(result_file_configuration["nodal_results"])
        self.gauss_point_variables = self._GenerateVariableListFromInput(result_file_configuration["gauss_point_results"])
        self.nodal_nonhistorical_variables = self._GenerateVariableListFromInput(result_file_configuration["nodal_nonhistorical_results"])
        self.nodal_flags = self._GenerateFlagsListFromInput(result_file_configuration["nodal_flags_results"])
        self.nodal_flags_names =[]
        for i in range(result_file_configuration["nodal_flags_results"].size()):
            self.nodal_flags_names.append(result_file_configuration["nodal_flags_results"][i].GetString())
        self.elemental_conditional_flags = self._GenerateFlagsListFromInput(result_file_configuration["elemental_conditional_flags_results"])
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

        output_control_type = result_file_configuration["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_is_time = True
        elif output_control_type == "step":
            self.output_control_is_time = False
        else:
            msg = "{0} Error: Unknown value \"{1}\" read for parameter \"{2}\"".format(self.__class__.__name__,output_file_label,"file_label")
            raise Exception(msg)

        self.output_frequency = result_file_configuration["output_frequency"].GetDouble()

        # get .post.lst files
        additional_list_file_data = result_file_configuration["additional_list_files"]
        additional_list_files = [ additional_list_file_data[i].GetInt() for i in range(0,additional_list_file_data.size()) ]

        # Set current time parameters
        if(  self.model_part.ProcessInfo[IS_RESTARTED] == True ):
            self.step_count = self.model_part.ProcessInfo[STEP]
            self.printed_step_count = self.model_part.ProcessInfo[PRINTED_STEP]

            if self.output_control_is_time:
                self.next_output = self.model_part.ProcessInfo[TIME]
            else:
                self.next_output = self.model_part.ProcessInfo[STEP]

                # Remove post results
            if self.output_label_is_time:
                label = self.model_part.ProcessInfo[TIME]
            else:
                label = self.printed_step_count

            self._remove_post_results_files(label)

            # Restart .post.lst files
            self._restart_list_files(additional_list_files)
        else:
            # Create .post.lst files
            self._initialize_list_files(additional_list_files)

        # Process point recording data
        if self.point_output_process is not None:
            self.point_output_process.ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):
        if self.IsOutputStep():
            self.PrintOutput()

        if self.point_output_process is not None:
            self.point_output_process.ExecuteFinalizeSolutionStep()