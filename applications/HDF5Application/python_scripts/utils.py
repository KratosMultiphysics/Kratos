'''HDF5 process utilities.

license: HDF5Application/license.txt
'''
import KratosMultiphysics


def IsDistributed():
    return KratosMultiphysics.DataCommunicator.GetDefault().IsDistributed()


def InsertSingleSetting(dest, key, parameters):
    if dest.Has(key):
        dest[key] = parameters
    else:
        dest.AddValue(key, parameters)


def InsertSettings(src, dest):
    for key, parameters in src.items():
        InsertSingleSetting(dest, key, parameters)


def InsertArrayOfSettings(src, dest):
    for current_src, current_dest in zip(src, dest):
        InsertSettings(current_src, current_dest)


def CheckForDeprecatedFilename(deprecated_settings, module_name, *new_io_settings):
    if deprecated_settings.Has("file_name"):
        depr_msg = '\nDEPRECATION-WARNING: "file_name" should be specified in "file_settings". This will be removed in the future!\n'
        KratosMultiphysics.Logger.PrintWarning(module_name, depr_msg)
        file_name = deprecated_settings["file_name"].GetString().replace(
            ".h5", "")
        for io_settings in new_io_settings:
            new_file_name = io_settings["file_name"].GetString().replace(
                "<identifier>", file_name)
            io_settings["file_name"].SetString(new_file_name)


def CheckForDeprecatedTemporalSettings(deprecated_settings, module_name, *new_controller_settings):
    if deprecated_settings.Has("time_tag_precision"):
        depr_msg = '\nDEPRECATION-WARNING: "time_tag_precision" is ignored, please use "time_format" in "file_settings"!\n'
        KratosMultiphysics.Logger.PrintWarning(module_name, depr_msg)

    if deprecated_settings.Has("output_time_frequency"):
        depr_msg = '\nDEPRECATION-WARNING: "output_time_frequency" is deprecated, please use "time_frequency"!\n'
        KratosMultiphysics.Logger.PrintWarning(module_name, depr_msg)
        for controller_settings in new_controller_settings:
            InsertSingleSetting(controller_settings, "time_frequency",
                                deprecated_settings["output_time_frequency"])

    if deprecated_settings.Has("output_step_frequency"):
        depr_msg = '\nDEPRECATION-WARNING: "output_step_frequency" is deprecated, please use "step_frequency"!\n'
        KratosMultiphysics.Logger.PrintWarning(module_name, depr_msg)
        for controller_settings in new_controller_settings:
            InsertSingleSetting(controller_settings, "step_frequency",
                                deprecated_settings["output_step_frequency"])
