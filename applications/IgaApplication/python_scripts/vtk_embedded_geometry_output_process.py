import KratosMultiphysics
from KratosMultiphysics.IgaApplication.map_nurbs_volume_results_to_embedded_geometry_process import MapNurbsVolumeResultsToEmbeddedGeometryProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return VtkEmbeddeGeometryOutputProcess(model, settings["Parameters"])


class VtkEmbeddeGeometryOutputProcess(KratosMultiphysics.OutputProcess):
    def __init__(self, model, settings):
        super().__init__()

        mapping_settings = settings["mapping_parameters"]
        vtk_settings = settings["vtk_parameters"]

        # Check for invalid variable_names.
        self.VariableCheck(vtk_settings, "nodal_flags")
        self.VariableCheck(vtk_settings, "element_data_value_variables")
        self.VariableCheck(vtk_settings, "element_flags")
        self.VariableCheck(vtk_settings, "condition_data_value_variables")
        self.VariableCheck(vtk_settings, "condition_flags")
        self.VariableCheck(vtk_settings, "gauss_point_variables_extrapolated_to_nodes")

        # Map all values to from main_model_part to embedded_model_part.
        process_params = KratosMultiphysics.Parameters(
        """ {
                "main_model_part_name"                    : "",
                "nurbs_volume_name"                       : "",
                "embedded_model_part_name"                : "",
                "nodal_results": [],
                "gauss_point_results" : []
        } """ )

        process_params["main_model_part_name"].SetString(mapping_settings["main_model_part_name"].GetString())
        process_params["nurbs_volume_name"].SetString(mapping_settings["nurbs_volume_name"].GetString())
        process_params["embedded_model_part_name"].SetString(mapping_settings["embedded_model_part_name"].GetString())
        process_params["nodal_results"].SetStringArray(vtk_settings["nodal_solution_step_data_variables"].GetStringArray())
        process_params["gauss_point_results"].SetStringArray(vtk_settings["nodal_data_value_variables"].GetStringArray())
        self.process = MapNurbsVolumeResultsToEmbeddedGeometryProcess(model, process_params)

        # Create vtk output process
        self.vtk_io = VtkOutputProcess(model, vtk_settings) # this also validates the settings

    def VariableCheck(self, settings, variable_name):
        if( settings[variable_name].size() > 0 ):
            err_msg = "Variable '" + variable_name + "' is not available. Please use nodal_solution_step_data_variables"
            err_msg += "or nodal_data_value_variables"
            raise Exception(err_msg)

    def PrintOutput(self):
        self.process.ExecuteBeforeOutputStep()
        self.vtk_io.PrintOutput()


