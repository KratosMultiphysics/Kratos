"""Defaults used in hdf file settings."""

hdf5_default_settings = """
            {
                "model_part_name" : "MainModelPart",
                "file_settings" : {},
                "model_part_output_settings" : {},
                "nodal_results_settings" : {},
                "element_results_settings" : {},
                "output_time_settings" : {}
            }
            """

model_part_output_default_settings = """
            {
                "prefix" : "/ModelData"
            }
            """

temporal_default_settings = """
        {
            "prefix" : "/ResultsData",
            "list_of_variables": []
        }
        """