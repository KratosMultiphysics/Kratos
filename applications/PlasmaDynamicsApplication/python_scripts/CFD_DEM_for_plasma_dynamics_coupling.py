from KratosMultiphysics.SwimmingDEMApplication.CFD_DEM_coupling import ProjectionModule as SDEMProjectionModule


# TODO: 
# Complete the class
class ProjectionModule(SDEMProjectionModule):

    def __init__(self,
                fluid_model_part,
                balls_model_part,
                FEM_DEM_model_part,
                project_parameters,
                coupling_dem_vars,
                coupling_fluid_vars,
                time_filtered_vars,
                flow_field=None,
                domain_size=3):

        super().__init__(fluid_model_part,
                        balls_model_part,
                        FEM_DEM_model_part,
                        project_parameters,
                        coupling_dem_vars,
                        coupling_fluid_vars,
                        time_filtered_vars,
                        flow_field,
                        domain_size)

	
