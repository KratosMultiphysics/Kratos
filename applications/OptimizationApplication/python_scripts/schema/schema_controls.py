from pydantic import BaseModel, Field

class FilterSettings(BaseModel):
    filter_type: str = Field(..., description="Type of filter")
    filter_radius: float = Field(..., description="Radius of the filter")

class Control(BaseModel):
    name: str = Field(..., description="Name of the control")
    type: str = Field(..., description="Type of the control")
    settings: dict

class ThicknessControlSettings(BaseModel):
    controlled_model_part_names: list[str] = Field(..., description="List of model part names to control")
    filter_settings: FilterSettings = Field(..., description="Filter settings")
    output_all_fields: bool = Field(..., description="Flag to output all fields")
    physical_thicknesses: list[float] = Field(..., description="List of physical thicknesses")
    class ProjectionSettings(BaseModel):
        type: str = Field(..., description="Type of projection")
        initial_value: float = Field(..., description="Initial value for projection")
        max_value: float = Field(..., description="Maximum value for projection")
        increase_fac: float = Field(..., description="Increase factor for projection")
        update_period: int = Field(..., description="Update period for projection")
        penalty_factor: float = Field(..., description="Penalty factor for projection")
    thickness_projection_settings: ProjectionSettings

class ThicknessControl(Control):
    settings: ThicknessControlSettings


