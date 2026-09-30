from pydantic import BaseModel, Field
from pydantic.types import Literal

class KratosAnalysisBase(BaseModel):
    name: str = Field(..., description="Name of the analysis")
    type: Literal["kratos_analysis_execution_policy"] = Field(..., description="Type of the analysis")
    class KratosAnalysisSettings(BaseModel):
        model_part_names: list[str] = Field(..., description="List of model part names")
        analysis_module: str = Field(..., description="Analysis module")
        analysis_type: str = Field(..., description="Analysis type")
        analysis_settings: dict = Field(..., description="Analysis settings")
        class AnalysisOutputSettings(BaseModel):
            nodal_solution_step_data_variables: list[str] = Field(list, description="List of nodal solution step data variables")
            nodal_data_value_variables: list[str] = Field(list, description="List of nodal data value variables")
            element_data_value_variables: list[str] = Field(list, description="List of element data value variables")
            condition_data_value_variables: list[str] = Field(list, description="List of condition data value variables")
        analysis_output_settings: AnalysisOutputSettings = Field(dict, description="Analysis output settings")
    settings: KratosAnalysisSettings
