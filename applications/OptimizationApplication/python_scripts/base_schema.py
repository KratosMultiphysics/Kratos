import pydantic, json
from pydantic import BaseModel, Field
from pydantic import ValidationError
from pydantic.types import Literal
from typing import List, Union

from KratosMultiphysics.OptimizationApplication.schema.schema_controls import ThicknessControl, Control
from KratosMultiphysics.OptimizationApplication.schema.schema_responses import Response, MassResponse, LinearStrainEnergyResponse
from KratosMultiphysics.OptimizationApplication.schema.schema_analysis import KratosAnalysisBase


class ProblemDataSchema(BaseModel):
        parallel_type: Literal["OpenMP"] = Field("OpenMP", description="Type of parallelism")
        echo_level: int = Field(0, description="Verbosity level")

class ModelPartsSchema(BaseModel):
    class Settings(BaseModel):
        model_part_name: str = Field(..., description="Name of the model part")
        domain_size: Literal[2, 3] = Field(..., description="Dimension of the model part")
        input_filename: str = Field(..., description="Name of the input file")
    settings: Settings

class BaseSchema(BaseModel):
    problem_data: ProblemDataSchema = Field(..., description="Problem data")
    model_parts: List[ModelPartsSchema] = Field(..., description="List of model parts")
    analyses: List[KratosAnalysisBase] = Field(..., description="List of responses")
    responses: List[Union[Response, MassResponse, LinearStrainEnergyResponse]] = Field(..., description="List of responses")
    controls: List[Union[ThicknessControl, Control]] = Field(..., description="List of controls")
    algorithm_settings: dict = Field(..., description="Algorithm settings")
    processes: dict = Field(..., description="List of processes")

def ValidateJSON(data):
    valid = False
    try:
        valid_data = BaseSchema(**data)
        valid = True
    except ValidationError as e:
        print(e)

    return valid, valid_data
    
    