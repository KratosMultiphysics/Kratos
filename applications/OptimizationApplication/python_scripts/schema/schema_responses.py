from pydantic import BaseModel, Field
from pydantic.types import PositiveInt, Literal


class MassResponseSettings(BaseModel):
    evaluated_model_part_names: list[str] = Field(..., description="List of model part names to evaluate")

class LinearStrainEnergyResponseSettings(BaseModel):
    evaluated_model_part_names: list[str] = Field(..., description="List of model part names to evaluate")
    primal_analysis_name: str = Field(..., description="Name of the primal analysis")
    perturbation_size: float = Field(..., description="Perturbation size for the response")

class Response(BaseModel):
    name: str = Field(..., description="Name of the response")
    type: str = Field(..., description="Type of the response")
    settings: dict


class MassResponse(Response):
    type: Literal["mass_response_function"] = Field("mass_response_function", description="Type of the response")
    settings: MassResponseSettings

class LinearStrainEnergyResponse(Response):
    type: Literal["linear_strain_energy_response_function"] = Field("linear_strain_energy_response_function", description="Type of the response")
    settings: LinearStrainEnergyResponseSettings
