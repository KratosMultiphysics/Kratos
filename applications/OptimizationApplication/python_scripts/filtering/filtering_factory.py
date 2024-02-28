import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import FilterTypes

def Factory(model: Kratos.Model, filtering_model_part_name: str, variable: SupportedSensitivityFieldVariableTypes, data_location: Kratos.Globals.DataLocation, settings: Kratos.Parameters) -> FilterTypes:
    if not settings.Has("filter_type"):
        raise RuntimeError(f"\"filter_type\" not specified in settings:\n{settings}")

    if not settings.Has("filtering_model_part_name"):
        settings.AddString("filtering_model_part_name", filtering_model_part_name)

    if settings["filtering_model_part_name"].GetString() != filtering_model_part_name:
        raise RuntimeError("Model part name mismatch.")

    filter_type = settings["filter_type"].GetString()

    filter_dict: 'dict[str, dict[Kratos.Globals.DataLocation, typing.Callable[[Kratos.Model, SupportedSensitivityFieldVariableTypes, Kratos.Globals.DataLocation, Kratos.Parameters], FilterTypes]]]' = {
        "explicit_vertex_morphing": {
            Kratos.Globals.DataLocation.NodeHistorical   : lambda a,_,__,d: KratosOA.NodalExplicitFilter(a, d),
            Kratos.Globals.DataLocation.NodeNonHistorical: lambda a,_,__,d: KratosOA.NodalExplicitFilter(a, d),
            Kratos.Globals.DataLocation.Condition        : lambda a,_,__,d: KratosOA.ConditionExplicitFilter(a,d),
            Kratos.Globals.DataLocation.Element          : lambda a,_,__,d: KratosOA.ElementExplicitFilter(a,d)
        }
    }

    if filter_type not in filter_dict.keys():
        raise RuntimeError(f"Unsupported filter_type=\"{filter_type}\" requested. Followings are supported:\n\t" + "\n\t".join(filter_dict.keys()))

    if data_location not in filter_dict[filter_type].keys():
        raise RuntimeError(f"Unsupported data_location=\"{data_location.name}\" for filter_type=\"{filter_type}\" requested. Followings are supported:\n\t" + "\n\t".join([v.name for v in filter_dict[filter_type].keys()]))

    return filter_dict[filter_type][data_location](model, variable, data_location, settings)
