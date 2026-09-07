import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

def FilterRadiusFactory(model_part: Kratos.ModelPart, container_type: Kratos.Globals.DataLocation, filter_radius_settings: Kratos.Parameters) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
    if not filter_radius_settings.Has("filter_radius_type"):
        raise RuntimeError(f"\"filter_radius_type\" not specified in the following settings:\n{filter_radius_settings}")

    filter_radius_type = filter_radius_settings["filter_radius_type"].GetString()
    if filter_radius_type == "constant":
        defaults = Kratos.Parameters("""{
            "filter_radius_type": "constant",
            "filter_radius"     : 0.2
        }""")
        filter_radius_settings.ValidateAndAssignDefaults(defaults)

        if container_type == Kratos.Globals.DataLocation.NodeHistorical or container_type == Kratos.Globals.DataLocation.NodeNonHistorical:
            return Kratos.TensorAdaptors.DoubleTensorAdaptor(model_part.Nodes, Kratos.DoubleNDData([model_part.NumberOfNodes()], filter_radius_settings["filter_radius"].GetDouble()), copy=False)
        elif container_type == Kratos.Globals.DataLocation.Condition:
            return Kratos.TensorAdaptors.DoubleTensorAdaptor(model_part.Conditions, Kratos.DoubleNDData([model_part.NumberOfConditions()], filter_radius_settings["filter_radius"].GetDouble()), copy=False)
        elif container_type == Kratos.Globals.DataLocation.Element:
            return Kratos.TensorAdaptors.DoubleTensorAdaptor(model_part.Elements, Kratos.DoubleNDData([model_part.NumberOfElements()], filter_radius_settings["filter_radius"].GetDouble()), copy=False)
        else:
            raise RuntimeError(f"Unsupported container_type.")
    else:
        raise RuntimeError(f"Unsupported filter_radius_type = \"{filter_radius_type}\".")

def ExplicitFilterDampingFactory(model: Kratos.Model, data_location: Kratos.Globals.DataLocation, stride: int, parameters: Kratos.Parameters) -> 'typing.Union[KratosOA.NodeExplicitDamping, KratosOA.ConditionExplicitDamping, KratosOA.ElementExplicitDamping]':
        if not parameters.Has("damping_type"):
            raise RuntimeError(f"\"damping_type\" not specified in the following settings:\n{parameters}")

        damping_types_dict: 'dict[Kratos.Globals.DataLocation, dict[str, typing.Type[typing.Union[KratosOA.NodeExplicitDamping, KratosOA.ConditionExplicitDamping, KratosOA.ElementExplicitDamping]]]]' = {
            Kratos.Globals.DataLocation.NodeHistorical: {
                "nearest_entity": KratosOA.NearestNodeExplicitDamping,
                "integrated_nearest_entity": KratosOA.IntegratedNearestNodeExplicitDamping
            },
            Kratos.Globals.DataLocation.NodeNonHistorical: {
                "nearest_entity": KratosOA.NearestNodeExplicitDamping,
                "integrated_nearest_entity": KratosOA.IntegratedNearestNodeExplicitDamping
            },
            Kratos.Globals.DataLocation.Condition: {
                "nearest_entity": KratosOA.NearestConditionExplicitDamping,
                "integrated_nearest_entity": KratosOA.IntegratedNearestConditionExplicitDamping
            },
            Kratos.Globals.DataLocation.Element: {
                "nearest_entity": KratosOA.NearestElementExplicitDamping,
                "integrated_nearest_entity": KratosOA.IntegratedNearestElementExplicitDamping
            }
        }

        damping_type = parameters["damping_type"].GetString()
        if damping_type not in damping_types_dict[data_location].keys():
            raise RuntimeError(f"Unsupported damping_type = \"{damping_type}\". Followings are supported:\n\t" + "\n\t".join(damping_types_dict[data_location].keys()))
        return damping_types_dict[data_location][damping_type](model, parameters, stride)
