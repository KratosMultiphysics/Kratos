from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.panel_sjtu_buckling import PanelBuckling
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.panel_sjtu_puck import PuckFirstPlyFailure

_METHOD_TYPES = {
    "panel_buckling": PanelBuckling,
    "puck": PuckFirstPlyFailure
}

def CreateHandbookMethod(name: str):
    try:
        return _METHOD_TYPES[name]()
    except KeyError as exc:
        raise ValueError(f"Unsupported handbook method: {name!r}") from exc