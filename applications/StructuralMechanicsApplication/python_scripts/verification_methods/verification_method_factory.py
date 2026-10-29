from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.panel_sjtu_puck import PuckFirstPlyFailure
from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.plate_buckling import PlateBuckling

_METHOD_TYPES = {
    "plate_buckling": PlateBuckling,
    "puck": PuckFirstPlyFailure
}

def CreateVerificationMethod(name: str):
    try:
        return _METHOD_TYPES[name]()
    except KeyError as exc:
        raise ValueError(f"Unsupported verification method: {name!r}") from exc