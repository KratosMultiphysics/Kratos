from KratosMultiphysics.StructuralMechanicsApplication.verification_methods.plate_buckling import PlateBuckling

_METHOD_TYPES = {
    "plate_buckling": PlateBuckling
}

def CreateVerificationMethod(name: str):
    try:
        return _METHOD_TYPES[name]()
    except KeyError as exc:
        raise ValueError(f"Unsupported verification method: {name!r}") from exc