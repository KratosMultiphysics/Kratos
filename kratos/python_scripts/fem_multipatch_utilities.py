"""Python access to the C++ FEM multipatch geometry operations."""

import KratosMultiphysics as KM


def get_gauss_points(model_part, active_only=True):
    """Return (element ID, global coordinates) from the Kratos quadrature API.

    The C++ utility queries Geometry.IntegrationPoints using each element's
    integration method and maps them with GlobalCoordinates. Python does not
    traverse elements or approximate quadrature coordinates.
    """
    return KM.FemMultipatchUtilities.GetGaussPoints(model_part, active_only)
