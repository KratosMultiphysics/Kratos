from typing import Union

import KratosMultiphysics as Kratos

# used union types for type hinting in python and for automated documentation generation
SupportedSensitivityFieldVariableTypes = Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]