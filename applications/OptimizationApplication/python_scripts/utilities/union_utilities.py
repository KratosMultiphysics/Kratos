from typing import Union

import KratosMultiphysics as Kratos

SupportedControlVariableTypes = Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]
SupportedSensitivityVariableTypes = Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]