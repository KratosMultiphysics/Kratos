import KratosMultiphysics as Kratos
from typing import Sequence

class ExpressionBoundingManager:
    """
    Minimal helper to map expressions between physical bounds [lo, hi]
    and a normalized phi in [0, 1].

    Provides:
    - GetBoundedExpression(x): clamps x to [lo, hi]
    - GetUnboundedExpression(phi): maps phi in [0,1] to physical domain [lo,hi]
    - GetBoundGap(): returns (hi - lo) as a scalar expression multiplier
    """
    def __init__(self, bounds_vector: 'Kratos.Vector | Sequence[float]') -> None:
        if isinstance(bounds_vector, Kratos.Vector):
            self._lo = float(bounds_vector[0])
            self._hi = float(bounds_vector[1])
        else:
            self._lo = float(bounds_vector[0])
            self._hi = float(bounds_vector[1])
        if self._hi < self._lo:
            raise RuntimeError(f"Invalid bounds: hi({self._hi}) < lo({self._lo})")

        self._gap = self._hi - self._lo

    def GetBoundedExpression(self, physical_expression: 'Kratos.Expression.Expression') -> 'Kratos.Expression.Expression':
        return Kratos.Expression.Utils.Clip(physical_expression, self._lo, self._hi)

    def GetUnboundedExpression(self, phi_expression: 'Kratos.Expression.Expression') -> 'Kratos.Expression.Expression':
        # map phi in [0,1] to physical: x = lo + phi * (hi - lo)
        return Kratos.Expression.Utils.Collapse(phi_expression * self._gap + self._lo)

    def GetBoundGap(self) -> float:
        return self._gap
