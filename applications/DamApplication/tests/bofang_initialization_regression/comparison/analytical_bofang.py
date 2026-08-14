import math
import sys


class BofangAnalytical:
    """Analytical Bofang reservoir temperature, mirroring the C++ implementation.

    References: applications/DamApplication/custom_processes/dam_bofang_condition_temperature_process.hpp
    Formula (for nodes with ``aux = water_level - coordinate[direction] >= 0``)::

        aux1 = (bottom_temp - surface_temp*exp(-0.04*height)) / (1 - exp(-0.04*height))
        T    = aux1 + (surface_temp - aux1)*exp(-0.04*aux)
              + amplitude*exp(-0.018*aux)
                * cos(freq*(month - day_max_temp/30.0 - 2.15
                            + 1.30*exp(-0.085*aux)))

    with ``freq = 0.52323``.
    """

    FREQ = 0.52323

    def __init__(self, surface_temp=10.0, bottom_temp=5.0, height=20.0,
                 amplitude=3.0, day_max_temp=15, water_level=14.0, month=6.5,
                 gravity_direction="Y"):
        self.surface_temp = surface_temp
        self.bottom_temp = bottom_temp
        self.height = height
        self.amplitude = amplitude
        self.day_max_temp = day_max_temp
        self.water_level = water_level
        self.month = month
        self.gravity_direction = gravity_direction

    def _coordinate(self, node):
        if self.gravity_direction == "X":
            return node["x"]
        elif self.gravity_direction == "Y":
            return node["y"]
        else:
            return node["z"]

    def temperature(self, node):
        coord = self._coordinate(node)
        aux = self.water_level - coord
        if aux < 0.0:
            return None  # above water level: the process does not assign T
        aux1 = (self.bottom_temp - self.surface_temp * math.exp(-0.04 * self.height)) / (
            1 - math.exp(-0.04 * self.height))
        return (
            aux1
            + (self.surface_temp - aux1) * math.exp(-0.04 * aux)
            + self.amplitude * math.exp(-0.018 * aux)
            * math.cos(self.FREQ * (self.month - self.day_max_temp / 30.0 - 2.15
                                   + 1.30 * math.exp(-0.085 * aux)))
        )


def main():
    analytical = BofangAnalytical()

    # Face nodes of the case bofang_small.mdpa (upstream face, X = 0)
    nodes = [
        {"id": 1, "x": 0.0, "y": 20.0},
        {"id": 4, "x": 0.0, "y": 16.0},
        {"id": 7, "x": 0.0, "y": 12.0},
        {"id": 10, "x": 0.0, "y": 8.0},
        {"id": 13, "x": 0.0, "y": 4.0},
        {"id": 16, "x": 0.0, "y": 0.0},
    ]

    print("node_id, elevation, depth, analytical_temperature")
    for node in nodes:
        t = analytical.temperature(node)
        depth = analytical.water_level - node["y"] if t is not None else None
        t_str = "%.10f" % t if t is not None else "n/a (above water level)"
        depth_str = "%.4f" % depth if depth is not None else "n/a"
        print("%d, %.4f, %s, %s" % (node["id"], node["y"], depth_str, t_str))


if __name__ == "__main__":
    sys.exit(main())
