"""Isolated scalar Gap-SBM harmonic velocity extension on a Cartesian grid."""
import bisect
import math
import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFD
from KratosMultiphysics import python_linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication.two_fluid_level_set_geometry import create_complete_phase_view
from KratosMultiphysics.FluidDynamicsApplication.two_fluid_phase_history import PhaseHistory


class WaterHarmonicVelocityExtension:
    """Solve in air and optionally continue vertically above a horizontal top.

    The artificial top has zero normal derivative. Its position changes the
    harmonic solution below it; increasing the layer count moves it away.
    """

    def __init__(self, part, contour, water_sampler, water_patch, tolerance,
                 solid_walls=("left", "right", "bottom"), harmonic_extension_layers=None):
        if (harmonic_extension_layers is not None and
                (isinstance(harmonic_extension_layers, bool) or not isinstance(harmonic_extension_layers, int)
                 or harmonic_extension_layers <= 0)):
            raise ValueError("harmonic_extension_layers must be None or a positive integer")
        # Import registers the Laplace element only when the strategy is used.
        import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusion

        self.sampler = water_sampler
        self.model = KM.Model()
        root = self.model.CreateModelPart("HarmonicExtension")
        for variable in (KM.TEMPERATURE, KM.CONDUCTIVITY, KM.HEAT_FLUX, KM.VELOCITY):
            root.AddNodalSolutionStepVariable(variable)
        root.ProcessInfo[KM.DOMAIN_SIZE] = 2
        settings = KM.ConvectionDiffusionSettings()
        settings.SetUnknownVariable(KM.TEMPERATURE)
        settings.SetDiffusionVariable(KM.CONDUCTIVITY)
        settings.SetVolumeSourceVariable(KM.HEAT_FLUX)
        root.ProcessInfo[KM.CONVECTION_DIFFUSION_SETTINGS] = settings
        properties = root.CreateNewProperties(1)
        properties[KM.PENALTY_COEFFICIENT] = 40.0
        nodes = list(part.Nodes)
        xs, ys = sorted({node.X for node in nodes}), sorted({node.Y for node in nodes})
        if len(nodes) != len(xs) * len(ys):
            raise ValueError("Harmonic extension requires a complete Cartesian level-set grid")
        top_index = len(ys) - 1
        bottom_index = 0
        if harmonic_extension_layers is not None:
            y_max = max(point[1] for point in contour)
            spacings = [b - a for a, b in zip(ys, ys[1:])]
            if not spacings:
                raise ValueError("Harmonic extension requires at least two grid rows")
            if all(math.isclose(h, spacings[0], rel_tol=1e-10, abs_tol=tolerance) for h in spacings):
                target = y_max + harmonic_extension_layers * spacings[0]
                top_index = min(bisect.bisect_left(ys, target - tolerance), len(ys) - 1)
            else:
                first = bisect.bisect_left(ys, y_max - tolerance)
                top_index = min(first + harmonic_extension_layers, len(ys) - 1)
            if top_index < len(ys) - 1:
                bottom_index = max(0, bisect.bisect_right(ys, min(point[1] for point in contour) + tolerance) - 1)
        self.top = ys[top_index]
        auxiliary_ys = ys[bottom_index:top_index + 1]
        by_coordinate = {(node.X, node.Y): node for node in nodes}
        auxiliary_nodes = [node for node in nodes if auxiliary_ys[0] <= node.Y <= self.top]
        for node in auxiliary_nodes:
            auxiliary = root.CreateNewNode(node.Id, node.X, node.Y, node.Z)
            auxiliary.AddDof(KM.TEMPERATURE)
            auxiliary.SetSolutionStepValue(KM.CONDUCTIVITY, 1.0)
            auxiliary.SetSolutionStepValue(KM.HEAT_FLUX, 0.0)
        patch = root.CreateSubModelPart("Grid")
        patch.AddNodes([node.Id for node in auxiliary_nodes])
        for j in range(len(auxiliary_ys) - 1):
            for i in range(len(xs) - 1):
                ids = [by_coordinate[(x, y)].Id for x, y in ((xs[i], auxiliary_ys[j]), (xs[i+1], auxiliary_ys[j]),
                       (xs[i+1], auxiliary_ys[j+1]), (xs[i], auxiliary_ys[j+1]))]
                patch.CreateNewElement("LaplacianElement2D4N", j * (len(xs)-1) + i + 1, ids, properties)
        phase = create_complete_phase_view(patch, "Air", contour, False, tolerance)
        matrix = KM.Matrix(len(contour), 2)
        for i, point in enumerate(contour):
            matrix[i, 0], matrix[i, 1] = point
        computing = root.CreateSubModelPart("Computing")
        ConvectionDiffusion.LaplacianGapGeometryUtility.Create(phase, computing, matrix, properties, tolerance)
        gap_geometry = root.CreateSubModelPart("GapGeometry")
        KM.FemMultipatchGapUtilities.CreateLevelSetPhaseGeometry(phase, gap_geometry, matrix, False, tolerance)
        interface = computing.GetSubModelPart("PhysicalInterfaces")
        self.values = {}
        self._boundary_conditions = tuple(interface.Conditions)
        boundary_points = [condition.CalculateOnIntegrationPoints(KM.INTEGRATION_COORDINATES, root.ProcessInfo)
                           for condition in self._boundary_conditions]
        traces = [[tuple(water_sampler(point[0], point[1], 0)[:2]) for point in points]
                  for points in boundary_points]
        walls = {"left": (0, xs[0]), "right": (0, xs[-1]),
                 "bottom": (1, ys[0]), "top": (1, ys[-1])}
        if any(wall not in walls for wall in solid_walls):
            raise ValueError("Unknown harmonic extension solid wall")
        fixed = [walls[wall] for wall in solid_walls]
        solver = python_linear_solver_factory.ConstructSolver(
            KM.Parameters('{"solver_type":"amgcl", "tolerance":1e-10, "max_iteration":1000}'))
        components = []
        for component in range(2):
            for node in computing.Nodes:
                node.Free(KM.TEMPERATURE)
                node.SetSolutionStepValue(KM.TEMPERATURE, 0.0)
                if any(axis == component and abs((node.X, node.Y)[axis] - coordinate) <= tolerance
                       for axis, coordinate in fixed):
                    node.Fix(KM.TEMPERATURE)
            for condition, values in zip(self._boundary_conditions, traces):
                condition.SetValue(ConvectionDiffusion.GAP_DIRICHLET_VALUES,
                                   KM.Vector([value[component] for value in values]))
            strategy = KM.ResidualBasedLinearStrategy(computing, KM.ResidualBasedIncrementalUpdateStaticScheme(),
                KM.ResidualBasedBlockBuilderAndSolver(solver), False, False, False, False)
            strategy.SetEchoLevel(0)
            try:
                strategy.Initialize()
                strategy.Check()
                strategy.Solve()
                components.append({node.Id: node.GetSolutionStepValue(KM.TEMPERATURE) for node in computing.Nodes})
            finally:
                strategy.Clear()
        for node in computing.Nodes:
            node.SetSolutionStepValue(KM.VELOCITY, [components[0][node.Id], components[1][node.Id], 0.0])
        history = PhaseHistory("air extension")
        history.capture(phase, gap_geometry, 0.0)
        air_sampler = CFD.TwoFluidPhaseVelocitySampler(history.slots[0][1], "air extension", tolerance)
        for node in nodes:
            if node.GetSolutionStepValue(KM.DISTANCE) > tolerance:
                self.values[(node.X, node.Y)] = tuple(air_sampler.Sample(node.X, min(node.Y, self.top))[:2])
            elif abs(node.GetSolutionStepValue(KM.DISTANCE)) <= tolerance:
                self.values[(node.X, node.Y)] = tuple(water_sampler(node.X, node.Y, 0)[:2])

    def __call__(self, x, y, phase):
        if phase == 0:
            return self.sampler(x, y, 0)[:2]
        return self.values[(x, y)]

    def diagnostics(self, part, tolerance):
        return {"strategy": "water_harmonic_extension", "samples": None}
