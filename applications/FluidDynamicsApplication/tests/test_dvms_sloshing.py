"""Short fixed-parameter one- and two-phase QSVMS/DVMS sloshing comparisons."""
import math

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import dvms_gap_history
from KratosMultiphysics.FluidDynamicsApplication.two_fluid_level_set_geometry import (
    TwoFluidLevelSetTransport, extract_height_contour)
from test_one_phase_moving_gap import make_lifecycle


def run_short_case(fluid_model, fluid_element, mass_source_strategy=None):
    lifecycle = make_lifecycle(fluid_model, fluid_element=fluid_element,
                               mass_source=mass_source_strategy is not None,
                               mass_source_strategy=mass_source_strategy or "instantaneous")
    dt = 0.01
    lifecycle.dt = dt
    lifecycle.root.ProcessInfo[KM.DELTA_TIME] = dt
    lifecycle.root.ProcessInfo[KM.BDF_COEFFICIENTS] = KM.Vector([1.0/dt, -1.0/dt, 0.0])
    transport = TwoFluidLevelSetTransport(lifecycle.model, "SloshingTransport", (0., 0.), (1., 1.),
                                         (8, 8), KM.SkylineLUFactorizationSolver())
    height = lambda x: 0.62 + 0.02 * math.cos(math.pi*x)
    contour = transport.initialize(lambda x, y: y - height(x))
    lifecycle.transport = transport
    for phase, part in enumerate((lifecycle.water, lifecycle.air)):
        rho = (2.0, 0.1)[phase]
        for node in part.Nodes:
            pressure = rho * 9.81 * (height(node.X) - node.Y)
            if fluid_model == "two_phase":
                pressure += 0.1 * 9.81 * (1.0 - height(node.X))
            for slot in range(3):
                node.SetSolutionStepValue(KM.PRESSURE, slot, pressure)
    lifecycle._build_geometry(contour)
    if fluid_element == "dvms":
        for element in lifecycle.computing.Elements:
            element.Initialize(lifecycle.root.ProcessInfo)
        lifecycle._dvms_snapshot = dvms_gap_history.capture(lifecycle.computing, lifecycle.water, lifecycle.air)
    lifecycle._capture(0.0)
    settings = KM.Parameters('''{"relative_tolerance":1e-8,"absolute_tolerance":1e-10,
        "maximum_iterations":50,"maximum_free_dof_residual":1e-7,"echo_level":0}''')
    initial_area = lifecycle._mass_source_areas()[0]
    try:
        for step in range(1, 3):
            lifecycle.prepare_solution_step(step * dt)
            residual = lifecycle.solve_solution_step(KM.SkylineLUFactorizationSolver(), settings)
        velocity = max(math.hypot(node.GetSolutionStepValue(KM.VELOCITY)[0],
                                  node.GetSolutionStepValue(KM.VELOCITY)[1])
                       for node in lifecycle.computing.Nodes)
        contour = extract_height_contour(transport.part, lifecycle.tolerance)
        pressure = max(abs(node.GetSolutionStepValue(KM.PRESSURE)) for node in lifecycle.computing.Nodes)
        energy = 0.0
        for eid, point, weight in KM.FemMultipatchUtilities.GetGaussPointsWithWeights(lifecycle.computing):
            element = lifecycle.computing.GetElement(eid)
            geometry = element.GetGeometry()
            density = element.Properties[KM.DENSITY]
            a, b, _, d = geometry
            xi = (point[0] - a.X) / (b.X - a.X)
            eta = (point[1] - a.Y) / (d.Y - a.Y)
            shape = ((1-xi)*(1-eta), xi*(1-eta), xi*eta, (1-xi)*eta)
            u = [sum(shape[i]*geometry[i].GetSolutionStepValue(KM.VELOCITY)[j] for i in range(4))
                 for j in range(2)]
            energy += weight*density*(0.5*(u[0]**2+u[1]**2)+9.81*point[1])
        source = lifecycle.mass_source_diagnostics
        return {"velocity": velocity, "surface": contour[0][1] - 0.62,
                "area": lifecycle._mass_source_areas()[0] - initial_area, "residual": residual,
                "pressure": pressure, "energy": energy,
                "source_work": source.get("mass_source_cumulative_work", 0.0) if source else 0.0}
    finally:
        transport.close()


class TestDVMSSloshing(KratosUnittest.TestCase):
    def test_short_instantaneous_and_relaxed_comparisons(self):
        for model in ("one_phase", "two_phase"):
            for element in ("qsvms", "dvms"):
                results = {strategy: run_short_case(model, element, strategy)
                           for strategy in ("instantaneous", "relaxed")}
                for strategy, metrics in results.items():
                    print(model, element, strategy, metrics)
                    self.assertTrue(all(math.isfinite(value) for value in metrics.values()))
                    self.assertLess(abs(metrics["area"]), 0.05)
                    self.assertLess(metrics["residual"], 1e-7)

    def test_short_one_and_two_phase_comparisons(self):
        for model in ("one_phase", "two_phase"):
            results = {element: run_short_case(model, element) for element in ("qsvms", "dvms")}
            for element, metrics in results.items():
                print(model, element, metrics)
                self.assertTrue(all(math.isfinite(value) for value in metrics.values()))
                self.assertLess(abs(metrics["area"]), 0.05)
                self.assertLess(metrics["residual"], 1e-7)


if __name__ == "__main__":
    KratosUnittest.main()
