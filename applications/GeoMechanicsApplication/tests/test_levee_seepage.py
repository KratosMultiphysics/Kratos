import csv
import math
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
import helper_utilities

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils

from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)


class KratosGeoMechanicsLeveeSeepageTests(KratosUnittest.TestCase):
    """Compare steady-state levee seepage with commercial FE reference profiles."""

    PHREATIC_LINE_TOLERANCE = 0.15  # m
    FLUX_RELATIVE_TOLERANCE = 0.085 # percentage
    SATURATION_TOLERANCE = 0.075

    def _read_reference(self, file_path, filename, columns):
        path = os.path.join(file_path, filename)
        self.assertTrue(os.path.isfile(path), f"Missing commercial FE reference: {path}")
        with open(path, newline="") as csv_file:
            reader = csv.DictReader(csv_file, skipinitialspace=True)
            self.assertTrue(
                set(columns).issubset(reader.fieldnames or []),
                f"{filename} must contain columns {columns}",
            )
            points = [tuple(float(row[column]) for column in columns) for row in reader]
        self.assertTrue(points, f"{filename} must contain reference values")
        self.assertTrue(
            all(math.isfinite(value) for point in points for value in point),
            f"{filename} must contain finite values",
        )
        return sorted(points)

    def _section_profile(self, model_part, nodal_values, x):
        # The case uses linear triangles: edge intersections define a piecewise
        # linear vertical profile without introducing a different triangulation.
        points = {}
        for element in model_part.Elements:
            nodes = list(element.GetGeometry())
            self.assertEqual(len(nodes), 3, "Section interpolation requires 2D3N elements")
            for first, second in zip(nodes, nodes[1:] + nodes[:1]):
                if first.X == second.X:
                    if math.isclose(first.X, x, abs_tol=1e-10, rel_tol=0.0):
                        for node in (first, second):
                            points[round(node.Y, 10)] = nodal_values[node.Id]
                    continue
                fraction = (x - first.X) / (second.X - first.X)
                if 0.0 <= fraction <= 1.0:
                    y = first.Y + fraction * (second.Y - first.Y)
                    points[round(y, 10)] = tuple(
                        a + fraction * (b - a)
                        for a, b in zip(nodal_values[first.Id], nodal_values[second.Id])
                    )
        self.assertGreaterEqual(len(points), 2, f"No vertical mesh section at x = {x}")
        return sorted(points.items())

    def _interpolate_at_y(self, profile, y):
        self.assertGreaterEqual(y, profile[0][0], f"Reference y = {y} is below the mesh")
        self.assertLessEqual(y, profile[-1][0], f"Reference y = {y} is above the mesh")
        for (lower_y, lower), (upper_y, upper) in zip(profile, profile[1:]):
            if lower_y <= y <= upper_y:
                fraction = (y - lower_y) / (upper_y - lower_y)
                return tuple(a + fraction * (b - a) for a, b in zip(lower, upper))
        self.fail(f"Cannot interpolate the reference elevation y = {y}")

    def _phreatic_elevation(self, profile, x):
        elevations = {y for y, values in profile if values[0] == 0.0}
        for (lower_y, lower), (upper_y, upper) in zip(profile, profile[1:]):
            lower_pressure, upper_pressure = lower[0], upper[0]
            self.assertFalse(
                lower_pressure == upper_pressure == 0.0,
                f"The p = 0 elevation is ambiguous at x = {x}: a whole interval is zero",
            )
            if lower_pressure * upper_pressure < 0.0:
                elevations.add(
                    lower_y - lower_pressure * (upper_y - lower_y)
                    / (upper_pressure - lower_pressure)
                )
        self.assertEqual(
            len(elevations), 1,
            f"Expected one p = 0 crossing at x = {x}, found {sorted(elevations)}",
        )
        return elevations.pop()

    def _create_profile_plot(self, file_path, filename, actual, expected, xlabel, ylabel):
        plot_utils._make_plot(
            [
                plot_utils.DataSeries(actual, label="Kratos", line_style="-", marker=""),
                plot_utils.DataSeries(
                    expected, label="Commercial FE package", line_style="", marker="x"
                ),
            ],
            os.path.join(file_path, filename),
            xlabel=xlabel,
            ylabel=ylabel,
            yaxis_inverted=False,
        )

    def test_levee_seepage(self):
        file_path = test_helper.get_file_path('test_levee_seepage')
        expected_phreatic = self._read_reference(
            file_path, "expected_phreatic_line.csv", ("X", "Y")
        )
        expected_flux = self._read_reference(
            file_path, "expected_flux_at_x_20.csv", ("Y", "Flux")
        )
        expected_saturation = self._read_reference(
            file_path, "expected_saturation_at_x_4.csv", ("Y", "Effective Saturation")
        )
        self.assertTrue(all(flux >= 0.0 for _, flux in expected_flux))
        self.assertTrue(all(0.0 <= value <= 1.0 for _, value in expected_saturation))

        with KratosUnittest.WorkFolderScope(file_path, __file__):
            with open("ProjectParameters.json", "r") as f:
                project_parameters = Kratos.Parameters(f.read())
            project = helper_utilities.run_orchestrator(project_parameters)

        model_part = project.GetModel().GetModelPart(
            "PorousDomain.porous_computational_model_part"
        )
        output_data = GiDOutputFileReader().read_output_from(
            os.path.join(file_path, "levee_seepage_find_equilibrium.post.res")
        )
        node_ids = [node.Id for node in model_part.Nodes]
        pressures, fluxes, saturations = (
            GiDOutputFileReader.nodal_values_at_time(variable, 1.0, output_data, node_ids)
            for variable in ("WATER_PRESSURE", "FLUID_FLUX_VECTOR", "EFFECTIVE_SATURATION")
        )
        nodal_values = {
            node_id: (pressure, *flux[:3], saturation)
            for node_id, pressure, flux, saturation in zip(
                node_ids, pressures, fluxes, saturations
            )
        }
        phreatic = [
            (x, self._phreatic_elevation(
                self._section_profile(model_part, nodal_values, x), x
            ))
            for x, _ in expected_phreatic
        ]
        flux_profile = self._section_profile(model_part, nodal_values, 20.0)
        saturation_profile = self._section_profile(model_part, nodal_values, 4.0)
        actual_flux = [
            (y, math.sqrt(sum(
                v * v for v in self._interpolate_at_y(flux_profile, y)[1:4]
            )) * 86400)
            for y, _ in expected_flux
        ]
        actual_saturation = [
            (y, self._interpolate_at_y(saturation_profile, y)[4])
            for y, _ in expected_saturation
        ]

        if test_helper.want_test_plots():
            self._create_profile_plot(
                file_path, "phreatic_line.svg", phreatic, expected_phreatic, "X [m]", "Y [m]"
            )
            self._create_profile_plot(
                file_path, "fluxes_for_x_20.svg",
                [
                    (math.sqrt(sum(v * v for v in values[1:4])) * 86400, y)
                    for y, values in flux_profile
                ],
                [(flux, y) for y, flux in expected_flux], "|Flux| [m/d]", "Y [m]",
            )
            self._create_profile_plot(
                file_path, "saturation_for_x_4.svg",
                [(values[4], y) for y, values in saturation_profile],
                [(saturation, y) for y, saturation in expected_saturation],
                "Effective saturation [-]", "Y [m]",
            )

        for (x, actual), (_, expected) in zip(phreatic, expected_phreatic):
            with self.subTest(profile="phreatic line", x=x):
                self.assertAlmostEqual(actual, expected, delta=self.PHREATIC_LINE_TOLERANCE)
        for (y, actual), (_, expected) in zip(actual_flux, expected_flux):
            with self.subTest(profile="flux magnitude", x=20.0, y=y):
                self.assertAlmostEqual(
                    actual, expected, delta=self.FLUX_RELATIVE_TOLERANCE * abs(expected)
                )
        for (y, actual), (_, expected) in zip(actual_saturation, expected_saturation):
            with self.subTest(profile="effective saturation", x=4.0, y=y):
                self.assertAlmostEqual(actual, expected, delta=self.SATURATION_TOLERANCE)


if __name__ == '__main__':
    KratosUnittest.main()