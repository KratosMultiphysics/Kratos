"""Solve the one-hexahedron FE benchmark; no material-point calls or reference fit."""

import json
import math

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication  # Registers the element.
import KratosMultiphysics.ConstitutiveLawsApplication as CLA

from validation_utils import (
    COLUMNS, DIRECTORY, read_reference_parameters, read_results, write_csv,
)


def von_mises(stress):
    xx, yy, zz, xy, yz, xz = stress
    return math.sqrt(
        0.5 * ((xx - yy)**2 + (yy - zz)**2 + (zz - xx)**2)
        + 3.0 * (xy**2 + yz**2 + xz**2)
    )


def main():
    parameters = read_reference_parameters()
    read_results(DIRECTORY / "ansys/swift_validation_ansys.csv", parameters)
    material_path = DIRECTORY / "kratos/Materials.json"
    material = json.loads(material_path.read_text())["properties"][0]["Material"]
    parameter_names = {
        "YOUNG_MODULUS": "EMOD", "POISSON_RATIO": "NU",
        "SWIFT_COEFFICIENT": "KSWIFT", "SWIFT_INITIAL_STRAIN": "EPS0",
        "SWIFT_HARDENING_EXPONENT": "NSWIFT",
    }
    for variable, parameter in parameter_names.items():
        if material["Variables"][variable] != parameters[parameter]:
            raise ValueError(f"Material {variable} differs from the APDL input")

    model = KM.Model()
    model_part = model.CreateModelPart("Structure")
    model_part.SetBufferSize(2)
    model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
    model_part.ProcessInfo[KM.IS_RESTARTED] = False
    for variable in (KM.DISPLACEMENT, KM.REACTION, KM.VOLUME_ACCELERATION):
        model_part.AddNodalSolutionStepVariable(variable)

    # Exactly the APDL SOLID185 node order, also the Hexahedra3D8 local order.
    coordinates = (
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
    )
    length = parameters["LENGTH"]
    for node_id, xyz in enumerate(coordinates, 1):
        model_part.CreateNewNode(node_id, *(length * value for value in xyz))
    for displacement, reaction in (
        (KM.DISPLACEMENT_X, KM.REACTION_X),
        (KM.DISPLACEMENT_Y, KM.REACTION_Y),
        (KM.DISPLACEMENT_Z, KM.REACTION_Z),
    ):
        KM.VariableUtils().AddDof(displacement, reaction, model_part)
    settings = KM.Parameters(json.dumps({
        "Parameters": {"materials_filename": str(material_path)}
    }))
    KM.ReadMaterialsUtility(settings, model)
    element = model_part.CreateNewElement(
        "SmallDisplacementElement3D8N", 1, list(range(1, 9)),
        model_part.GetProperties()[1],
    )
    loaded_nodes = []
    for node in model_part.Nodes:
        if node.X0 == 0.0:
            node.Fix(KM.DISPLACEMENT_X)
        if node.Y0 == 0.0:
            node.Fix(KM.DISPLACEMENT_Y)
        if node.Z0 == 0.0:
            node.Fix(KM.DISPLACEMENT_Z)
        if node.X0 == length:
            node.Fix(KM.DISPLACEMENT_X)
            loaded_nodes.append(node)

    scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
    builder = KM.ResidualBasedBlockBuilderAndSolver(KM.SkylineLUFactorizationSolver())
    criterion = KM.ResidualCriteria(1e-10, 1e-12)
    criterion.SetEchoLevel(0)
    strategy = KM.ResidualBasedNewtonRaphsonStrategy(
        model_part, scheme, criterion, builder, 30, True, False, False,
    )  # Compute reactions, reuse DOF set, keep reference mesh fixed.
    strategy.SetEchoLevel(0)
    strategy.Check()
    strategy.Initialize()

    results, integration_points, diagnostics = [], [], []
    stress_names = ("sigma_xx", "sigma_yy", "sigma_zz", "sigma_xy", "sigma_yz", "sigma_xz")
    count = int(parameters["NSUB"])
    for step in range(1, count + 1):
        model_part.CloneTimeStep(step / count)
        model_part.ProcessInfo[KM.STEP] = step
        for node in loaded_nodes:
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, parameters["UMAX"] * step / count)
        strategy.InitializeSolutionStep()
        strategy.Predict()
        if not strategy.SolveSolutionStep():
            raise RuntimeError(f"FE equilibrium did not converge at substep {step}")
        strategy.FinalizeSolutionStep()  # Commit the accepted integration-point history.

        stresses = element.CalculateOnIntegrationPoints(
            KM.CAUCHY_STRESS_VECTOR, model_part.ProcessInfo,
        )
        plastic_strains = element.CalculateOnIntegrationPoints(
            CLA.ACCUMULATED_PLASTIC_STRAIN, model_part.ProcessInfo,
        )
        if len(stresses) != 8 or len(plastic_strains) != 8:
            raise RuntimeError("Expected eight integration-point stresses and histories")
        point_rows = []
        for point, (stress, p) in enumerate(zip(stresses, plastic_strains), 1):
            row = {"substep": step, "integration_point": point}
            row.update(zip(stress_names, stress))
            row.update(seqv=von_mises(stress), p_eq=p)
            if not all(math.isfinite(value) for value in row.values()):
                raise RuntimeError(f"Non-finite integration-point result at substep {step}")
            point_rows.append(row)
        integration_points.extend(point_rows)
        means = {key: math.fsum(row[key] for row in point_rows) / 8
                 for key in (*stress_names, "seqv", "p_eq")}
        spreads = {f"spread_{key}": max(row[key] for row in point_rows)
                   - min(row[key] for row in point_rows) for key in means}
        # Require homogeneity before averaging. These absolute tolerances allow
        # roundoff in the FE assembly: 1e-9 MPa for stress, 1e-14 for dimensionless p.
        for key, spread in spreads.items():
            tolerance = 1e-14 if key == "spread_p_eq" else 1e-9
            if spread > tolerance:
                raise RuntimeError(f"Nonhomogeneous {key} at substep {step}: {spread}")
        ux = model_part.GetNode(2).GetSolutionStepValue(KM.DISPLACEMENT_X)
        values = (step, ux / length, means["sigma_xx"], means["sigma_yy"],
                  means["sigma_zz"], means["seqv"], means["p_eq"], ux)
        results.append(dict(zip(COLUMNS, values)))
        diagnostics.append({
            "substep": step,
            "nonlinear_iterations": model_part.ProcessInfo[KM.NL_ITERATION_NUMBER],
            "residual_norm": model_part.ProcessInfo[KM.RESIDUAL_NORM],
            "convergence_ratio": model_part.ProcessInfo[KM.CONVERGENCE_RATIO],
            **spreads,
        })

    # Write only after all 100 actual FE solves have converged; no resampling.
    output = DIRECTORY / "kratos"
    write_csv(output / "swift_validation_kratos.csv", results)
    write_csv(output / "swift_validation_integration_points.csv", integration_points)
    write_csv(output / "swift_validation_diagnostics.csv", diagnostics)
    read_results(output / "swift_validation_kratos.csv", parameters)
    print(f"Converged {count}/{count} increments; final strain = {results[-1]['eps_xx']}")
    max_iterations = max(row["nonlinear_iterations"] for row in diagnostics)
    print(f"Maximum nonlinear iterations: {max_iterations}")
    print("Maximum integration-point spreads:")
    for key in spreads:
        print(f"  {key}: {max(row[key] for row in diagnostics):.12g}")


if __name__ == "__main__":
    main()
