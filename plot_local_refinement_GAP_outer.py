from __future__ import annotations

import argparse
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt

import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication  # noqa: F401
import KratosMultiphysics.StructuralMechanicsApplication  # noqa: F401


def run_modelers(current_model: KM.Model, modelers_list: KM.Parameters) -> None:
    from KratosMultiphysics.modeler_factory import KratosModelerFactory

    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()


def load_project_parameters(project_parameters_path: Path) -> KM.Parameters:
    return KM.Parameters(project_parameters_path.read_text())


def get_root_model_part(model: KM.Model, parameters: KM.Parameters) -> KM.ModelPart:
    solver_settings = parameters["solver_settings"]
    model_part_name = solver_settings["model_part_name"].GetString()
    domain_size = solver_settings["domain_size"].GetInt()
    buffer_size = solver_settings["buffer_size"].GetInt()

    root_model_part = (
        model.GetModelPart(model_part_name)
        if model.HasModelPart(model_part_name)
        else model.CreateModelPart(model_part_name)
    )
    root_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)
    root_model_part.SetBufferSize(buffer_size)
    root_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    root_model_part.AddNodalSolutionStepVariable(KM.REACTION)
    return root_model_part


def find_modeler_parameters(
    parameters: KM.Parameters,
    modeler_name: str,
) -> KM.Parameters:
    modelers = parameters["modelers"]
    for i in range(modelers.size()):
        if modelers[i].Has("modeler_name") and modelers[i]["modeler_name"].GetString() == modeler_name:
            return modelers[i]["Parameters"]
    raise RuntimeError(f"Modeler '{modeler_name}' was not found in the project parameters.")


def get_patch_prefix(local_refinement_parameters: KM.Parameters) -> str:
    return (
        local_refinement_parameters["child_patch_prefix"].GetString()
        if local_refinement_parameters.Has("child_patch_prefix")
        else "Patch"
    )


def get_skin_root_name(local_refinement_parameters: KM.Parameters) -> str:
    geometry_parameters = local_refinement_parameters["geometry_parameters"]
    if geometry_parameters.Has("skin_model_part_name"):
        candidate_name = geometry_parameters["skin_model_part_name"].GetString()
        if candidate_name:
            return candidate_name
    return local_refinement_parameters["model_part_name"].GetString()


def get_element_submodel_part_suffixes(local_refinement_parameters: KM.Parameters) -> list[str]:
    suffixes: list[str] = []
    if not local_refinement_parameters.Has("analysis_parameters"):
        return suffixes

    analysis_parameters = local_refinement_parameters["analysis_parameters"]
    if not analysis_parameters.Has("element_condition_list"):
        return suffixes

    element_condition_list = analysis_parameters["element_condition_list"]
    for i in range(element_condition_list.size()):
        item = element_condition_list[i]
        if not item.Has("type") or item["type"].GetString() != "element":
            continue
        if item.Has("iga_model_part"):
            suffixes.append(item["iga_model_part"].GetString())
    return suffixes


def get_patch_sort_key(patch_name: str, patch_prefix: str) -> tuple[int, str]:
    suffix = patch_name[len(patch_prefix):]
    if suffix.isdigit():
        return int(suffix), patch_name
    return 10**9, patch_name


def iter_patch_model_parts(
    model: KM.Model,
    root_model_part_name: str,
    patch_prefix: str,
) -> Iterable[KM.ModelPart]:
    root_model_part = model.GetModelPart(root_model_part_name)
    patch_sub_model_parts = [
        root_model_part.GetSubModelPart(sub_model_part.Name())
        for sub_model_part in root_model_part.SubModelParts
        if sub_model_part.Name().startswith(patch_prefix)
    ]
    patch_sub_model_parts.sort(key=lambda mp: get_patch_sort_key(mp.Name(), patch_prefix))
    return patch_sub_model_parts


def get_patch_skin_model_part_name(skin_root_name: str, patch_name: str) -> str:
    return f"{skin_root_name}_{patch_name}_skin"


def get_patch_box_coordinates(patch_model_part: KM.ModelPart) -> tuple[float, float, float, float]:
    if patch_model_part.NumberOfNodes() == 0:
        raise RuntimeError(f"Patch model part '{patch_model_part.FullName()}' has no nodes.")

    xs = [node.X for node in patch_model_part.Nodes]
    ys = [node.Y for node in patch_model_part.Nodes]
    return min(xs), max(xs), min(ys), max(ys)


def plot_patch_box(
    ax: plt.Axes,
    patch_model_part: KM.ModelPart,
    color: tuple[float, float, float, float],
    label: str,
) -> None:
    x_min, x_max, y_min, y_max = get_patch_box_coordinates(patch_model_part)
    ax.plot(
        [x_min, x_max, x_max, x_min, x_min],
        [y_min, y_min, y_max, y_max, y_min],
        color=color,
        linewidth=2.0,
        label=label,
    )
    ax.text(
        0.5 * (x_min + x_max),
        0.5 * (y_min + y_max),
        patch_model_part.Name(),
        color=color,
        fontsize=10,
        ha="center",
        va="center",
        bbox={"facecolor": "white", "edgecolor": color, "alpha": 0.8, "boxstyle": "round,pad=0.2"},
    )


def plot_skin_loop_segments(
    ax: plt.Axes,
    skin_model_part: KM.ModelPart,
    loop_name: str,
    color: tuple[float, float, float, float],
    label: str | None,
) -> None:
    if not skin_model_part.HasSubModelPart(loop_name):
        return

    loop_model_part = skin_model_part.GetSubModelPart(loop_name)
    first_segment = True
    for condition in loop_model_part.Conditions:
        geometry = condition.GetGeometry()
        ax.plot(
            [geometry[0].X, geometry[1].X],
            [geometry[0].Y, geometry[1].Y],
            color=color,
            linestyle="--" if loop_name == "inner" else "-.",
            linewidth=1.2,
            alpha=0.9,
            label=label if first_segment else None,
        )
        first_segment = False


def plot_element_centers(
    ax: plt.Axes,
    element_model_part: KM.ModelPart,
    color: tuple[float, float, float, float],
    marker: str,
    label: str,
    size: float,
) -> None:
    if element_model_part.NumberOfElements() == 0:
        return

    xs = []
    ys = []
    for element in element_model_part.Elements:
        center = element.GetGeometry().Center()
        xs.append(center[0])
        ys.append(center[1])

    ax.scatter(xs, ys, s=size, marker=marker, color=color, label=label, zorder=3)


def deduplicate_legend(ax: plt.Axes) -> None:
    handles, labels = ax.get_legend_handles_labels()
    unique: dict[str, object] = {}
    for handle, label in zip(handles, labels):
        if label and label not in unique:
            unique[label] = handle
    ax.legend(unique.values(), unique.keys(), loc="best")


@contextmanager
def pushd(directory: Path):
    previous_directory = Path.cwd()
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(previous_directory)


def build_and_plot(project_parameters_path: Path, output_path: Path | None, show: bool) -> None:
    parameters = load_project_parameters(project_parameters_path)
    local_refinement_parameters = find_modeler_parameters(parameters, "LocalRefinementModeler")

    model = KM.Model()
    root_model_part = get_root_model_part(model, parameters)

    with pushd(project_parameters_path.parent):
        run_modelers(model, parameters["modelers"])

    patch_prefix = get_patch_prefix(local_refinement_parameters)
    skin_root_name = get_skin_root_name(local_refinement_parameters)
    element_suffixes = get_element_submodel_part_suffixes(local_refinement_parameters)

    fig, ax = plt.subplots(figsize=(12, 10))
    cmap = plt.get_cmap("tab10")

    for patch_index, patch_model_part in enumerate(
        iter_patch_model_parts(model, root_model_part.FullName(), patch_prefix)
    ):
        color = cmap(patch_index % 10)
        patch_name = patch_model_part.Name()
        patch_full_name = patch_model_part.FullName()

        plot_patch_box(ax, patch_model_part, color, f"{patch_name} border")

        skin_model_part_name = get_patch_skin_model_part_name(skin_root_name, patch_name)
        if model.HasModelPart(skin_model_part_name):
            skin_model_part = model.GetModelPart(skin_model_part_name)
            plot_skin_loop_segments(ax, skin_model_part, "outer", color, f"{patch_name} outer skin")
            plot_skin_loop_segments(ax, skin_model_part, "inner", color, f"{patch_name} inner skin")

        for element_suffix in element_suffixes:
            full_name = f"{patch_full_name}.{element_suffix}"
            if model.HasModelPart(full_name):
                plot_element_centers(
                    ax,
                    model.GetModelPart(full_name),
                    color,
                    marker="o",
                    label=f"{patch_name} {element_suffix} centers",
                    size=18.0,
                )

        gap_elements_name = f"{patch_full_name}.GapElements"
        if model.HasModelPart(gap_elements_name):
            plot_element_centers(
                ax,
                model.GetModelPart(gap_elements_name),
                color,
                marker="x",
                label=f"{patch_name} GapElements centers",
                size=28.0,
            )

    ax.set_title("Local refinement GAP-SBM patches and element centers")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linestyle=":", linewidth=0.7)
    deduplicate_legend(ax)
    fig.tight_layout()

    if output_path is None:
        output_path = project_parameters_path.with_name(
            f"{project_parameters_path.stem}_patch_plot.png"
        )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {output_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def parse_arguments() -> argparse.Namespace:
    default_project_parameters = Path(__file__).with_name("ProjectParameters_local_refinement_half_circle.json")

    parser = argparse.ArgumentParser(
        description="Plot local refinement GAP-SBM patch borders and element centers."
    )
    parser.add_argument(
        "project_parameters",
        nargs="?",
        default=str(default_project_parameters),
        help="Path to the ProjectParameters JSON file.",
    )
    parser.add_argument(
        "--output",
        dest="output",
        default=None,
        help="Optional output image path.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Save the figure without opening a matplotlib window.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    project_parameters_path = Path(args.project_parameters).resolve()
    if not project_parameters_path.is_file():
        raise FileNotFoundError(
            f"Project parameters file not found: {project_parameters_path}"
        )

    output_path = Path(args.output).resolve() if args.output else None
    build_and_plot(project_parameters_path, output_path, show=not args.no_show)


if __name__ == "__main__":
    main()
