"""Kratos's own VTK output as physicsnemo training data.

Every Kratos solve can already write `.vtu` files through the core
`VtkOutputProcess`, and physicsnemo 2.2 ships a `VTKReader` that consumes
one sample per subdirectory. Joining the two gives a training path that
needs no export process at all: point an existing simulation campaign's
output at a model.

What this path is NOT: it carries no PROVENANCE and offers no
scatter-back. A `.vtu` file is a rendering of the mesh, not the model part,
so nothing here can write a prediction onto Kratos entities - use the mesh
bridge for that. Treat it as a read-only path for training data that
already exists on disk.

WHICH READER. physicsnemo's `VTKReader` is NOT a general VTK reader: it
recognizes a fixed vocabulary of external-aerodynamics keys
(`stl_coordinates`, `surface_normals`, `volume_mesh_centers`,
`volume_fields`, ...) and returns an EMPTY sample for anything else - a
`.vtu` carrying a field called `PRESSURE` reads as nothing at all, with no
error. It is the right reader for a DrivAer-style dataset and the wrong one
for ordinary solver output.

The general path is `physicsnemo.mesh.io.from_pyvista`, which reads any
file pyvista can and auto-triangulates polyhedra, so hexahedral and mixed
meshes come through as simplices without a tessellation step.
`CreateVtkMeshDataset` is that path as a torch Dataset.

pyvista and physicsnemo are imported lazily.
"""

import shutil
from pathlib import Path

import KratosMultiphysics as Kratos


def _TryImportPyVista():
    try:
        import pyvista
        return pyvista
    except ImportError as e:
        raise ImportError(
            "Reading VTK output requires pyvista, which could not be imported. "
            "Install it with e.g. 'pip install pyvista'.") from e


def _TryImportVtkReader():
    try:
        from physicsnemo.datapipes.readers import VTKReader
        return VTKReader
    except ImportError as e:
        raise ImportError(
            "VTKReader requires physicsnemo >= 2.2, which could not be imported. "
            "Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def _TryImportFromPyVista():
    try:
        from physicsnemo.mesh.io import from_pyvista
        return from_pyvista
    except ImportError as e:
        raise ImportError(
            "Converting a pyvista mesh requires physicsnemo >= 2.2, which could not be "
            "imported. Install it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


_READER_EXTENSIONS = (".vtu", ".vtp", ".stl")


def ArrangeVtkOutputForReader(vtk_output_directory, destination, pattern: str = "*.vtk",
                              sample_prefix: str = "sample", move: bool = False) -> int:
    """Lays a VtkOutputProcess directory out the way VTKReader expects it.

    Two mismatches to bridge, not one:

    - The reader treats each SUBDIRECTORY as one sample, while
      VtkOutputProcess writes one flat file per step per model part
      (``<Part>_<rank>_<step>.vtk``).
    - Kratos writes the LEGACY ``.vtk`` format, and physicsnemo's VTKReader
      reads only ``.vtu``, ``.vtp`` and ``.stl``. Pointing the reader at a
      Kratos output directory therefore finds nothing at all, which is why
      this function converts through pyvista rather than just copying.

    Args:
        vtk_output_directory: Where the solver wrote its files.
        destination: The directory to build (created if missing).
        pattern: Which files to take. The default takes Kratos's own
            ``.vtk``; narrow it (e.g. ``"Main_0_*.vtk"``) to pick one model
            part out of several.
        sample_prefix: Subdirectory name prefix.
        move: Delete the source file after converting it.

    Returns:
        How many samples were arranged.
    """
    source = Path(str(vtk_output_directory))
    destination = Path(str(destination))
    files = sorted(source.glob(pattern))
    if not files:
        raise ValueError(
            f"No files matching \"{pattern}\" in \"{source}\"; VtkOutputProcess writes "
            "one file per step per model part, named <Part>_<rank>_<step>.vtk.")
    destination.mkdir(parents=True, exist_ok=True)

    needs_conversion = any(path.suffix not in _READER_EXTENSIONS for path in files)
    pyvista = _TryImportPyVista() if needs_conversion else None

    for index, path in enumerate(files):
        sample_directory = destination / f"{sample_prefix}_{index:04d}"
        sample_directory.mkdir(exist_ok=True)
        if path.suffix in _READER_EXTENSIONS:
            shutil.copy2(str(path), str(sample_directory / path.name))
        else:
            pyvista.read(str(path)).save(str(sample_directory / (path.stem + ".vtu")))
        if move:
            path.unlink()
    return len(files)


def CreateVtkReaderDataset(directory, keys_to_read=None, **options):
    """A physicsnemo VTKReader over an ArrangeVtkOutputForReader layout.

    Args:
        directory: The arranged directory.
        keys_to_read: Which arrays to read; None reads what it finds.
        options: Forwarded upstream (exclude_patterns, pin_memory, ...).

    Returns:
        A VTKReader, indexable as reader[i] -> (TensorDict, metadata).
    """
    VTKReader = _TryImportVtkReader()
    return VTKReader(str(directory), keys_to_read=keys_to_read, **options)


def MeshFromVtkFile(path, manifold_dim="auto", **options):
    """One `.vtu`/`.vtp`/`.stl` file as a physicsnemo Mesh.

    Polyhedral cells are auto-triangulated on the way in, so a hexahedral
    Kratos mesh arrives as simplices with no tessellation step here.
    """
    pyvista = _TryImportPyVista()
    from_pyvista = _TryImportFromPyVista()
    return from_pyvista(pyvista.read(str(path)), manifold_dim=manifold_dim, **options)


def CreateVtkMeshDataset(directory, pattern: str = "*", field_names=None,
                         manifold_dim="auto"):
    """A torch Dataset of physicsnemo Meshes over a directory of VTK files.

    The general counterpart of CreateVtkReaderDataset: it reads whatever
    fields the files carry, by their own names, instead of the upstream
    reader's fixed external-aero vocabulary. Files are taken in sorted
    order, both flat and one-per-subdirectory layouts.

    Args:
        directory: Where the files are (searched recursively).
        pattern: Glob for the file stems, e.g. "Main_0_*".
        field_names: Restrict to these point-data arrays; None keeps all.
        manifold_dim: Forwarded to from_pyvista.

    Returns:
        A torch Dataset yielding physicsnemo Meshes.
    """
    torch = _TryImportTorch()
    _TryImportPyVista()
    _TryImportFromPyVista()

    directory = Path(str(directory))
    files = sorted(
        path for extension in (".vtu", ".vtp", ".stl", ".vtk")
        for path in directory.rglob(pattern + extension))
    if not files:
        raise ValueError(
            f"No VTK files matching \"{pattern}\" under \"{directory}\".")

    class _VtkMeshDataset(torch.utils.data.Dataset):
        def __init__(self, paths):
            self.paths = paths

        def __len__(self):
            return len(self.paths)

        def __getitem__(self, index):
            mesh = MeshFromVtkFile(self.paths[index], manifold_dim=manifold_dim)
            if field_names is not None:
                for key in list(mesh.point_data.keys()):
                    if key not in field_names:
                        del mesh.point_data[key]
            return mesh

    return _VtkMeshDataset(files)


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.vtk_bridge's dataset requires torch, which could "
            "not be imported. Install it with e.g. 'pip install torch'.") from e
