from pathlib import Path

import KratosMultiphysics as Kratos

def SplitRestartFileName(restart_file_name: str) -> 'tuple[Path, str]':
    """Splits "dir/restart_<step>.pkl" into (Path("dir"), "restart_<step>.pkl"); no dir gives Path(".")."""
    path = Path(restart_file_name)
    return path.parent, path.name

def ParseStepFromFileName(restart_file_name: str, file_name: str) -> 'int | None':
    """Inverse of restart_file_name.replace("<step>", str(step)).

    Returns None if file_name does not match restart_file_name with a non-negative integer step
    (e.g. "restart_5_old.pkl" for "restart_<step>.pkl").
    """
    prefix, suffix = restart_file_name.split("<step>")
    if not file_name.startswith(prefix) or not file_name.endswith(suffix) or len(file_name) <= len(prefix) + len(suffix):
        return None
    step_str = file_name[len(prefix):len(file_name) - len(suffix)]
    return int(step_str) if step_str.isdigit() else None

def ListCheckpoints(restart_files_path: Path, restart_file_name: str) -> 'list[tuple[int, Path]]':
    """Returns (step, path) of every checkpoint in restart_files_path, sorted by step.

    Files that match the glob but whose step cannot be parsed are ignored.
    """
    if not restart_files_path.is_dir():
        return []
    checkpoints = []
    for path in restart_files_path.glob(restart_file_name.replace("<step>", "*")):
        step = ParseStepFromFileName(restart_file_name, path.name)
        if step is not None:
            checkpoints.append((step, path))
    return sorted(checkpoints, key=lambda checkpoint: checkpoint[0])

def ResolveRestartLoadStep(restart_files_path: Path, restart_file_name: str, restart_load_step_param: Kratos.Parameters) -> 'int | None':
    """Resolves "restart_load_step" to a step.

    "latest" is the highest parsed step (not the newest mtime, which breaks on copied checkpoints)
    and gives None when no checkpoint exists yet; an int is returned as is. Any other string raises.
    """
    if restart_load_step_param.IsString():
        value = restart_load_step_param.GetString()
        if value != "latest":
            raise RuntimeError(f"\"restart_load_step\" must be \"latest\" or an int [ restart_load_step = \"{value}\" ].")
        checkpoints = ListCheckpoints(restart_files_path, restart_file_name)
        return checkpoints[-1][0] if checkpoints else None
    return restart_load_step_param.GetInt()
