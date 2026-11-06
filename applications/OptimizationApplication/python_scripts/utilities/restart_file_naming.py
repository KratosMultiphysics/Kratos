import re
from pathlib import Path

import KratosMultiphysics as Kratos

def SplitRestartFileName(restart_file_name: str) -> 'tuple[Path, str]':
    """Splits "dir/restart_<step>.pkl" into (Path("dir"), "restart_<step>.pkl"); no dir gives Path(".")."""
    path = Path(restart_file_name)
    return path.parent, path.name

def ParseStepFromFileName(restart_file_name: str, file_name: str) -> 'int | None':
    """Inverse of restart_file_name.replace("<step>", str(step)).

    Returns None unless the whole file_name is restart_file_name with "<step>" replaced by ASCII
    digits, e.g. "restart_5_old.pkl" and the temp file "restart_5.pkl.tmp" give None for "restart_<step>.pkl".
    """
    prefix, suffix = restart_file_name.split("<step>")
    match = re.fullmatch(re.escape(prefix) + "([0-9]+)" + re.escape(suffix), file_name)
    return int(match.group(1)) if match else None

def ListCheckpoints(restart_files_path: Path, restart_file_name: str) -> 'list[tuple[int, Path]]':
    """Returns (step, path) of every checkpoint in restart_files_path, sorted by step.

    Only regular files whose full name parses to a step count (no temp files or directories).
    """
    if not restart_files_path.is_dir():
        return []
    checkpoints = []
    for path in restart_files_path.glob(restart_file_name.replace("<step>", "*")):
        step = ParseStepFromFileName(restart_file_name, path.name)
        if step is not None and path.is_file():
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
