from pathlib import Path

import KratosMultiphysics as Kratos

def SplitRestartFileName(restart_file_name: str) -> 'tuple[Path, str]':
    """Splits a "restart_file_name" setting into its directory and bare filename template, e.g.
    "Optimization_Restart/restart_<step>.pkl" -> (Path("Optimization_Restart"), "restart_<step>.pkl").
    A bare template with no directory (e.g. "restart_<step>.pkl") splits to (Path("."), ...)."""
    path = Path(restart_file_name)
    return path.parent, path.name

def GetModelPartRestartFileBaseName(restart_file_name: str, step: int, model_part_name: str) -> str:
    """Derives the base file name (no ".rest" -- Kratos.FileSerializer appends it itself) for a
    given model part's restart file at the given step, from the checkpoint payload's own
    "<step>"-templated restart_file_name, e.g. "restart_<step>.pkl" -> "restart_5_Structure"."""
    stem = Path(restart_file_name.replace("<step>", str(step))).stem
    return f"{stem}_{model_part_name}"

def ParseStepFromFileName(restart_file_name: str, file_name: str) -> int:
    """Inverse of restart_file_name.replace("<step>", str(step))."""
    prefix, suffix = restart_file_name.split("<step>")
    end = len(file_name) - len(suffix) if suffix else len(file_name)
    return int(file_name[len(prefix):end])

def ResolveRestartLoadStep(restart_files_path: Path, restart_file_name: str, restart_load_step_param: Kratos.Parameters) -> 'int | None':
    """Resolves "restart_load_step" ("latest" or an explicit int) to a concrete step number.

    Returns None only for "latest" when no checkpoint exists yet (fresh run).
    """
    if restart_load_step_param.IsString():
        checkpoint_glob = restart_file_name.replace("<step>", "*")
        checkpoints = list(restart_files_path.glob(checkpoint_glob)) if restart_files_path.is_dir() else []
        if not checkpoints:
            return None
        latest = max(checkpoints, key=lambda p: p.stat().st_mtime)
        return ParseStepFromFileName(restart_file_name, latest.name)
    return restart_load_step_param.GetInt()
