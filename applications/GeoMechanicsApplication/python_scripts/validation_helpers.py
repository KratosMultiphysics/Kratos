from pathlib import Path
from typing import List

def _enforce_base_boundary(resolved_path: Path, base_dir: Path) -> None:
    """Ensures resolved_path stays within base_dir or raises TypeError."""
    try:
        resolved_path.relative_to(base_dir)
    except ValueError as exc:
        raise TypeError(
            f"Path escapes allowed directory '{base_dir}': {resolved_path}"
        ) from exc

def _resolve_safe_path(
    user_input: str | Path, 
    base_dir: Path | None = None,
    *,
    strict: bool = False
) -> Path:
    """
    Expands user shortcuts and resolves paths safely against an optional base.
    """
    candidate = Path(user_input).expanduser()
    if not candidate.is_absolute() and base_dir is not None:
        candidate = base_dir / candidate
    return candidate.resolve(strict=strict)

def validated_stage_file_paths(
    input_path: str | Path,
    n_stages: int,
    filename_pattern: str,
    base_path: str | Path | None = None
) -> List[Path]:
    """
    Build stage file paths and ensure they stay within the configured input directory.
    Accepts both string and Path objects for flexibility.
    """
    safe_base = _resolve_safe_path(base_path) if base_path else None
    input_dir = _resolve_safe_path(input_path, base_dir=safe_base)
        
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input path does not exist or is not a directory: {input_dir}")

    result = []
    for i in range(n_stages):
        candidate = (input_dir / filename_pattern.format(i + 1)).resolve()
        _enforce_base_boundary(candidate, input_dir)
        result.append(candidate)

    return result

def validated_parameter_path(filename: str | Path) -> Path:
    """
    Validates that 'filename' resolves within the current working directory 
    and has a .json extension.
    """
    safe_base = Path.cwd().resolve()
    candidate = _resolve_safe_path(filename, base_dir=safe_base, strict=False)

    _enforce_base_boundary(candidate, safe_base)

    if candidate.suffix.lower() != ".json":
        raise TypeError(f"Invalid file type '{candidate.suffix}'. Only .json files are accepted.")

    if not candidate.exists():
        raise FileNotFoundError(f"Parameter file not found: {candidate}")

    return candidate