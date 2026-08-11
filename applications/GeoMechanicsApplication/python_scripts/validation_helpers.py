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
    safe_base = Path(base_path).expanduser().resolve() if base_path else None
    input_dir = Path(input_path).expanduser()
    
    # Resolve relative inputs against the provided base_path
    if not input_dir.is_absolute() and safe_base is not None:
        input_dir = safe_base / input_dir
        
    input_dir = input_dir.resolve()
    
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
    candidate = Path(filename).expanduser()
    
    if not candidate.is_absolute():
        candidate = safe_base / candidate
        
    candidate = candidate.resolve(strict=False)
    _enforce_base_boundary(candidate, safe_base)

    if candidate.suffix.lower() != ".json":
        raise TypeError(f"Invalid file type '{candidate.suffix}'. Only .json files are accepted.")

    if not candidate.exists():
        raise FileNotFoundError(f"Parameter file not found: {candidate}")

    return candidate