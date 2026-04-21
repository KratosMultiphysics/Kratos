# Kratos Wheel Generation: Understanding PyProject files

## TLDR

Here is a quick summary of how to generate Kratos wheels depending on your deployment target:

### 1. Generating a local release from the build dir
After building Kratos with CMake, simply navigate to your newly generated release folder (which will contain the `libs` sub-folder with your `.so`/`.dll` files) and run the standard python build module from there:
```bash
cd bin/Release/
python3 -m build --wheel
```

### 2. Generating a bundled release (Unified Wheel)
To build a single "fat" wheel containing the Kratos Core alongside all compiled Applications, use the provided package orchestration script. Inside `scripts/wheels/build_wheel.py`, ensure the flag `CURRENT_CONFIG['UNIFIED_WHEEL'] = True`, then execute it from the project root: 
```bash
python3 scripts/wheels/build_wheel.py
```

### 3. Generating a regular release (Different wheels for Core and Apps)
To build an ecosystem of separate wheels (one main wheel for the Core, and one standalone wheel per Application), ensure exactly the opposite: `CURRENT_CONFIG['UNIFIED_WHEEL'] = False` (which is typically the default). Run the orchestrator:
```bash
python3 scripts/wheels/build_wheel.py
```

## Introduction

In Kratos Multiphysics, wheel packaging is governed by several `pyproject.toml` files using standard PEP-517 tools with `hatchling` as the build system. This document outlines the purpose and critical features of the `pyproject.toml` files found in the core and scripts directories. 

## Context: Target Execution Path

> [!WARNING]
> While these files are physically authored at `kratos/pyproject.toml` and `scripts/wheels/pyproject.toml` (or similar), **they are not executed from these source locations**. 

During the CMake configuration and build step, the required `pyproject.toml` files are orchestrated and copied to build-staging directories or final release binary directories. Wheels must be built from these target staging locations, as that is where all the compiled dynamic libraries (`.so` / `.dll`) and Python wrappers securely reside at build time.

## 1. Core vs Unified Wheel Generators

Kratos defines distinct strategies to package either the minimal core or a "fat" unified package containing all applications:

### The Core `pyproject.toml` (`./kratos`)
**Purpose:** This configuration generates a wheel corresponding exclusively to the Kratos kernel. It bundles only the core files, the MPI core wrappers, and fundamental system orchestrators.
**Location in Source:** `kratos/pyproject.toml`

### The Applications `pyproject.toml` (`./applications/[APP_NAME]/pyproject.toml`)
**Purpose:** This configuration generates a wheel corresponding to the selected application. It is ment to be used along with the core pyproject.toml file to generate application wheels that depend on the core wheel.
**Location in Source:** `kratos/pyproject.toml`

### The Unified `pyproject.toml` (`./scripts`)
**Purpose:** This project builds Kratos and **all compiled applications** from your local environment as a single bundled package. 
**Location in Source:** `scripts/wheels/pyproject.toml`

## 2. `pyproject.toml` Structure

Both configs share a common anatomy dictated by the `hatchling` backend, but handle certain specifics differently.

### Standard Metadata (`[build-system]` and `[project]`)
Defines Kratos Multiphysics as the package name, the backend used (`hatch.build`), target python versions, and dependencies:
```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "KratosMultiphysics"
dependencies = ["numpy>=1.20"]
dynamic = ["version"]
```

### The Kratos Specific Section `[kratos]` (Core only)
Because `kratos/pyproject.toml` handles the core natively alongside application builds, it contains a `[kratos]` table explicitly defining the C++ libraries needed for the kernel build.

```toml
[kratos]
libs = [
    "Kratos.*",
    "KratosMPI.*",
    "KratosCore.*",
    # ...
]
```
> [!NOTE]
> This section is parsed by the Kratos wrapping orchestration script `build_wheel.py` to filter precisely which binaries belong to the core package versus applications and is not part of a standard pyproject.toml file. The unified pyproject omits this because it intentionally scoops up *everything* in the binary tree.

### Build Hooks `[tool.hatch.build.targets.wheel.hooks.custom]`
Kratos wheels are strictly platform-specific and are not purely python code. These hooks execute a local `hatch_build.py` script which forces the final wheel to be recognized as system-specific. 
The hook script typically sets:
- `build_data["infer_tag"] = True`
- `build_data["pure_pyton"] = False`

This ensures the wheel filenames correspond exactly to the current OS architecture (e.g., `linux_x86_64` or `win_amd64`) instead of `any`. It also injects the `KRATOS_VERSION` environment variable.

## 3. `.py`, `.dll`, and `.so` Files

### Python source and Hinting files inclusion
Python files native to the framework are injected explicitly via standard `include` directives.
```toml
[tool.hatch.build.targets.wheel]
include = [
    "KratosMultiphysics/**/*.py",
    "KratosMultiphysics/**/*.pyi"
]
```

### Including Compiled C++ Artifacts (`.dll` / `.so`)
Packaging `.dll` (Windows), `.so` (Linux) and `.dylib` (Mac) files side-by-side with Python code is structurally the most critical operation of Kratos wheels. Kratos relies on `hatchling`'s `force-include` mapping technique mapping:

```toml
[tool.hatch.build.targets.wheel.force-include]
libs = "KratosMultiphysics/.libs"
```

**How it works:**
1. A staging folder named `libs` is set up containing all required `.so` or `.dll` binaries.
2. The key string (`libs`) in the TOML corresponds to this physical layout.
3. The assigned value (`"KratosMultiphysics/.libs"`) forces `hatchling` to place the files into that destination directory structure inside the final `.whl` archive.
4. Hence, all collected compiled `.so`/`.dll`/`.dylib` packages are pushed into the hidden `.libs` directory, correctly satisfying dynamic links at runtime.

---

## 4. The `build_wheel.py` Script

While it is possible to build wheels manually by fulfilling the `force-include` prerequisites, Kratos provides a robust automation script to handle this: `scripts/wheels/build_wheel.py`.

### Purpose
The `build_wheel.py` script is the main driver for generating Kratos Python Wheels. It orchestrates the CMake build of Kratos Python interfaces and coordinates the PEP-517 wheel packaging across multiple operating systems and Python versions.

### Execution Flow & Usage

**1. Managing the Build Target (Unified vs Separate Packages)**
The script contains a configuration rule `CURRENT_CONFIG['UNIFIED_WHEEL']`:
- **If `False` (Packages Target):** The script sets up the build for the "Core" using `kratos/pyproject.toml`. It then dynamically detects all compiled directories inside `applications/*Application*` and generates a standalone `.whl` extension package for *each application*.
- **If `True` (Unified Target):** The script executes a single package generation utilizing `scripts/wheels/pyproject.toml`. This creates one massive package containing both the core and all application libraries.

**2. Staging Environment (`setupWheelDir`)**
It dynamically generates an isolated `WHEEL_ROOT` staging directory. For each iteration, the script copies:
- The compiled Python wrappers (`.py` files).
- The raw binary libraries (`.dll` / `.so` mapped to a folder called `libs`).
- The relevant `pyproject.toml`.
- The `README.md` and `hatch_build.py` hook.

**3. Binary Filtering Using `[kratos]` Config**
Because applications and core share a master binary release folder (`bin/Release/libs`), the script performs intelligent filtering. 
It uses Python's `toml.load` to read the target `pyproject.toml`. If the file contains a `[kratos] libs = [...]` list (like `kratos/pyproject.toml`), the orchestration script aggressively iterates over the physics binaries within the staging `libs` directory and **deletes any binaries that do not match the declared patterns**. 
This strictly ensures that the Core wheel gets only Core binaries, and Application wheels get only their respective subsets.

**4. Executing Hatchling & Wheel Output**
After staging and filtering, the orchestrator triggers the Python build module natively inside the temporary root:
```python
subprocess.run([python_interpreter, "-m", "build", "--wheel", ...])
```
Resulting `.whl` files are exported to the target `--outdir` (often `WHEEL_OUT`), and the staging configuration gets safely destroyed via `cleanupWheelDir`.
