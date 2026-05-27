# VS Code Templates

Copy the files in this directory to your `.vscode/` folder inside the repository root.

## Setup

1. **Copy template files**
   ```bash
   mkdir -p .vscode
   cp scripts/vscode_templates/tasks.json .vscode/tasks.json
   ```

2. **Create your personalized build script**

   Copy the appropriate configure template and customize it for your machine (compilers, Python path, desired applications):
   ```bash
   # Linux
   cp scripts/standard_configure.sh build/configure.sh
   # Windows
   copy scripts\standard_configure.bat build\configure.bat
   ```

3. **Point the Build tasks to your script**

   In `.vscode/tasks.json`, update the `Build` and `MPI Build` commands to use your personalized script:
   ```json
   // Linux
   "command": "sh",
   "args": ["./build/configure.sh"]
   // Windows
   "command": ".\\build\\configure.bat"
   ```

4. **Windows: adjust the Visual Studio path** in the top-level `windows.options.shell.args` to match your VS installation (version and edition).

## Available Tasks

| Task | Group | Description |
|------|-------|-------------|
| `Build` | build (default) | Configure + compile |
| `MPI Build` | build | Configure + compile with `KRATOS_MPI_BUILD=ON` |
| `Run Tests` | test | Run all Python test suites |
| `Run Tests MPI` | test | Run Python MPI test suites |
| `Run CPP Tests` | test | Run all C++ GTest suites |
| `Run CPP Tests MPI` | test | Run C++ GTest suites under MPI |
| `Run CurrentFile` | test | Run the active Python file |
| `Run CurrentFile MPI` | test | Run the active Python file under MPI |
| `Run CurrentFile Valgrind` | test | Run the active file under Valgrind |
| `Run CurrentFile MPI Valgrind` | test | Run the active file under MPI + Valgrind |
| `Run CurrentFile VTune` | test | Profile the active file with Intel VTune |
| `Run CurrentFile cProfile` | test | Profile the active file with Python cProfile |
| `Run Current Benchmark file to JSON` | test | Run a Google Benchmark binary, output JSON |

Use `Shift+Ctrl+B` to trigger the default `Build` task.

## Inputs

Each task that needs a build type will prompt you to select one:

| Input | Options | Default |
|-------|---------|---------|
| `BuildType` | `Release`, `RelWithDebInfo`, `FullDebug` | `Release` |
| `CmakeGenerator` | `Unix Makefiles`, `Ninja` | `Unix Makefiles` |
| `NumberOfCores` | free text | `4` |
