import os
import re
import toml
import build
import shutil
import fnmatch
import logging
import platform
import subprocess

from pathlib import Path

KRATOS_VERSION = "10.4.0"
PLATFORM_CONFIG = {
    'Linux': {
        'PYTHONS': ["38", "39", "310", "311", "312", "313", "314"],
        'PYTHON_BUILD_VER': "314",
        'BUILD_SCRIPT': "scripts/wheels/linux/configure.sh",
        'KRATOS_ROOT': "/workspace/kratos/Kratos",
        'WHEEL_ROOT': "/workspace/wheel",
        'WHEEL_OUT': "/data_swap_guest"
    },
    'Windows': {
        'PYTHONS': ["38", "39", "310", "311", "312", "313", "314"],
        'PYTHON_BUILD_VER': "39",
        'BUILD_SCRIPT': "scripts/wheels/windows/configure.bat",
        'BASE_LD_LIBRARY_PATH': "",
        'KRATOS_ROOT': "C:/kratos/Kratos",
        'WHEEL_ROOT': "C:/dist/wheel",
        'WHEEL_OUT': "C:/data_swap_guest/"
    },
    'Darwin': {
        'PYTHONS': ["39"],
        'PYTHON_BUILD_VER': "39",
        'BUILD_SCRIPT': "scripts/wheels/darwin/configure.sh",
        'KRATOS_ROOT': "/workspace/kratos/Kratos",
        'WHEEL_ROOT': "/workspace/wheel",
        'WHEEL_OUT': "/workspace/dist",
    }
}

def setupLogging(file):
    logging.basicConfig(filename=file, level=logging.DEBUG)
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger().addHandler(console)

def logResult(res):
    if res.stdout:
        for line in res.stdout.splitlines():
            logging.info(line)
    if res.stderr:
        for line in res.stderr.splitlines():
            logging.error(line)


def getAppList(kts_apps_dir):

    applications = []

    try:
        applications = [
            entry.name for entry in os.scandir(kts_apps_dir) if entry.is_dir() and "Application" in entry.name
        ]
    except FileNotFoundError as e:
        logging.critical(f"Error detecting apps in {kts_apps_dir}: {e}")

    return applications

def getPythonInterpreter(platform: str, python_ver: str):
    if platform == "Windows":
        return Path("C:/python") / f"{python_ver}" / "python.exe"
    elif platform == "Linux":
        return Path("/opt") / "python" / f"cp{python_ver}-cp{python_ver}" / "bin" / "python"
    elif platform == "Darwin":
        return Path("/opt") / "python" / f"cp{python_ver}-cp{python_ver}" / "bin" / "python"
    else:
        logging.critical(f"Cannot retrieve python interpreter for platform: {platform}")

def configure(CURRENT_CONFIG: dict, platform: str, python_ver: str):
    python_interpreter = getPythonInterpreter(platform, python_ver)

    if platform == "Windows":
        res = subprocess.run(
            [
                "powershell", 
                "-Command", 
                Path(CURRENT_CONFIG['KRATOS_ROOT']) / Path(CURRENT_CONFIG['BUILD_SCRIPT']),
                python_interpreter,
                Path(CURRENT_CONFIG['KRATOS_ROOT']),
                Path(CURRENT_CONFIG['KRATOS_ROOT']) / "bin" / "Release" / f"Python-{python_ver}",
                "24"
            ], 
            check=True
        )
        logResult(res)
    elif platform == "Linux":
        res = subprocess.run(
            [
                "bash", 
                Path(CURRENT_CONFIG['KRATOS_ROOT']) / "Scripts" / "Wheels" / platform.lower() / "configure.sh",
                python_interpreter,
                Path(CURRENT_CONFIG['KRATOS_ROOT']) / "bin" / "Release" / f"Python-{python_ver}" / "libs"
            ],
            check=True
        )
        logResult(res)
    elif platform == "Darwin":
        res = subprocess.run(
            [
                "bash", 
                Path(CURRENT_CONFIG['KRATOS_ROOT']) / "Scripts" / "Wheels" / platform.lower() / "configure.sh",
                python_interpreter,
                Path(CURRENT_CONFIG['KRATOS_ROOT']) / "bin" / "Release" / f"Python-{python_ver}" / "libs"
            ],
            check=True
        )
        logResult(res)
    else:
        logging.critical(f"Cannot build for platform: {platform}")
      
def buildKernel(CURRENT_CONFIG: dict, platform: str, python_ver: str):
    if platform == "Windows":
        # Currently windows insist on rebuilding all sources when we change python and not reuse the kernel. 
        pass
    else:
        configure(CURRENT_CONFIG, platform, python_ver)
        res = subprocess.run(['cmake', '--build', Path(CURRENT_CONFIG['KRATOS_ROOT']) / "build" / "Release", "--target", "KratosKernel"], check=True)
        logResult(res)

def buildInterface(CURRENT_CONFIG: dict, platform: str, python_ver: str):
    configure(CURRENT_CONFIG, platform, python_ver)

    if platform == "Windows":
        res = subprocess.run(['cmake', '--build', Path(CURRENT_CONFIG['KRATOS_ROOT']) / "build" / "Release", "--target", "install", "--", "/property:configuration=Release", "/p:Platform=x64"], check=True)
        logResult(res)
    else:
        res = subprocess.run(['cmake', '--build', Path(CURRENT_CONFIG['KRATOS_ROOT']) / "build" / "Release", "--target", "KratosPythonInterface"], check=True)
        logResult(res)
        res = subprocess.run(['cmake', '--build', Path(CURRENT_CONFIG['KRATOS_ROOT']) / "build" / "Release", "--target", "install"], check=True)
        logResult(res)

def buildWheel (CURRENT_CONFIG: dict, paths: dict, platform: str, python_ver: str):
    """ Builds the wheel for the current target.

    """
    
    wheel_staging_path = CURRENT_CONFIG['WHEEL_ROOT']
    
    setupWheelDir(wheel_staging_path)

    # Copy the python module contents
    shutil.copytree(
        src = paths['project_path'], 
        dst = Path(CURRENT_CONFIG['WHEEL_ROOT']) / paths['whl_dst_path'],
        dirs_exist_ok=True
    )

    # Copy the libraries
    shutil.copytree(
        src = paths['project_libs'], 
        dst = Path(CURRENT_CONFIG['WHEEL_ROOT']) / "libs",
        dirs_exist_ok=True
    )

    # Copy the correct pyproject file
    shutil.copy(
        src = paths['project_toml'],
        dst = Path(CURRENT_CONFIG['WHEEL_ROOT']) / os.path.basename(paths['project_toml'])
    )

    # Copy the build and metadata hooks
    shutil.copy(
        src = paths['project_hook'],
        dst = Path(CURRENT_CONFIG['WHEEL_ROOT']) / os.path.basename(paths['project_hook'])
    )

    logging.debug(Path(CURRENT_CONFIG['WHEEL_ROOT']))
    logging.debug(paths['project_read'])
    
    shutil.copy(
        src = paths['project_read'],
        dst = Path(CURRENT_CONFIG['WHEEL_ROOT']) / os.path.basename(paths['project_read'])
    )

    # Filter files if we are building separate wheels
    if paths['project_toml'] is not None:
        meta = toml.load(paths['project_toml'])
        libs = Path(CURRENT_CONFIG['WHEEL_ROOT']) / "libs"

        # Check if the file matches any of the patterns
        for item in libs.iterdir():
            if item.is_file() and not any(fnmatch.fnmatch(item.name, pattern) for pattern in meta['kratos']['libs']):
                logging.debug(f"Removing excluded file: {item.name}")
                os.remove(item)

    # Build
    python_interpreter = getPythonInterpreter(platform, python_ver)

    res = subprocess.run([python_interpreter, "-m", "build", "--wheel"], check=True, cwd=Path(CURRENT_CONFIG['WHEEL_ROOT']))
    logResult(res)

    cleanupWheelDir(wheel_staging_path)

def setupWheelDir(wheel_root_path: str):
    """
    Sets up the required directory structure for building wheels in the staging area.
    This creates the main wheel root and the KratosMultiphysics structure inside it.
    """
    
    logging.debug(f"Setting up Wheel Directory at: {wheel_root_path}")

    # Define the required directories to create inside the wheel root
    # This structure is relative to wheel_root_path
    target_structure = [
        '', 'KratosMultiphysics', os.path.join('KratosMultiphysics', '.libs')
    ]

    # Create all directories using the full path
    for relative_path in target_structure:
        dir_path = os.path.join(wheel_root_path, relative_path)
        try:
            # os.makedirs is Python's equivalent of 'mkdir -p'
            os.makedirs(dir_path, exist_ok=True)
            logging.debug(f"Created/Ensured directory: {dir_path}")
        except Exception as e:
            logging.critical(f"Error creating directory {dir_path}: {e}")
            return None # Return None if setup fails

    logging.debug("\nSetup complete.")
    return wheel_root_path

def cleanupWheelDir(wheel_root_path: str):
    """
    Removes the entire wheel building directory recursively.
    This is equivalent to 'rm -r $WHEEL_ROOT'.
    """

    logging.debug(f"Cleaning up Wheel Directory: {wheel_root_path}")
    
    if os.path.exists(wheel_root_path):
        try:
            # shutil.rmtree is Python's equivalent of 'rm -r'
            shutil.rmtree(wheel_root_path)
            logging.debug(f"Successfully removed directory: {wheel_root_path}")
        except Exception as e:
            logging.critical(f"Error removing directory {wheel_root_path}: {e}")
    else:
        logging.warning(f"Directory not found, nothing to clean: {wheel_root_path}")
    
    logging.info("Cleanup complete.\n")

if __name__ == "__main__":
    # Set some defines
    os.environ['KRATOS_VERSION'] = KRATOS_VERSION

    # TODO: Detect OS
    OS = platform.system()
    
    CURRENT_CONFIG = PLATFORM_CONFIG[OS]
    CURRENT_CONFIG['UNIFIED_WHEEL'] = False
    CURRENT_CONFIG['CLEAN_COMPILE'] = True

    # Setup Logging
    setupLogging(Path(CURRENT_CONFIG['WHEEL_OUT']) / "wheelbuild.log") 

    logging.info(f"Build mode is: {'Unified' if CURRENT_CONFIG['UNIFIED_WHEEL'] else 'Packages'}")

    logging.debug(f"Using KRATOS_ROOT: {CURRENT_CONFIG['KRATOS_ROOT']}")
    logging.debug(f"Targeting Python versions: {', '.join(CURRENT_CONFIG['PYTHONS'])}")

    # Complile the kernel if necessary.
    if CURRENT_CONFIG['CLEAN_COMPILE']:
        buildKernel(CURRENT_CONFIG, OS, CURRENT_CONFIG['PYTHON_BUILD_VER'])

    # Build for each python version
    for PYTHON in CURRENT_CONFIG['PYTHONS']:
        logging.info(f"Begining to build wheel for python version {PYTHON}")
        logging.info(f"Building Core...")

        main_whl_tree = Path(CURRENT_CONFIG['KRATOS_ROOT']) / "bin" / "Release" / f"Python-{PYTHON}" / "KratosMultiphysics"
        
        # Complile the interface.
        if CURRENT_CONFIG['CLEAN_COMPILE']:
            buildInterface(CURRENT_CONFIG, OS, PYTHON)

        # Set Project paths
        paths = {
            'project_path': main_whl_tree,
            'whl_dst_path': "KratosMultiphysics",
            'project_libs': Path(CURRENT_CONFIG['KRATOS_ROOT']) / "bin" / "Release" / f"Python-{PYTHON}" / "libs",
            'project_read': Path(CURRENT_CONFIG['KRATOS_ROOT']) / "README.md" ,
            'project_hook': Path(CURRENT_CONFIG['KRATOS_ROOT']) / "scripts" / "wheels" / "hatch_build.py"
        }

        if not CURRENT_CONFIG['UNIFIED_WHEEL']:
            paths['project_toml'] = Path(CURRENT_CONFIG['KRATOS_ROOT']) / "kratos" / "pyproject.toml"
        else:
            paths['project_toml'] = Path(CURRENT_CONFIG['KRATOS_ROOT']) / "scripts" / "wheels" / "pyproject.toml"

        # Create the core wheel
        buildWheel(CURRENT_CONFIG, paths, OS, PYTHON)
        
        # Build applications (if not unified build)
        if not CURRENT_CONFIG['UNIFIED_WHEEL']: 
            for APP in getAppList(Path(CURRENT_CONFIG['KRATOS_ROOT']) / "bin" / "Release" / f"Python-{PYTHON}" / "KratosMultiphysics"):
                logging.info(f"Building {APP}...")

                paths['project_path'] = main_whl_tree / APP
                paths['whl_dst_path'] = Path("KratosMultiphysics") / APP
                paths['project_toml'] = Path(CURRENT_CONFIG['KRATOS_ROOT']) / "applications" / APP / "pyproject.toml"
                paths['project_read'] = Path(CURRENT_CONFIG['KRATOS_ROOT']) / "applications" / APP / "README.md"

                buildWheel(CURRENT_CONFIG, paths)
            

