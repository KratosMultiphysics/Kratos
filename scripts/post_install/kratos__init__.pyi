from dataclasses import dataclass
import Kratos.python_registry as python_registry
import Kratos.kratos_globals as kratos_globals

@dataclass
class KratosPaths:
    kratos_install_path: str
    kratos_libs: str
    kratos_module_libs: str
    kratos_applications: str
    kratos_scripts: str
    kratos_tests: str


KratosGlobals: kratos_globals.KratosGlobalsImpl
Registry: python_registry.PythonRegistry
RegisterPrototype = python_registry.RegisterPrototype
python_version: str
kratos_version_info: str

def IsDistributedRun() -> bool: ...
