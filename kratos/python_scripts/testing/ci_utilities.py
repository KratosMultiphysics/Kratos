import json
from functools import lru_cache
from os import getenv
from pathlib import Path
from pprint import pprint
from typing import List, Set


def check_valid_environment_configuration_exists() -> None:
    if not getenv("KRATOS_CI_CHANGED_FILES"):
        raise RuntimeError("Invalid CI-environment: KRATOS_CI_CHANGED_FILES")

    kratos_ci_applications: Optional[str] = getenv("KRATOS_CI_APPLICATIONS")

    if kratos_ci_applications == "ONLY_CORE":
        return

    if kratos_ci_applications is None:
        raise RuntimeError("Invalid CI-environment: KRATOS_CI_APPLICATIONS")

    if not Path(kratos_ci_applications).exists():
        raise RuntimeError("Invalid CI-environment: KRATOS_CI_APPLICATIONS file does not exist")


@lru_cache
def changed_files() -> List[Path]:
    check_valid_environment_configuration_exists()
    return [Path(f) for f in json.loads(getenv("KRATOS_CI_CHANGED_FILES"))]


def ci_applications() -> List[str]:
    check_valid_environment_configuration_exists()

    kratos_ci_applications: str = getenv("KRATOS_CI_APPLICATIONS")

    if kratos_ci_applications == "ONLY_CORE":
        return []

    with open(kratos_ci_applications) as ci_apps_file:
        return json.load(ci_apps_file)


def get_changed_applications() -> Set[str]:
    return {f.parts[1] for f in changed_files() if f.parts[0] == "applications"}


def is_core_changed() -> bool:
    return any([f.parts[0] == "kratos" for f in changed_files()])


def is_mpi_core_changed() -> bool:
    return any([f.parts[0] == "kratos" and f.parts[1] == "mpi" for f in changed_files()])


def get_changed_files_extensions() -> Set[str]:
    return {f.suffix for f in changed_files()}


def are_only_python_files_changed() -> bool:
    return get_changed_files_extensions() == {".py"}


def print_ci_information() -> None:
    """This function prints an overview of the CI related information"""
    pprint(sorted(map(lambda p : p.as_posix(), changed_files())))
    pprint(sorted(ci_applications()))
    print(f"{sorted(get_changed_files_extensions())=}")
    print(f"{are_only_python_files_changed()=}")
    print(f"{sorted(get_changed_applications())=}")
    print(f"{is_core_changed()=}")
    print(f"{is_mpi_core_changed()=}\n")


if __name__ == "__main__":
    print_ci_information()
