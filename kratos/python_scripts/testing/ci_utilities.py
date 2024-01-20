import json
from argparse import ArgumentParser
from functools import lru_cache
from os import getenv
from pathlib import Path
from pprint import pprint
from typing import List, Set, Optional

from KratosMultiphysics.kratos_utilities import GetKratosMultiphysicsPath


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
    pprint(sorted(map(lambda p: p.as_posix(), changed_files())))
    pprint(sorted(ci_applications()))
    print(f"{sorted(get_changed_files_extensions())=}")
    print(f"{are_only_python_files_changed()=}")
    print(f"{sorted(get_changed_applications())=}")
    print(f"{is_core_changed()=}")
    print(f"{is_mpi_core_changed()=}\n")


def write_compiled_apps_to_file() -> None:
    """This function add the applications that are to be compiled to the environment
    For now this adds everything, but in the future this will be depending on the actual changes of the PR
    """
    with open("ci_compiled_apps.txt", "w") as ci_apps_file:
        kratos_path = Path(GetKratosMultiphysicsPath())
        for app in ci_applications():
            ci_apps_file.write(f"{kratos_path / 'applications' / app}\;")


def write_tested_apps_to_file() -> None:
    with open("ci_tested_apps.txt", "w") as ci_apps_file:
        kratos_path = Path(GetKratosMultiphysicsPath())
        for app in ci_applications():
            ci_apps_file.write(f"{kratos_path / 'applications' / app}\;")


def prepare_ci():
    print_ci_information()
    write_compiled_apps_to_file()
    write_tested_apps_to_file()

if __name__ == "__main__":
    prepare_ci()
    # parser: ArgumentParser = ArgumentParser(description="Utilities for the continuous integration")

    # option_groups = parser.add_mutually_exclusive_group(required=True)

    # option_groups.add_argument(
    #     "--compilation", action="store_true", help="Set the environment for the apps to be compiled"
    # )
    # option_groups.add_argument("--testing", action="store_true", help="Set the environment for the apps to be tested")
    # option_groups.add_argument("--info", action="store_true", help="Print information about this CI run")

    # args = parser.parse_args()

    # if args.compilation:
    #     write_compiled_apps_to_file()
    # elif args.testing:
    #     # not yet implemented
    #     pass
    # elif args.info:
    #     print_ci_information()
