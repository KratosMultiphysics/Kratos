import json
from argparse import ArgumentParser
from functools import lru_cache
from os import getenv
from pathlib import Path
from typing import List, Set

CHANGED_FILES_CI = Path("ci_changed_files.json")


def write_changed_files_to_disk(changed_files: List[str]) -> None:
    with open(CHANGED_FILES_CI, "w", encoding="utf-8") as ci_file:
        json.dump(changed_files, ci_file, indent=4)


@lru_cache
def changed_files() -> List[Path]:
    with open(CHANGED_FILES_CI, encoding="utf-8") as ci_file:
        changed_files = json.load(ci_file)

    return [Path(f) for f in changed_files]


def get_changed_applications() -> List[str]:
    return list({f.parts[1] for f in changed_files() if f.parts[0] == "applications"})


def is_core_changed() -> bool:
    return any([f.parts[0] == "kratos" for f in changed_files()])


def is_mpi_core_changed() -> bool:
    return any([f.parts[0] == "kratos" and f.parts[1] == "mpi" for f in changed_files()])


def get_changed_files_extensions() -> Set[str]:
    return {f.suffix for f in changed_files()}


def are_only_python_files_changed() -> bool:
    return get_changed_files_extensions() == {".py"}


if __name__ == "__main__":
    parser: ArgumentParser = ArgumentParser()

    parser.add_argument("-w", "--write", help="Write changed files as json string format", nargs="*")

    args = parser.parse_args()

    if args.write:
        write_changed_files_to_disk(args.write)

    print(f"{get_changed_files_extensions()=}")
    print(f"{are_only_python_files_changed()=}")
    print(f"{get_changed_applications()=}")
    print(f"{is_core_changed()=}")
    print(f"{is_mpi_core_changed()=}")
    print(f"{getenv('CI_CHANGED_FILES')=}")
