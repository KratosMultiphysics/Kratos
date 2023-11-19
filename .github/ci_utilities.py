from pathlib import Path
import subprocess
from sys import argv, version
from os import getenv


def get_files_changed_in_pr(pr_number: int) -> list[Path] | None:
    try:
        process_output: str = subprocess.run(
            [
                "gh",
                "pr",
                "view",
                str(pr_number),
                "--json",
                "files",
                "--jq",
                ".files.[].path",
            ],
            check=True,
            capture_output=True,
        ).stdout.decode()

        modified_files: list[Path] = [Path(f) for f in process_output.splitlines()]

        return modified_files

    except:
        return None


if __name__ == "__main__":
    print(version)
    print(argv)
    print(f"{getenv('GITHUB_PR_NUMBER')=}")
    modified_files = get_files_changed_in_pr(argv[1])

    print(f"Modified files: {modified_files}")
