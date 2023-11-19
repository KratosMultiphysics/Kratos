from pathlib import Path
import subprocess


def get_files_changed_in_pr() -> list[Path] | None:
    try:
        subprocess.run(
            [
                "gh",
                "pr",
                "view",
                "${{ github.event.pull_request.number }}",
                "--json",
                "files",
                "--jq",
                "'.files.[].path'",
            ],
            check=True,
        )
    except:
        return None


if __name__ == "__main__":
    get_files_changed_in_pr()
