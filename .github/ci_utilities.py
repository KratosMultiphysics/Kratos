from pathlib import Path
import subprocess
from sys import argv


def get_files_changed_in_pr(pr_number: int) -> list[Path] | None:
    try:
        # subprocess.run(
        #     ["gh pr view ${{ github.event.pull_request.number }} --json files --jq '.files.[].path'"],
        #     check=True,
        # )
        p = subprocess.run(
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
            capture_output=True
        )

        print(f"Process output: {p.stdout.decode('ascii')}")
    except:
        return None


if __name__ == "__main__":
    print(argv)
    get_files_changed_in_pr(argv[1])
