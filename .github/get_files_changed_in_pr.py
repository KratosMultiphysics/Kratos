import json
import subprocess
from os import environ


def get_files_changed_in_pr(pr_number: int) -> list[str]:
    # from: https://stackoverflow.com/questions/25071579/list-all-files-changed-in-a-pull-request-in-git-github

    modified_files: list[str] = []

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
            text=True,
        ).stdout

        modified_files = process_output.splitlines()

    except Exception as e:
        print(f"An error occured while getting the modified files: {e}")

    return modified_files


if __name__ == "__main__":
    assert "GITHUB_TOKEN" in environ, 'missing environment variable "GITHUB_TOKEN"!'
    assert "GITHUB_PR_NUMBER" in environ, 'missing environment variable "GITHUB_PR_NUMBER"!'

    pr_number: int = int(environ["GITHUB_PR_NUMBER"])
    modified_files = get_files_changed_in_pr(pr_number)

    print(json.dumps(modified_files))
