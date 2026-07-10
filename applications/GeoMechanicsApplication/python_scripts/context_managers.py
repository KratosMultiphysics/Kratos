import contextlib
import os
import pathlib


# Based on the following post: https://dev.to/teckert/changing-directory-with-a-python-context-manager-2bj8
@contextlib.contextmanager
def set_cwd_to(path: pathlib.Path):
    """Sets the current working directory within the present context"""
    cwd_on_entry = pathlib.Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd_on_entry)
