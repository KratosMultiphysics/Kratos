# python imports
from shutil import which

def GetPython3Command():
    """Return the name of the python command, can be used with subprocess."""
    if which("python3"):
        return "python3"
    elif which("python"):
        return "python"
    raise Exception("The python command could not be determined!")
