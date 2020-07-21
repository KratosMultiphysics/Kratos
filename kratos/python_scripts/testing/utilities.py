# python imports
import sys

def GetPython3Command():
    """Return the name of the python command, can be used with subprocess."""
    sys_executable = sys.executable
    if sys_executable: # was found
        return sys_executable
    raise Exception("The python command could not be determined!")
