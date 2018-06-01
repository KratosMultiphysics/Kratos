from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import ctypes
import io
import os, sys
import tempfile
import h5py
from contextlib import contextmanager

@contextmanager
def stdout_redirector( stream, c_stdout, libc ):
    # The original fd stdout points to. Usually 1 on POSIX systems.
    original_stdout_fd = sys.stdout.fileno()

    def _redirect_stdout(to_fd):
        """Redirect stdout to the given file descriptor."""
        # Flush the C-level buffer stdout
        libc.fflush(c_stdout)
        # Flush and close sys.stdout - also closes the file descriptor (fd)
        sys.stdout.close()
        # Make original_stdout_fd point to the same file as to_fd
        os.dup2(to_fd, original_stdout_fd)
        # Create a new sys.stdout that points to the redirected fd
        sys.stdout = io.TextIOWrapper(os.fdopen(original_stdout_fd, 'wb'))

    # Save a copy of the original stdout fd in saved_stdout_fd
    saved_stdout_fd = os.dup(original_stdout_fd)
    try:
        # Create a temporary file and redirect stdout to it
        tfile = tempfile.TemporaryFile(mode='w+b')
        _redirect_stdout(tfile.fileno())
        # Yield to caller, then redirect stdout back to the saved fd
        yield
        _redirect_stdout(saved_stdout_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(str(tfile.read()).replace('\\n', '\n' ))
    finally:
        tfile.close()
        os.close(saved_stdout_fd)

def CreateResponseFunction(response_id, response_settings, model_part):
    return AnalysisDriverBasedResponseFunction(response_id, response_settings, model_part)

class ResponseFunctionBase(object):
    """The base class for structural response functions. Each response function
    is able to calculate its response value and gradient.
    All the necessary steps have to be implemented, like e.g. initializing,
    solving of primal (and adjoint) analysis ...
    """

    def RunCalculation(self, calculate_gradient):
        self.Initialize()
        self.InitializeSolutionStep()
        self.CalculateValue()
        if calculate_gradient:
            self.CalculateGradient()
        self.FinalizeSolutionStep()
        self.Finalize()

    def Initialize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def CalculateValue(self):
        raise NotImplementedError("CalculateValue needs to be implemented by the derived class")

    def CalculateGradient(self):
        raise NotImplementedError("CalculateGradient needs to be implemented by the derived class")

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        pass

    def GetValue(self):
        raise NotImplementedError("GetValue needs to be implemented by the derived class")

    def GetShapeGradient(self):
        raise NotImplementedError("GetShapeGradient needs to be implemented by the derived class")

class AnalysisDriverBasedResponseFunction(ResponseFunctionBase):

    def __init__(self, response_id, response_settings, model_part):
        self.model_part = model_part
        self.response_settings = response_settings
        self.response_id = response_id

        self.model_part_filename = response_settings["optimization_model_part_name"].GetString()
        self.analysis_driver_name = response_settings["analysis_driver"].GetString()
        self.log_file = "%s.log" % self.analysis_driver_name
        self.results_file = "%s.results.h5" % self.analysis_driver_name

        self.analysis_driver = __import__(self.analysis_driver_name)

        self.is_analysis_step_completed = False


    def InitializeSolutionStep(self):
        self.is_analysis_step_completed = self.__IsAnalysisCompleted()
        
        if not self.is_analysis_step_completed:
            model_part_io = ModelPartIO(self.model_part_filename, ModelPartIO.WRITE)
            model_part_io.WriteModelPart(self.model_part)
        else:
            # read results from the file.
            print("> Iteration is already completed. Reading results...")
            self.__ReadResults()

    def CalculateValue(self):
        if not self.is_analysis_step_completed:
            self.response_data = {}
            libc = ctypes.CDLL(None)
            c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')
            std_io_out = io.StringIO()           
            with stdout_redirector( std_io_out, c_stdout, libc ):
                self.analysis_driver.Run(self.model_part, self.response_data)
            std_file_out = open(self.log_file , "w")
            std_file_out.write(std_io_out.getvalue())
            std_file_out.close()

    def GetValue(self):
        return self.response_data["response_value"]

    def GetShapeGradient(self):
        return self.response_data["response_gradients"]           

    def FinalizeSolutionStep(self):
        if not self.is_analysis_step_completed:
            # write the final results to have restart capabilities
            self.__WriteResults()
            log_file = open(self.log_file, "a")
            log_file.write("\n\nKRATOS RESPONSE VALUE SUCCESSFULLY CALCULATED.\n")
            log_file.close()

    def __ReadResults(self):
        results_file = h5py.File(self.results_file, "r")
        self.response_data = {}
        self.response_data["response_value"] = results_file["ResponseValue"]
        _temp_results = results_file["ResponseGradients"]
        print(_temp_results)
        self.response_data["response_gradients"] = {}
        for i in range(0, len(_temp_results)):
            _node = int(_temp_results[i][0])
            _gradient = Vector(3)
            _gradient[0] = _temp_results[i][1]
            _gradient[1] = _temp_results[i][2]
            _gradient[2] = _temp_results[i][3]
            self.response_data["response_gradients"][_node] = _gradient
        results_file.close()

    def __WriteResults(self):
        results_file = h5py.File(self.results_file, "w")
        results_file["ResponseValue"] = self.GetValue()
        _response_gradients = []
        for _node in self.GetShapeGradient():
            _response_gradients.append([_node, self.GetShapeGradient()[_node][0], self.GetShapeGradient()[_node][1], self.GetShapeGradient()[_node][2]])
        results_file["ResponseGradients"] = _response_gradients
        results_file.close()

    def __IsAnalysisCompleted(self):
        if os.path.isfile(self.log_file):
            with open(self.log_file, "r") as file_input:
                _lines = file_input.readlines()
            file_input.close()

            for i in range(0, len(_lines)):
                _line = _lines[-i-1].strip()
                if _line != "":
                    if _line == "KRATOS RESPONSE VALUE SUCCESSFULLY CALCULATED.":
                        return True
                    else:
                        return False
        return False



