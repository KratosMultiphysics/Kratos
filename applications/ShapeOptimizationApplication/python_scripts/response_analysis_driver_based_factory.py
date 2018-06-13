# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
#
# ==============================================================================
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
import os, sys, localimport

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
        self.results_file = "%s.results" % self.analysis_driver_name
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
            #  TODO: Transfer iteration specific log entries to a file in the iteration folder to seperate them from optimization logging entries
            self.analysis_driver.Run(self.model_part, self.response_data)

    def GetValue(self):
        return self.response_data["response_value"]

    def GetShapeGradient(self):
        return self.response_data["response_gradients"]           

    def FinalizeSolutionStep(self):
        if not self.is_analysis_step_completed:
            # write the final results to have restart capabilities
            self.__WriteResults()

    def __ReadResults(self):
        with open(self.results_file, "r") as results_file:
            lines = results_file.readlines()
        results_file.close()

        self.response_data = {}
        found_gradients = False
        index = 0
        while (index < len(lines)):
            line = lines[index]
            if line=="Response Value:\n":
                self.response_data["response_value"] = float(lines[index+1][:-1])
                index+=1

            if line=="Response Gradients:\n":
                found_gradients = True
                index += 2
                self.response_data["response_gradients"] = {}
                continue
            
            if line=="KRATOS RESPONSE VALUE SUCCESSFULLY CALCULATED.":
                found_gradients==False
                index += 1
                continue

            if found_gradients:
                line = lines[index].strip().split()
                _node = int(line[0])
                _gradient = Vector(3)
                _gradient[0] = float(line[1])
                _gradient[1] = float(line[2])
                _gradient[2] = float(line[3])
                self.response_data["response_gradients"][_node] = _gradient
            
            index += 1

    def __WriteResults(self):
        results_file = open(self.results_file, "w")
        results_file.write("Response Value:\n")
        results_file.write("%23.18e" % self.GetValue())
        results_file.write("\nResponse Gradients:\n")
        results_file.write("  NodeId       ResponseGradient_X       ResponseGradient_Y       ResponseGradient_Z")
        for _node in self.GetShapeGradient():
            results_file.write("\n%8d %23.18e %23.18e %23.18e" % (_node, self.GetShapeGradient()[_node][0], self.GetShapeGradient()[_node][1], self.GetShapeGradient()[_node][2]))
        results_file.write("\nKRATOS RESPONSE VALUE SUCCESSFULLY CALCULATED.")
        results_file.close()

    def __IsAnalysisCompleted(self):
        if os.path.isfile(self.results_file):
            with open(self.results_file, "r") as file_input:
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



