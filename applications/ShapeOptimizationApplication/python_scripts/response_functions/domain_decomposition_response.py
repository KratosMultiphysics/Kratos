# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Aditya Ghantasala, https://github.com/armingeiser
#
# ==============================================================================

import time as timer
import KratosMultiphysics as KM
from KratosMultiphysics import Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
try:
    import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
except ModuleNotFoundError:
    print("DomainDecompositionResponse requires StructuralMechanicsApplication. Make sure this is compiled.")

try:
    import KratosMultiphysics.IgaApplication as IgaApplication
except ModuleNotFoundError:
    print("DomainDecompositionResponse requires IgaApplication. Make sure this is compiled.")

def _GetModelPart(model, mp_name, dimension):
    if not model.HasModelPart(mp_name):
        model_part = model.CreateModelPart(mp_name, 2)
        if dimension < 0:
            raise Exception('Please specify a "domain_size" >= 0!')
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, dimension)
    else:
        model_part = model.GetModelPart(mp_name)

    return model_part


class DomainDecompositionResponse(ResponseFunctionInterface):
    """
    A class for optimizing the domain decomposition of a part based on different criteria.
    By default, stress along the specified line is used.

    Important settings:
    coming up soon
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        response_settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.response_settings = response_settings
        self.model = model

        self.cad_model_part = None
        self.trimming_curve = None

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = KM.Parameters(parameters_file.read())

        self.primal_analysis = StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_analysis_done = False

        self.value = None

        self.gradient = {}

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "response_type"             : "UNKNOWN_TYPE",
            "type"                      : "stress_aggregation",
            "cad_model_part_name"       : "UNKNOWN_NAME",
            "domain_size"               : 3,
            "cad_model_import_settings" : {
                "input_filename"    : "UNKNOWN_NAME"
            },
            "step_size"             : 1e-08 ,
            "trimming_curve_id"     : 10,
            "primal_settings"       : "primal_parameters.json"
        }""")
        return this_defaults

    def Initialize(self):
        self.primal_analysis.Initialize()
        if self.response_settings["cad_model_import_settings"]["input_type"].GetString() == "json":
            file_name = self.response_settings["cad_model_import_settings"]["input_filename"].GetString()
            mp_name = self.response_settings["cad_model_part_name"].GetString()
            dimension = self.response_settings["domain_size"].GetInt()
            self.cad_model_part = _GetModelPart(self.model, mp_name, dimension)
            ### Here import the cad model from the file.
            cad_model_part_io = KM.CadJsonInput(file_name)
            cad_model_part_io.ReadModelPart(self.cad_model_part)

            cad_geo = self.cad_model_part.Geometries[2]
            print("##################  ",cad_geo.Center())
            adsafsdf

            ### TODO: Get the trimming curve
        else:
            RuntimeError("Cad Geometry can only be imported from JSON file. https://github.com/orbingol/rw3dm can be helpful.")

    def InitializeSolutionStep(self):
        self.primal_analysis.time = self.primal_analysis._GetSolver().AdvanceInTime(self.primal_analysis.time)
        self.primal_analysis.InitializeSolutionStep()
        self.value = None
        self.gradient = {}

    def CalculateValue(self):
        Logger.PrintInfo("\n> Starting primal analysis for response", self.identifier)

        if not self.primal_analysis_done:
            startTime = timer.time()
            self.primal_analysis._GetSolver().Predict()
            self.primal_analysis._GetSolver().SolveSolutionStep()
            Logger.PrintInfo("DomainDecompositionResponse", "Time needed for solving the primal analysis",round(timer.time() - startTime,2),"s")
            self.primal_analysis_done = True

        startTime = timer.time()

        ## Calculate the value here

        self.value = 0.0

        Logger.PrintInfo("> Time needed for calculating the response value = ", round(timer.time() - startTime,2), "s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)

        startTime = timer.time()

        if not self.directions or not self.signed_distances:
            self._CalculateDistances()

        for i, node in enumerate(self.model_part.Nodes):
            gradient = [0.0,0.0,0.0] ## Calculate finite differences here !!
            self.gradient[node.Id] = gradient

        Logger.PrintInfo("> Time needed for calculating gradients = ", round(timer.time() - startTime,2), "s")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if not self.gradient:
            raise RuntimeError("Gradient was not calculated")
        if variable != KM.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        return self.gradient

    def FinalizeSolutionStep(self):
        pass

    def Finalize(self):
        if self.primal_analysis_done:
            self.primal_analysis.FinalizeSolutionStep()
            self.primal_analysis.OutputSolutionStep()
            self.primal_analysis.Finalize()