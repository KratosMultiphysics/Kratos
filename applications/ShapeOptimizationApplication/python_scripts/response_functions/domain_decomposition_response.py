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
import math
import KratosMultiphysics as KM
from KratosMultiphysics import Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface
from KratosMultiphysics import ShapeOptimizationApplication
try:
    import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
except ModuleNotFoundError:
    print("DomainDecompositionResponse requires StructuralMechanicsApplication. Make sure this is compiled.")

from .import cad_dd_utilities as cad_util

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
        print(response_settings)

        response_settings.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        self.response_settings = response_settings
        self.model = model

        self.evaluation_model_part = None
        self.trimming_curve = None

        with open(response_settings["primal_settings"].GetString()) as parameters_file:
            ProjectParametersPrimal = KM.Parameters(parameters_file.read())

        self.primal_analysis = StructuralMechanicsAnalysis(model, ProjectParametersPrimal)
        self.primal_analysis_done = False

        self.value = None

        self.gradient = {}

        self.updates = None

        self.sa_object = None # Stress Aggregation

        self.step = 0

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "response_type"             : "UNKNOWN_TYPE",
            "type"                      : "stress_aggregation",
            "domain_size"               : 3,
            "cad_model_io_settings" : {
                "type"              : "json",
                "input_filename"    : "UNKNOWN_NAME",
                "output_filename"   : "UNKNOWN_NAME"
            },
            "step_size"             : 1e-08 ,
            "trimming_curve_index"     : 0,
            "averaging_parameter"       :1.0,
            "primal_settings"       : "primal_parameters.json"
        }""")
        return this_defaults

    def Initialize(self):
        self.primal_analysis.Initialize()
        if self.response_settings["cad_model_io_settings"]["type"].GetString() == "json":
            file_name = self.response_settings["cad_model_io_settings"]["input_filename"].GetString()
            self.cad_geom = cad_util.GetNurbsGeometry(file_name)
            self.trimming_curve = cad_util.GetTrimmingCurve(self.cad_geom, self.response_settings["trimming_curve_index"].GetInt())
        else:
            RuntimeError("Cad Geometry can only be imported from JSON file. https://github.com/orbingol/rw3dm can be helpful.")

    def InitializeSolutionStep(self):
        self.__CreateModelPart()
        if not self.primal_analysis_done:
            self.primal_analysis.time = self.primal_analysis._GetSolver().AdvanceInTime(self.primal_analysis.time)
            self.primal_analysis.InitializeSolutionStep()
        self.value = None
        self.gradient = {} # node 1 [u,v] node2 [u, v] .... in paremetric space
        evl_pts = self.trimming_curve.evalpts
        cad_util.MakeModelPart(self.cad_geom, evl_pts, self.evaluation_model_part)
        KM.VtkOutput(self.evaluation_model_part).PrintOutput("trim_curve_points_"+str(self.step))
        self.step = self.step+1

    def CalculateValue(self):
        if not self.primal_analysis_done:
            Logger.PrintInfo("\n> Starting primal analysis for response. This will be done only once.", self.identifier)
            startTime = timer.time()
            self.primal_analysis._GetSolver().Predict()
            self.primal_analysis._GetSolver().SolveSolutionStep()
            Logger.PrintInfo("DomainDecompositionResponse", "Time needed for solving the primal analysis",round(timer.time() - startTime,2),"s")
            self.primal_analysis_done = True

        startTime = timer.time()
        rho = self.response_settings["averaging_parameter"].GetDouble()
        ## Calculate the value here using Kreisselmeier–Steinhauser (KS) Function
        ## value = (1/rho) * log(sigma(exp(rho*f(x))))
        ## f(x) here is the stress value at the control points of the trimming curve.
        ctrlpts = cad_util.GetControlPoints(self.trimming_curve)
        [physical_coordinates, derivatives] = cad_util.GetPointCoordinatesAndDerivatives(self.cad_geom, ctrlpts, 1)
        self.ks_util = ShapeOptimizationApplication.KreisselmeierSteinhauserAggregationUtility(self.primal_analysis._GetSolver().GetComputingModelPart())

        self.value = self.ks_util.CalculateValue(self.evaluation_model_part, StructuralMechanicsApplication.VON_MISES_STRESS, 0.1)

        Logger.PrintInfo("> Time needed for calculating the response value = ", round(timer.time() - startTime,2), "s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)
        startTime = timer.time()

        original_ctrlpts = cad_util.GetControlPoints(self.trimming_curve)
        original_ctrlpts_cpy = cad_util.GetControlPoints(self.trimming_curve)

        self.gradient = self.__GetGradients(original_ctrlpts.copy())
        self.trimming_curve.ctrlpts = original_ctrlpts

        # for i, ctrlpt in enumerate(original_ctrlpts):
        #     self.gradient[i] = [0.01,0.01]

        Logger.PrintInfo("> Time needed for calculating gradients = ", round(timer.time() - startTime,2), "s")

    def GetValue(self):
        return self.value

    def GetGradient(self):
        if not self.gradient:
            raise RuntimeError("Gradient was not calculated")
        return self.gradient

    def SetUpdates(self, updates):
        self.updates = updates

    def FinalizeSolutionStep(self):
        # Here update the u and v values.
        ctrlpts = cad_util.GetControlPoints(self.trimming_curve)
        for i, ctrlpt in enumerate(ctrlpts):
            ctrlpts[i] = [x + y for x, y in zip(ctrlpt, self.updates[i])]

            # u and v should be with in [0,1]
            if(ctrlpts[i][0] > 1):
                ctrlpts[i][0] = 1
            if(ctrlpts[i][0] < 0):
                ctrlpts[i][0] = 0

            if(ctrlpts[i][1] > 1):
                ctrlpts[i][1] = 1
            if(ctrlpts[i][1] < 0):
                ctrlpts[i][1] = 0

        self.trimming_curve.ctrlpts = ctrlpts

    def Finalize(self):
        if self.primal_analysis_done:
            self.primal_analysis.FinalizeSolutionStep()
            self.primal_analysis.OutputSolutionStep()
            self.primal_analysis.Finalize()
        cad_util.OutputCadToJson(self.cad_geom, self.response_settings["cad_model_io_settings"]["output_filename"].GetString())
        cad_util.VisualizeSurface(self.cad_geom)

    def __CreateModelPart(self, mp_name="dd_evaluation_model_part"):
        if self.model.HasModelPart(mp_name):
            self.model.DeleteModelPart(mp_name)
            self.evaluation_model_part = self.model.CreateModelPart(mp_name, 2)
        else:
            self.evaluation_model_part = self.model.CreateModelPart(mp_name, 2)

    def __GetGradients(self, points):
        delta = 1e-2
        grad = {}
        for i, point in enumerate(points):
            self.__CreateModelPart()
            u, v = points[i]
            cad_util.MakeModelPart(self.cad_geom, self.trimming_curve.evalpts, self.evaluation_model_part)
            u_val1 = self.ks_util.CalculateValue(self.evaluation_model_part, StructuralMechanicsApplication.VON_MISES_STRESS, 0.1)
            if not u == 1.0:
                points[i][0] = u+delta
            self.trimming_curve.ctrlpts = points
            cad_util.MakeModelPart(self.cad_geom, self.trimming_curve.evalpts, self.evaluation_model_part)
            u_val2 = self.ks_util.CalculateValue(self.evaluation_model_part, StructuralMechanicsApplication.VON_MISES_STRESS, 0.1)
            grad_u = -1*(u_val1-u_val2)/delta

            points[i][0] = u ## Reset u and now change v
            cad_util.MakeModelPart(self.cad_geom, self.trimming_curve.evalpts, self.evaluation_model_part)
            v_val1 = self.ks_util.CalculateValue(self.evaluation_model_part, StructuralMechanicsApplication.VON_MISES_STRESS, 0.1)
            if not v == 1.0:
                points[i][1] = v+delta
            self.trimming_curve.ctrlpts = points
            cad_util.MakeModelPart(self.cad_geom, self.trimming_curve.evalpts, self.evaluation_model_part)
            v_val2 = self.ks_util.CalculateValue(self.evaluation_model_part, StructuralMechanicsApplication.VON_MISES_STRESS, 0.1)
            grad_v = -1*(v_val1-v_val2)/delta

            points[i][1] = v ## Reset v
            norm = math.sqrt(grad_u*grad_u + grad_v*grad_v)
            grad[i] = [ grad_u/norm, grad_v/norm ]

        return grad

