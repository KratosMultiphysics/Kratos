# importing the Kratos Library
from numpy import gradient
from . import base_response
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.OptimizationApplication.responses.base_response import BaseResponseFunction
import KratosMultiphysics.StructuralMechanicsApplication as KSM
from KratosMultiphysics.IgaApplication.map_nurbs_volume_results_to_embedded_geometry_process import MapNurbsVolumeResultsToEmbeddedGeometryProcess

import time as timer
import numpy as np

# ==============================================================================
class StressResponseFunction(BaseResponseFunction):
    """Stress response function. It triggers the primal analysis and
    uses the primal analysis results to calculate response value and gradient.

    Attributes
    ----------

    """

    def __init__(self,response_name, response_settings,response_analysis,model):

        self.type = "stress"
        self.variable = "STRESS"
        super().__init__(response_name, response_settings, model, response_analysis)

        if not self.response_settings.Has("gradient_settings"):
            self.gradient_settings = KM.Parameters()
            self.gradient_settings.AddString("gradient_mode","semi_analytic")
            self.gradient_settings.AddDouble("step_size",1e-6)
        else:
            self.gradient_settings = self.response_settings["gradient_settings"]

        self.supported_control_types = ["shape","material"]
        self.gradients_variables = {"shape":"D_STRESS_D_X","material":"D_STRESS_D_FD"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("StressResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("StressResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        # add vars
        for control_type in self.control_types:
            self.analysis_model_part.AddNodalSolutionStepVariable(KSM.ADJOINT_DISPLACEMENT)
            self.analysis_model_part.AddNodalSolutionStepVariable(KOA.ADJOINT_RHS)
            if control_type == "shape":
                self.response_settings["gradient_settings"].AddString("shape_gradient_field_name",self.gradients_variables[control_type])
                self.analysis_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))
            if control_type == "material":
                self.response_settings["gradient_settings"].AddString("material_gradient_field_name",self.gradients_variables[control_type])
                self.analysis_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

        ## Construct the linear solver
        import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
        self.adj_solver_settings = KM.Parameters("""{
                        "solver_type" : "amgcl",
                        "smoother_type":"ilu0",
                        "krylov_type": "gmres",
                        "coarsening_type": "aggregation",
                        "max_iteration": 200,
                        "provide_coordinates": false,
                        "gmres_krylov_space_dimension": 200,
                        "verbosity" : 0,
                        "tolerance": 1e-7,
                        "scaling": false,
                        "block_size": 1,
                        "use_block_matrices_if_possible" : true,
                        "coarse_enough" : 5000
                }""")
        self.linear_solvers = []
        root_model_parts = []
        for model_part_name in self.evaluated_model_parts:
            extracted_root_model_part_name = model_part_name.split(".")[0]
            if not extracted_root_model_part_name in root_model_parts:
                root_model_parts.append(extracted_root_model_part_name)
                self.linear_solvers.append(python_linear_solver_factory.ConstructSolver(self.adj_solver_settings))

        # create response
        self.response_function = KOA.StressOptResponse(response_name,model,self.response_settings,self.linear_solvers)

    def GetVariableName(self):
        return  self.variable

    def GetGradientsVariablesName(self):
        return self.gradients_variables

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("StressResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("StressResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("StressResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
        return self.value

    def GetValue(self):
        self.value = self.response_function.CalculateValue()
        return self.value

    def CalculateGradients(self):
        Logger.PrintInfo("StressResponseFunction", "Starting gradient calculation for response ", self.name)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("StressResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

    def CalculateGradientsForTypesAndObjects(self,control_types,controlled_objects,raise_error=True):

        if raise_error:
            for itr in range(len(controlled_objects)):
                controlled_object = controlled_objects[itr]
                control_type = control_types[itr]
                found = False
                for itr_2 in range(len(self.controlled_model_parts)):
                    controlled_model_part = self.controlled_model_parts[itr_2]
                    controlled_type = self.control_types[itr_2]
                    if controlled_type==control_type and controlled_model_part==controlled_object:
                        found = True
                        break
                if not found:
                    raise RuntimeError("StressResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("StressResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("StressResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

# ==============================================================================
class StrainEnergyResponseFunction(BaseResponseFunction):
    """Linear strain energy response function. It triggers the primal analysis and
    uses the primal analysis results to calculate response value and gradient.

    Attributes
    ----------
    primal_model_part : Model part of the primal analysis object
    primal_analysis : Primal analysis object of the response function
    response_function: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self,response_name, response_settings,response_analysis,model):

        self.type = "strain_energy"
        super().__init__(response_name, response_settings, model, response_analysis)

        if not self.response_settings.Has("gradient_settings"):
            self.gradient_settings = KM.Parameters()
            self.gradient_settings.AddString("gradient_mode","semi_analytic")
            self.gradient_settings.AddDouble("step_size",1e-6)
        else:
            self.gradient_settings = self.response_settings["gradient_settings"]

        if not self.analysis_model_part.HasNodalSolutionStepVariable(KM.KratosGlobals.GetVariable("D_STRAIN_ENERGY_1_D_X")):
            self.variable = "STRAIN_ENERGY_1"
        elif not self.analysis_model_part.HasNodalSolutionStepVariable(KM.KratosGlobals.GetVariable("D_STRAIN_ENERGY_2_D_X")):
            self.variable = "STRAIN_ENERGY_2"
        elif not self.analysis_model_part.HasNodalSolutionStepVariable(KM.KratosGlobals.GetVariable("D_STRAIN_ENERGY_3_D_X")):
            self.variable = "STRAIN_ENERGY_3"

        self.supported_control_types = ["shape","thickness","material"]
        self.gradients_variables = {"shape":"D_"+self.variable+"_D_X","thickness":"D_"+self.variable+"_D_FT","material":"D_"+self.variable+"_D_FD"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("StrainEnergyResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("StrainEnergyResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        # add vars
        for control_type in self.control_types:
            if control_type == "shape":
                self.response_settings["gradient_settings"].AddString("shape_gradient_field_name",self.gradients_variables[control_type])
                self.analysis_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))
            if control_type == "thickness":
                self.response_settings["gradient_settings"].AddString("thickness_gradient_field_name",self.gradients_variables[control_type])
                self.analysis_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))
            if control_type == "material":
                self.response_settings["gradient_settings"].AddString("material_gradient_field_name",self.gradients_variables[control_type])
                self.analysis_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

        # create response
        self.response_function = KOA.LinearStrainEnergyOptResponse(response_name,model,self.response_settings)

    def GetVariableName(self):
        return  self.variable

    def GetGradientsVariablesName(self):
        return self.gradients_variables

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("StrainEnergyResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("StrainEnergyResponse:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("StrainEnergyResponse:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
        return self.value

    def GetValue(self):
        self.value = self.response_function.CalculateValue()
        return self.value

    def CalculateGradients(self):
        Logger.PrintInfo("StrainEnergyResponse", "Starting gradient calculation for response ", self.name)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

    def CalculateGradientsForTypesAndObjects(self,control_types,controlled_objects,raise_error=True):

        if raise_error:
            for itr in range(len(controlled_objects)):
                controlled_object = controlled_objects[itr]
                control_type = control_types[itr]
                found = False
                for itr_2 in range(len(self.controlled_model_parts)):
                    controlled_model_part = self.controlled_model_parts[itr_2]
                    controlled_type = self.control_types[itr_2]
                    if controlled_type==control_type and controlled_model_part==controlled_object:
                        found = True
                        break
                if not found:
                    raise RuntimeError("StrainEnergyResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("StrainEnergyResponse", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")


class MassResponseFunction(BaseResponseFunction):

    def __init__(self,response_name, response_settings,model):

        self.type = "mass"
        self.variable = "MASS"
        super().__init__(response_name, response_settings, model)

        if not self.response_settings.Has("gradient_settings"):
            self.gradient_settings = KM.Parameters()
            self.gradient_settings.AddString("gradient_mode","semi_analytic")
            self.gradient_settings.AddDouble("step_size",1e-6)
        else:
            self.gradient_settings = self.response_settings["gradient_settings"]

        self.supported_control_types = ["shape","thickness","material"]
        self.gradients_variables = {"shape":"D_MASS_D_X","thickness":"D_MASS_D_FT","material":"D_MASS_D_FD"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("MassResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("MassResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))


        root_model_part_name = self.evaluated_model_parts[0].split(".")[0]
        for evaluated_model_part in self.evaluated_model_parts:
            if evaluated_model_part.split(".")[0] != root_model_part_name:
                raise RuntimeError("MassResponseFunction: evaluated_model_parts of mass response must have the same root model part !")

        self.root_model_part = self.model.GetModelPart(root_model_part_name)

        # add vars and response
        self.response_function = KOA.MassOptResponse(response_name,model,self.response_settings)
        for control_type in self.control_types:
            self.root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))


    def GetVariableName(self):
        return  self.variable

    def GetGradientsVariablesName(self):
        return self.gradients_variables

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("MassResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("MassResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("MassResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
        return self.value

    def CalculateGradientsForTypesAndObjects(self,control_types,controlled_objects,raise_error=True):

        if raise_error:
            for itr in range(len(controlled_objects)):
                controlled_object = controlled_objects[itr]
                control_type = control_types[itr]
                found = False
                for itr_2 in range(len(self.controlled_model_parts)):
                    controlled_model_part = self.controlled_model_parts[itr_2]
                    controlled_type = self.control_types[itr_2]
                    if controlled_type==control_type and controlled_model_part==controlled_object:
                        found = True
                        break
                if not found:
                    raise RuntimeError("MassResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("MassResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("MassResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

# ==============================================================================
class EmbeddedStrainEnergyResponseFunction(BaseResponseFunction):
    """Linear strain energy response function. It uses the primal analysis
    results to calculate response value and gradient.

    Attributes
    ----------
    primal_model_part : Model part of the primal analysis object
    primal_analysis : Primal analysis object of the response function
    response_function: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self,response_name, response_settings,response_analysis,model):

        self.type = "strain_energy"
        super().__init__(response_name, response_settings, model, response_analysis)

        if not self.response_settings.Has("gradient_settings"):
            self.gradient_settings = KM.Parameters()
            self.gradient_settings.AddString("gradient_mode","semi_analytic")
            self.gradient_settings.AddDouble("step_size",1e-6)
        else:
            self.gradient_settings = self.response_settings["gradient_settings"]

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("EmbeddedStrainEnergyResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        if len(self.controlled_model_parts) != 1:
            raise RuntimeError("EmbeddedStrainEnergyResponseFunction: 'controlled_objects' of response '{}' must have only one entry !".format(self.name))


        self.embedded_model_part = self.model.GetModelPart(self.controlled_model_parts[0])
        self.embedded_model_part_name = self.controlled_model_parts[0]
        self.nurbs_model_part = self.model.GetModelPart(self.evaluated_model_parts[0])
        self.nurbs_model_part_name = self.evaluated_model_parts[0]

        if not self.embedded_model_part.HasNodalSolutionStepVariable(KM.KratosGlobals.GetVariable("D_STRAIN_ENERGY_1_D_X")):
            self.variable = "STRAIN_ENERGY_1"
        elif not self.embedded_model_part.HasNodalSolutionStepVariable(KM.KratosGlobals.GetVariable("D_STRAIN_ENERGY_2_D_X")):
            self.variable = "STRAIN_ENERGY_2"
        elif not self.embedded_model_part.HasNodalSolutionStepVariable(KM.KratosGlobals.GetVariable("D_STRAIN_ENERGY_3_D_X")):
            self.variable = "STRAIN_ENERGY_3"

        self.supported_control_types = ["shape"]
        self.gradients_variables = {"shape":"D_"+self.variable+"_D_X"}
        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("EmbeddedStrainEnergyResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        # add vars
        for control_type in self.control_types:
            if control_type == "shape":
                self.response_settings["gradient_settings"].AddString("shape_gradient_field_name",self.gradients_variables[control_type])
                self.embedded_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))
                self.nurbs_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

        # create embedded response
        self.embedded_response_settings = self.response_settings.Clone()
        self.embedded_response_settings["controlled_objects"][0].SetString(self.embedded_response_settings["evaluated_objects"][0].GetString())
        self.response_function = KOA.LinearStrainEnergyOptResponse(response_name,model,self.embedded_response_settings)

    def GetVariableName(self):
        return  self.variable

    def GetGradientsVariablesName(self):
        return self.gradients_variables

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("EmbeddedStrainEnergyResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("EmbeddedStrainEnergyResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("EmbeddedStrainEnergyResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
        return self.value

    def GetValue(self):
        self.value = self.response_function.CalculateValue()
        return self.value

    def CalculateGradients(self):
        pass

    def CalculateGradientsForTypesAndObjects(self,control_types,controlled_objects,raise_error=True):

        if raise_error:
            for itr in range(len(controlled_objects)):
                controlled_object = controlled_objects[itr]
                control_type = control_types[itr]
                found = False
                for itr_2 in range(len(self.controlled_model_parts)):
                    controlled_model_part = self.controlled_model_parts[itr_2]
                    controlled_type = self.control_types[itr_2]
                    if controlled_type==control_type and controlled_model_part==controlled_object:
                        found = True
                        break
                if not found:
                    raise RuntimeError("EmbeddedStrainEnergyResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("EmbeddedStrainEnergyResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()

        # compute shape sens on the NurbsMesh
        self.response_function.CalculateGradient()

        # Map shape sensitivities from nurbs geometry to embedded geometry
        map_params = KM.Parameters(
        """ {
                "main_model_part_name"                    : \""""+self.nurbs_model_part_name+"""\",
                "nurbs_volume_name"                       : "NurbsVolume",
                "embedded_model_part_name"                : \""""+self.embedded_model_part_name+"""\",
                "nodal_results": [\""""+self.gradients_variables["shape"]+"""\"]
        } """ )
        process = MapNurbsVolumeResultsToEmbeddedGeometryProcess(self.model, map_params)
        process.ExecuteBeforeOutputStep()
        Logger.PrintInfo("EmbeddedStrainEnergyResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

# ==============================================================================
class EmbeddedMassResponseFunction(BaseResponseFunction):
    """Linear strain energy response function. It uses the primal analysis
    results to calculate response value and gradient.

    Attributes
    ----------
    primal_model_part : Model part of the primal analysis object
    primal_analysis : Primal analysis object of the response function
    response_function: Cpp utilities object doing the actual computation of response value and gradient.
    """

    def __init__(self,response_name, response_settings,model):

        self.type = "mass"
        self.variable = "MASS"
        super().__init__(response_name, response_settings, model)

        if not self.response_settings.Has("gradient_settings"):
            self.gradient_settings = KM.Parameters()
            self.gradient_settings.AddString("gradient_mode","semi_analytic")
            self.gradient_settings.AddDouble("step_size",1e-6)
        else:
            self.gradient_settings = self.response_settings["gradient_settings"]

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("EmbeddedMassResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        if len(self.controlled_model_parts) != 1:
            raise RuntimeError("EmbeddedMassResponseFunction: 'controlled_objects' of response '{}' must have only one entry !".format(self.name))


        self.embedded_model_part = self.model.GetModelPart(self.controlled_model_parts[0])
        self.embedded_model_part_name = self.controlled_model_parts[0]
        self.nurbs_model_part = self.model.GetModelPart(self.evaluated_model_parts[0])
        self.nurbs_model_part_name = self.evaluated_model_parts[0]

        self.supported_control_types = ["shape"]
        self.gradients_variables = {"shape":"D_MASS_D_X"}
        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("EmbeddedMassResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        # add vars
        for control_type in self.control_types:
            if control_type == "shape":
                self.response_settings["gradient_settings"].AddString("shape_gradient_field_name",self.gradients_variables[control_type])
                self.embedded_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))
                self.nurbs_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

        # create embedded response
        self.embedded_response_settings = self.response_settings.Clone()
        self.embedded_response_settings["controlled_objects"][0].SetString(self.embedded_response_settings["evaluated_objects"][0].GetString())
        self.response_function = KOA.MassOptResponse(response_name,model,self.embedded_response_settings)

    def GetVariableName(self):
        return  self.variable

    def GetGradientsVariablesName(self):
        return self.gradients_variables

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("EmbeddedMassResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("EmbeddedMassResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("EmbeddedMassResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
        return self.value

    def GetValue(self):
        self.value = self.response_function.CalculateValue()
        return self.value

    def CalculateGradients(self):
        pass

    def CalculateGradientsForTypesAndObjects(self,control_types,controlled_objects,raise_error=True):

        if raise_error:
            for itr in range(len(controlled_objects)):
                controlled_object = controlled_objects[itr]
                control_type = control_types[itr]
                found = False
                for itr_2 in range(len(self.controlled_model_parts)):
                    controlled_model_part = self.controlled_model_parts[itr_2]
                    controlled_type = self.control_types[itr_2]
                    if controlled_type==control_type and controlled_model_part==controlled_object:
                        found = True
                        break
                if not found:
                    raise RuntimeError("EmbeddedMassResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("EmbeddedMassResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()

        # compute shape sens on the NurbsMesh
        self.response_function.CalculateGradient()

        # Map shape sensitivities from nurbs geometry to embedded geometry
        map_params = KM.Parameters(
        """ {
                "main_model_part_name"                    : \""""+self.nurbs_model_part_name+"""\",
                "nurbs_volume_name"                       : "NurbsVolume",
                "embedded_model_part_name"                : \""""+self.embedded_model_part_name+"""\",
                "nodal_results": [\""""+self.gradients_variables["shape"]+"""\"]
        } """ )
        process = MapNurbsVolumeResultsToEmbeddedGeometryProcess(self.model, map_params)
        process.ExecuteBeforeOutputStep()
        Logger.PrintInfo("EmbeddedMassResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

class SelfIntersectionResponseFunction(BaseResponseFunction):

    def __init__(self,response_name, response_settings,model):

        self.type = "self_penetration"
        self.variable = "SELF_INTERSECT"
        super().__init__(response_name, response_settings, model)

        if not self.response_settings.Has("gradient_settings"):
            self.gradient_settings = KM.Parameters()
            self.gradient_settings.AddString("gradient_mode","semi_analytic")
            self.gradient_settings.AddDouble("step_size",1e-6)
        else:
            self.gradient_settings = self.response_settings["gradient_settings"]

        self.supported_control_types = ["shape"]
        self.gradients_variables = {"shape":"D_SELF_INTERSECT_D_X"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("SelfIntersectionResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("SelfIntersectionResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))


        root_model_part_name = self.evaluated_model_parts[0].split(".")[0]
        for evaluated_model_part in self.evaluated_model_parts:
            if evaluated_model_part.split(".")[0] != root_model_part_name:
                raise RuntimeError("SelfIntersectionResponseFunction: evaluated_model_parts of mass response must have the same root model part !")

        self.root_model_part = self.model.GetModelPart(root_model_part_name)

        # add vars
        for control_type in self.control_types:
            self.root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))


    def GetVariableName(self):
        return  self.variable

    def GetGradientsVariablesName(self):
        return self.gradients_variables

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("SelfIntersectionResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("SelfIntersectionResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        
        # run PyQuESo
        try:
            import QuESo_PythonApplication as QuESoApp
            from QuESo_PythonApplication.PyQuESo import PyQuESo
        except ImportError:
            raise Exception("QuESo python library is not available")

        pyqueso = PyQuESo("QUESOParameters.json")

        nodes = QuESoApp.PointVector()
        directions = QuESoApp.PointVector()
        for node in self.model.GetModelPart(self.evaluated_model_parts[0]).Nodes:
            nodal_normal = node.GetSolutionStepValue(KM.NORMAL)
            nodes.append( QuESoApp.Point(node.X0, node.Y0, node.Z0) )
            directions.append( QuESoApp.Point(-nodal_normal[0], -nodal_normal[1], -nodal_normal[2]) ) #This should be the direction of the shape

        pyqueso.Run()

        self.value = 0
        area = 0.0
        self.distances = pyqueso.ClosestDistances(nodes, directions)
        for node, distance in zip(self.model.GetModelPart(self.evaluated_model_parts[0]).Nodes, self.distances):
            nodal_area = node.GetSolutionStepValue(KM.NODAL_AREA)
            if abs(distance) < 20:
                self.value += (20-distance) * (20-distance) * nodal_area
            node.SetSolutionStepValue(KOA.AUXILIARY_FIELD, [distance, 0.0, 0.0])
            area += nodal_area

        Logger.PrintInfo("SelfIntersectionResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
        return self.value/area

    def CalculateGradientsForTypesAndObjects(self,control_types,controlled_objects,raise_error=True):

        if raise_error:
            for itr in range(len(controlled_objects)):
                controlled_object = controlled_objects[itr]
                control_type = control_types[itr]
                found = False
                for itr_2 in range(len(self.controlled_model_parts)):
                    controlled_model_part = self.controlled_model_parts[itr_2]
                    controlled_type = self.control_types[itr_2]
                    if controlled_type==control_type and controlled_model_part==controlled_object:
                        found = True
                        break
                if not found:
                    raise RuntimeError("SelfIntersectionResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("SelfIntersectionResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()

        for node, distance in zip(self.model.GetModelPart(self.evaluated_model_parts[0]).Nodes, self.distances):
            nodal_area = node.GetSolutionStepValue(KM.NODAL_AREA)
            nodal_normal = node.GetSolutionStepValue(KM.NORMAL)
            if abs(distance) < 20:
                sens = -1 * nodal_normal * (20-distance) * nodal_area
                node.SetSolutionStepValue(KOA.D_SELF_INTERSECT_D_X, sens)
            else:
                node.SetSolutionStepValue(KOA.D_SELF_INTERSECT_D_X, [0.0000,0.0000,0.0000])

        Logger.PrintInfo("SelfIntersectionResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")
