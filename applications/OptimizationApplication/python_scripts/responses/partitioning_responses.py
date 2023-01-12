import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.OptimizationApplication.responses.base_response import BaseResponseFunction

import time as timer
import numpy as np

class InterfaceResponseFunction(BaseResponseFunction):

    def __init__(self,response_name, response_settings,model):

        self.type = "interface"
        self.variable = "INTERFACE"
        super().__init__(response_name, response_settings, model)

        self.supported_control_types = ["material"]
        self.gradients_variables = {"material":"D_INTERFACE_D_FD"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("InterfaceResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("InterfaceResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))


        root_model_part_name = self.evaluated_model_parts[0].split(".")[0]
        for evaluated_model_part in self.evaluated_model_parts:
            if evaluated_model_part.split(".")[0] != root_model_part_name:
                raise RuntimeError("InterfaceResponseFunction: evaluated_model_parts of mass response must have the same root model part !")

        self.root_model_part = self.model.GetModelPart(root_model_part_name)

        # add vars and response
        self.response_function = KOA.InterfaceOptResponse(response_name,model,self.response_settings)
        for control_type in self.control_types:
            self.root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

    def GetVariableName(self):
        return  self.variable


    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("InterfaceResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("InterfaceResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("InterfaceResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
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
                    raise RuntimeError("InterfaceResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("InterfaceResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("InterfaceResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

class PartitionMassResponseFunction(BaseResponseFunction):

    def __init__(self,response_name, response_settings,model):

        self.type = "partition_mass"
        self.variable = "PARTITION_MASS"
        super().__init__(response_name, response_settings, model)

        self.supported_control_types = ["material"]
        self.gradients_variables = {"material":"D_PARTITION_MASS_D_FD"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("PartitionMassResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("PartitionMassResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))


        root_model_part_name = self.evaluated_model_parts[0].split(".")[0]
        for evaluated_model_part in self.evaluated_model_parts:
            if evaluated_model_part.split(".")[0] != root_model_part_name:
                raise RuntimeError("PartitionMassResponseFunction: evaluated_model_parts of mass response must have the same root model part !")

        self.root_model_part = self.model.GetModelPart(root_model_part_name)

        # add vars and response
        self.response_function = KOA.PartitionMassOptResponse(response_name,model,self.response_settings)
        for control_type in self.control_types:
            self.root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

    def GetVariableName(self):
        return  self.variable

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("PartitionMassResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("PartitionMassResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("PartitionMassResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
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
                    raise RuntimeError("PartitionMassResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("PartitionMassResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("PartitionMassResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

class PartitionMassResponseFunction(BaseResponseFunction):

    def __init__(self,response_name, response_settings,model):

        self.type = "partition_mass"
        self.variable = "PARTITION_MASS"
        super().__init__(response_name, response_settings, model)

        self.supported_control_types = ["material"]
        self.gradients_variables = {"material":"D_PARTITION_MASS_D_FD"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("PartitionMassResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("PartitionMassResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))


        root_model_part_name = self.evaluated_model_parts[0].split(".")[0]
        for evaluated_model_part in self.evaluated_model_parts:
            if evaluated_model_part.split(".")[0] != root_model_part_name:
                raise RuntimeError("PartitionMassResponseFunction: evaluated_model_parts of mass response must have the same root model part !")

        self.root_model_part = self.model.GetModelPart(root_model_part_name)

        # add vars and response
        self.response_function = KOA.PartitionMassOptResponse(response_name,model,self.response_settings)
        for control_type in self.control_types:
            self.root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

    def GetVariableName(self):
        return  self.variable

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("PartitionMassResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("PartitionMassResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("PartitionMassResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
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
                    raise RuntimeError("PartitionMassResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("PartitionMassResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("PartitionMassResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")

class PartitionInterfaceStressResponseFunction(BaseResponseFunction):

    def __init__(self,response_name, response_settings,model):

        self.type = "partition_interface_stress"
        self.variable = "STRESS"
        super().__init__(response_name, response_settings, model)

        self.supported_control_types = ["material"]
        self.gradients_variables = {"material":"D_STRESS_D_FD"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("PartitionInterfaceStressResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name))

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("PartitionInterfaceStressResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))


        root_model_part_name = self.evaluated_model_parts[0].split(".")[0]
        for evaluated_model_part in self.evaluated_model_parts:
            if evaluated_model_part.split(".")[0] != root_model_part_name:
                raise RuntimeError("PartitionInterfaceStressResponseFunction: evaluated_model_parts of mass response must have the same root model part !")

        self.root_model_part = self.model.GetModelPart(root_model_part_name)

        # add vars and response
        self.response_function = KOA.PartitionInterfaceStressOptResponse(response_name,model,self.response_settings)
        for control_type in self.control_types:
            self.root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

    def GetVariableName(self):
        return  self.variable

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("PartitionInterfaceStressResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types))

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()
        self.response_function.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("PartitionInterfaceStressResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = self.response_function.CalculateValue()
        Logger.PrintInfo("PartitionInterfaceStressResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")
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
                    raise RuntimeError("PartitionInterfaceStressResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("PartitionInterfaceStressResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()
        self.response_function.CalculateGradient()
        Logger.PrintInfo("PartitionInterfaceStressResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")