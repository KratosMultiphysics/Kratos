from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.DamApplication import *
CheckForPreviousImport()
from math import *

class ConditionsUtility:


    def __init__(self, delta_time, ConditionsOptions, model_part, time_unit_converter,evolution_type):

        self.time_unit_converter = time_unit_converter
        self.delta_time = delta_time

        self.imposed_displacement = ConditionsOptions.Imposed_Displacement
        self.imposed_pointload = ConditionsOptions.Imposed_PointLoad
        self.imposed_lineload = ConditionsOptions.Imposed_LineLoad
        self.imposed_surfaceload = ConditionsOptions.Imposed_SurfaceLoad
        self.imposed_normalload = ConditionsOptions.Imposed_NormalLoad
        self.imposed_tangentialload = ConditionsOptions.Imposed_TangentialLoad
        self.imposed_temperature = ConditionsOptions.Imposed_Temperature

        self.listofconditions=[]
        
        if(self.imposed_normalload=="Table"):
            if(evolution_type== "Exact"):
                self.listofconditions.append(ExactEvolutionConditionsLoadProcess(model_part, time_unit_converter))
            else:
                self.listofconditions.append(InterpolationEvolutionConditionsLoadProcess(model_part, time_unit_converter))
                
        if(self.imposed_temperature=="Table"):
            if(evolution_type== "Exact"):
                self.listofconditions.append(ExactEvolutionConditionsTemperatureProcess(model_part, time_unit_converter))
            else: 
                self.listofconditions.append(InterpolationEvolutionConditionsTemperatureProcess(model_part, time_unit_converter))
        elif(self.imposed_temperature=="Linearly_Incremented"):
            self.listofconditions.append(LinearEvolutionConditionsTemperatureProcess(model_part))

    def Initialize(self, model_part):

        self.NonLinearImposedConditions = False

        if(self.imposed_displacement=="Linearly_Incremented"):
            for node in model_part.Nodes:
                ImposedDisplacement = node.GetSolutionStepValue(IMPOSED_DISPLACEMENT)
                Displacement = node.GetSolutionStepValue(DISPLACEMENT)
                if(node.IsFixed(DISPLACEMENT_X)):
                    ImposedDisplacement[0] = Displacement[0]
                    Displacement[0] = 0
                if(node.IsFixed(DISPLACEMENT_Y)):
                    ImposedDisplacement[1] = Displacement[1]
                    Displacement[1] = 0
                if(node.IsFixed(DISPLACEMENT_Z)):
                    ImposedDisplacement[2] = Displacement[2]
                    Displacement[2] = 0
                node.SetSolutionStepValue(IMPOSED_DISPLACEMENT, ImposedDisplacement)
                node.SetSolutionStepValue(DISPLACEMENT, Displacement)
        elif(self.imposed_displacement=="Non-linearly_Modified"):
            self.NonLinearImposedConditions = True

        if(self.imposed_pointload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                ImposedPointLoad = node.GetSolutionStepValue(IMPOSED_POINT_LOAD)
                PointLoad = node.GetSolutionStepValue(POINT_LOAD)
                if(node.IsFixed(POINT_LOAD_X)):
                    ImposedPointLoad[0] = PointLoad[0]
                    PointLoad[0] = 0
                if(node.IsFixed(POINT_LOAD_Y)):
                    ImposedPointLoad[1] = PointLoad[1]
                    PointLoad[1] = 0
                if(node.IsFixed(POINT_LOAD_Z)):
                    ImposedPointLoad[2] = PointLoad[2]
                    PointLoad[2] = 0
                node.SetSolutionStepValue(IMPOSED_POINT_LOAD, ImposedPointLoad)
                node.SetSolutionStepValue(POINT_LOAD, PointLoad)
        elif(self.imposed_pointload=="Non-linearly_Modified"):
            self.NonLinearImposedConditions = True

        if(self.imposed_lineload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                ImposedLineLoad = node.GetSolutionStepValue(IMPOSED_LINE_LOAD)
                LineLoad = node.GetSolutionStepValue(LINE_LOAD)
                if(node.IsFixed(LINE_LOAD_X)):
                    ImposedLineLoad[0] = LineLoad[0]
                    LineLoad[0] = 0
                if(node.IsFixed(LINE_LOAD_Y)):
                    ImposedLineLoad[1] = LineLoad[1]
                    LineLoad[1] = 0
                if(node.IsFixed(LINE_LOAD_Z)):
                    ImposedLineLoad[2] = LineLoad[2]
                    LineLoad[2] = 0
                node.SetSolutionStepValue(IMPOSED_LINE_LOAD, ImposedLineLoad)
                node.SetSolutionStepValue(LINE_LOAD, LineLoad)
        elif(self.imposed_lineload=="Non-linearly_Modified"):
            self.NonLinearImposedConditions = True

        if(self.imposed_surfaceload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                ImposedSurfaceLoad = node.GetSolutionStepValue(IMPOSED_SURFACE_LOAD)
                SurfaceLoad = node.GetSolutionStepValue(SURFACE_LOAD)
                if(node.IsFixed(SURFACE_LOAD_X)):
                    ImposedSurfaceLoad[0] = SurfaceLoad[0]
                    SurfaceLoad[0] = 0
                if(node.IsFixed(SURFACE_LOAD_Y)):
                    ImposedSurfaceLoad[1] = SurfaceLoad[1]
                    SurfaceLoad[1] = 0
                if(node.IsFixed(SURFACE_LOAD_Z)):
                    ImposedSurfaceLoad[2] = SurfaceLoad[2]
                    SurfaceLoad[2] = 0
                node.SetSolutionStepValue(IMPOSED_SURFACE_LOAD, ImposedSurfaceLoad)
                node.SetSolutionStepValue(SURFACE_LOAD, SurfaceLoad)
        elif(self.imposed_surfaceload=="Non-linearly_Modified"):
            self.NonLinearImposedConditions = True

        if(self.imposed_normalload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                if(node.IsFixed(NORMAL_CONTACT_STRESS)):
                    NormalStress = node.GetSolutionStepValue(NORMAL_CONTACT_STRESS)
                    ImposedNormalStress = NormalStress
                    NormalStress = 0
                    node.SetSolutionStepValue(IMPOSED_NORMAL_STRESS, ImposedNormalStress)
                    node.SetSolutionStepValue(NORMAL_CONTACT_STRESS, NormalStress)
        elif(self.imposed_normalload=="Non-linearly_Modified"):
            self.NonLinearImposedConditions = True

        if(self.imposed_tangentialload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                if(node.IsFixed(TANGENTIAL_CONTACT_STRESS)):
                    TangentialStress = node.GetSolutionStepValue(TANGENTIAL_CONTACT_STRESS)
                    ImposedTangentialStress = TangentialStress
                    TangentialStress = 0
                    node.SetSolutionStepValue(IMPOSED_TANGENTIAL_STRESS, ImposedTangentialStress)
                    node.SetSolutionStepValue(TANGENTIAL_CONTACT_STRESS, TangentialStress)
        elif(self.imposed_tangentialload=="Non-linearly_Modified"):
            self.NonLinearImposedConditions = True

        if(self.imposed_temperature=="Linearly_Incremented"):
            for node in model_part.Nodes:
                if(node.IsFixed(TEMPERATURE)):
                    Temperature = node.GetSolutionStepValue(TEMPERATURE)
                    ImposedTemperature = Temperature
                    Temperature = 0
                    node.SetSolutionStepValue(IMPOSED_TEMPERATURE, ImposedTemperature)
                    node.SetSolutionStepValue(TEMPERATURE, Temperature)
        elif(self.imposed_temperature=="Non-linearly_Modified"):
            self.NonLinearImposedConditions = True

    def UpdateImposedConditions(self, model_part, current_step):
            
        for condition in self.listofconditions:
            condition.Execute()
        
        if(self.NonLinearImposedConditions):
            print("Customize your conditions")
        
        '''
        # Example of interpolated evolution conditions in python. This is currently performed in a process in C++.
        
        delta = self.delta_time/self.time_unit_converter
        step = current_step
        current_time= delta*current_step    

        # Update Reference Temperature for evolution time
        model_part.ProcessInfo[REFERENCE_TEMPERATURE] = model_part.GetTable(4).GetValue(current_time)
         
        # Bofang Conditions 
        if(model_part.GetMesh(1)):
            for node in model_part.GetMesh(1).Nodes:
                # Fixed Annual Information
                direction = model_part.GetMesh(1)[GRAVITY_DIRECTION]
                coordinate_base = model_part.GetMesh(1)[COORDINATE_BASE_DAM]
                surface_temp = model_part.GetMesh(1)[SURFACE_TEMP]
                bottom_temp = model_part.GetMesh(1)[BOTTOM_TEMP]
                height_dam = model_part.GetMesh(1)[HEIGHT_DAM]
                amplitude = model_part.GetMesh(1)[AMPLITUDE]
                freq = model_part.GetMesh(1)[FREQUENCY]
                day = model_part.GetMesh(1)[DAY_MAXIMUM]

                #Information provided by tables 
                time = model_part.GetTable(1).GetValue(current_time)             # Time Bofang
                water_level = model_part.GetTable(2).GetValue(current_time)      # Water Level
                outer_temp =  model_part.GetTable(3).GetValue(current_time)      # Outer temperature

                if(direction=="X"):
                    Temperature = node.GetSolutionStepValue(TEMPERATURE)
                    aux = (coordinate_base + water_level) - node.X
                    if(aux >= 0): 
                        aux1 = ((bottom_temp-(surface_temp*exp(-0.04*height_dam)))/(1-(exp(-0.04*height_dam))))
                        Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(time-(day/30)-2.15+(1.30*exp(-0.085*aux)))))))
                    else:
                        Temperature = outer_temp
                    node.SetSolutionStepValue(TEMPERATURE, Temperature)						

                elif(direction=="Y"):
                    Temperature = node.GetSolutionStepValue(TEMPERATURE)
                    aux = (coordinate_base + water_level) - node.Y
                    if(aux >= 0): 
                        aux1 = ((bottom_temp-(surface_temp*exp(-0.04*height_dam)))/(1-(exp(-0.04*height_dam))))
                        Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(time-(day/30)-2.15+(1.30*exp(-0.085*aux)))))))
                    else:
                        Temperature = outer_temp
                    node.SetSolutionStepValue(TEMPERATURE, Temperature)	

                elif(direction=="Z"):
                    Temperature = node.GetSolutionStepValue(TEMPERATURE)
                    aux = (coordinate_base + water_level) - node.Z
                    if(aux >= 0): 
                        aux1 = ((bottom_temp-(surface_temp*exp(-0.04*height_dam)))/(1-(exp(-0.04*height_dam))))
                        Temperature = (aux1+((surface_temp-aux1)*(exp(-0.04*aux)))+(amplitude*(exp(-0.018*aux))*(cos(freq*(time-(day/30)-2.15+(1.30*exp(-0.085*aux)))))))
                    else:
                        Temperature = outer_temp
                    node.SetSolutionStepValue(TEMPERATURE, Temperature)	

        # Uniform Temperature
        if(model_part.GetMesh(2)):
            for node in model_part.GetMesh(2).Nodes:
                #Information provided by tables
                #Temperature = node.GetSolutionStepValue(TEMPERATURE)
                Temperature = model_part.GetTable(3).GetValue(current_time)        # Outer temperature
                node.SetSolutionStepValue(TEMPERATURE, Temperature)

        # Hydrostatic pressure
        if(model_part.GetMesh(3)):
            for node in model_part.GetMesh(3).Nodes:
                # Fixed Annual Information
                direction = model_part.GetMesh(3)[GRAVITY_DIRECTION]
                coordinate_base = model_part.GetMesh(3)[COORDINATE_BASE_DAM]
                spe_weight = model_part.GetMesh(3)[SPECIFIC_WEIGHT]
                #Information provided by tables
                water_level = model_part.GetTable(2).GetValue(current_time)        # Water Level
   
                if(direction=="X"):
                    NormalStress = node.GetSolutionStepValue(NORMAL_CONTACT_STRESS)	
                    ref_coord = coordinate_base + water_level
                    pressure = spe_weight*(ref_coord- node.X)
                    if(pressure > 0):
                        node.SetSolutionStepValue(NORMAL_CONTACT_STRESS, pressure)
                    else: 
                        node.SetSolutionStepValue(NORMAL_CONTACT_STRESS, 0.0)

                elif(direction=="Y"):
                    NormalStress = node.GetSolutionStepValue(NORMAL_CONTACT_STRESS)	
                    ref_coord = coordinate_base + water_level
                    pressure = spe_weight*(ref_coord- node.Y)
                    if(pressure > 0):
                        node.SetSolutionStepValue(NORMAL_CONTACT_STRESS, pressure)
                    else: 
                        node.SetSolutionStepValue(NORMAL_CONTACT_STRESS, 0.0)

                elif(direction=="Z"):
                    NormalStress = node.GetSolutionStepValue(NORMAL_CONTACT_STRESS)	
                    ref_coord = coordinate_base + water_level
                    pressure = spe_weight*(ref_coord- node.Z)
                    if(pressure > 0):
                        node.SetSolutionStepValue(NORMAL_CONTACT_STRESS, pressure)
                    else: 
                        node.SetSolutionStepValue(NORMAL_CONTACT_STRESS, 0.0)
        '''
                            
