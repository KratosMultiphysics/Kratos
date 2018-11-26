from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.DamApplication import *

class DamConstructionUtility:

    def __init__(self,mechanical_model_part, thermal_model_part, parameters):

        self.mechanical_model_part = mechanical_model_part
        self.thermal_model_part = thermal_model_part

        # Getting neccessary files from constructive phase
        ambient_input_file_name = parameters["ambient_input_file_name"].GetString()
        self.construction_input_file_name = parameters["construction_input_file_name"].GetString()

        self.heat_source_parameters = Parameters("{}")
        if (parameters["source_type"].GetString() == "Adiabatic"):
            self.heat_source_type = "Noorzai"
            self.heat_source_parameters.AddValue("density",parameters["density"])
            self.heat_source_parameters.AddValue("specific_heat",parameters["specific_heat"])
            self.heat_source_parameters.AddValue("alpha",parameters["alpha"])
            self.heat_source_parameters.AddValue("t_max",parameters["tmax"])

        if (parameters["source_type"].GetString() == "NonAdiabatic"):
            self.heat_source_type = "Azenha"
            self.heat_source_parameters.AddValue("activation_energy",parameters["activation_energy"])
            self.heat_source_parameters.AddValue("gas_constant",parameters["gas_constant"])
            self.heat_source_parameters.AddValue("constant_rate",parameters["constant_rate"])
            self.heat_source_parameters.AddValue("alpha_initial",parameters["alpha_initial"])
            self.heat_source_parameters.AddValue("q_total",parameters["q_total"])
            self.heat_source_parameters.AddValue("A",parameters["A"])
            self.heat_source_parameters.AddValue("B",parameters["B"])
            self.heat_source_parameters.AddValue("C",parameters["C"])
            self.heat_source_parameters.AddValue("D",parameters["D"])

        # If there are parts alreadey existent (actived from the beginning of the model) getting these parts for its activation in the initialize. This parts can correspond to either the soil or the already built part of the dam.
        if (parameters["activate_soil_part"].GetBool()):
            name_soil_part = parameters["name_soil_part"].GetString()
            parameters.RemoveValue("name_soil_part")
            parameters.AddEmptyValue("mechanical_soil_part").SetString("sub_Parts_"+name_soil_part)
            parameters.AddEmptyValue("thermal_soil_part").SetString("sub_Thermal_"+name_soil_part)

        if (parameters["activate_existing_part"].GetBool()):
            name_existing_part = parameters["name_existing_part"].GetString()
            parameters.RemoveValue("name_existing_part")
            parameters.AddEmptyValue("mechanical_existing_part").SetString("sub_Parts_"+name_existing_part)
            parameters.AddEmptyValue("thermal_existing_part").SetString("sub_Thermal_"+name_existing_part)

        self.table_ambient = PiecewiseLinearTable()
        with open(ambient_input_file_name,'r') as file_name1:
            for j, line in enumerate(file_name1):
                file_1= line.split(" ")
                if (len(file_1)) > 1:
                    self.table_ambient.AddRow(float(file_1[0]), float(file_1[1]))

        # Construct the utility
        self.Construction = ConstructionUtility(self.mechanical_model_part,self.thermal_model_part, self.table_ambient, parameters)

    def Initialize(self):
        self.Construction.Initialize()

        # The function recieves the mame of submodel Part, the number of phase and the activation time
        print("Assigning time activation for each node")
        with open(self.construction_input_file_name,'r') as file_name2:
            #for line in reversed(list((file_name2))):
            for j, line in enumerate(file_name2):
                file_2 = line.split(" ")
                if (len(file_2)) > 1:
                    self.name_sub_thermal_part = "sub_Thermal_" + file_2[1]
                    self.Construction.AssignTimeActivation(self.name_sub_thermal_part,int(file_2[2]),float(file_2[0]), float(file_2[3]))

    def InitializeSolutionStep(self):

        time = self.mechanical_model_part.ProcessInfo[TIME]
        tol = self.mechanical_model_part.ProcessInfo[DELTA_TIME]*1e-10
        time_unit_converter = self.mechanical_model_part.ProcessInfo[TIME_UNIT_CONVERTER]

        # Activation according the input file
        with open(self.construction_input_file_name,'r') as file_name3:
            for j, line in enumerate(file_name3):
                file_3 = line.split(" ")
                if ((len(file_3)) > 1 and ((time <=(float(file_3[0])*time_unit_converter+tol)) and (time >=(float(file_3[0])*time_unit_converter-tol)))):
                    print("New phase has been activated...")
                    self.name_sub_thermal_part = "sub_Thermal_" + file_3[1]
                    self.name_sub_mechanical_part = "sub_Parts_" + file_3[1]
                    self.Construction.InitializeSolutionStep(self.name_sub_thermal_part,self.name_sub_mechanical_part,int(file_3[2]))

        # Detection of fluxes (it must be done each step)
        print("Searching free surfaces...")
        self.Construction.SearchingFluxes()

        # Contribution of heat source
        print("Assigning heat source...")
        if (self.heat_source_type == 'Noorzai'):
            self.Construction.ActiveHeatFluxNoorzai(self.heat_source_parameters)

        elif (self.heat_source_type == 'Azenha'):
            self.Construction.ActiveHeatFluxAzenha(self.heat_source_parameters)


    def AfterOutputStep(self):
        self.Construction.AfterOutputStep()
