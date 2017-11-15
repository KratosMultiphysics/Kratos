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
        
        # Getting the part corresponding to soil part for its activation in the initialize
        part_number = parameters["soil_part"].GetString()
        parameters.RemoveValue("soil_part")
        parameters.AddEmptyValue("mechanical_soil_part").SetString("Parts_Parts_Auto"+part_number[10:])
        parameters.AddEmptyValue("thermal_soil_part").SetString("Thermal_Part_Auto_"+part_number[10:])
        
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
            for j, line in enumerate(file_name2):
                file_2 = line.split(" ")
                if (len(file_2)) > 1:
                    self.Construction.AssignTimeActivation(file_2[1],int(file_2[3]),float(file_2[0]))

    def InitializeSolutionStep(self):

        time = self.mechanical_model_part.ProcessInfo[TIME]
        tol = self.mechanical_model_part.ProcessInfo[DELTA_TIME]*1e-10
        
        # Activation according the input file
        with open(self.construction_input_file_name,'r') as file_name3:
            for j, line in enumerate(file_name3):
                file_3 = line.split(" ")
                if ((len(file_3)) > 1 and ((time <=(float(file_3[0])+tol)) and (time >=(float(file_3[0])-tol)))):
                    print("New phase has been activated...")
                    self.Construction.InitializeSolutionStep(file_3[1],file_3[2],int(file_3[3]))

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