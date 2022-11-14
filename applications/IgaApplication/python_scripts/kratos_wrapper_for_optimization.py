import KratosMultiphysics as KM
KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
import KratosMultiphysics.StructuralMechanicsApplication as KSM
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import numpy as np
import math

from time import time
import sys, os
sys.path.append(os.pardir) # needed to import from parent dir

verbose = True
max_iterations = 300
model_parts_to_compare = ["IgaModelPart"]

class KratosWrapperForOptimization:
    def IsThisModelPartPartOfTheOptimization(self, mp_name):
        return self.IsStructuralAnalysisModelPart(mp_name)

    def IsStructuralAnalysisModelPart(self, mp_name):
        return "StructuralAnalysis" in mp_name

    def CountNumberOfOptimizationModelParts(self):
        counter = 0
        for model_part in self.main_model_part.SubModelParts:
            if self.IsThisModelPartPartOfTheOptimization(model_part.Name):
                counter += 1
        return counter

    def __init__(self):
        with open("ProjectParameters.json",'r') as parameter_file:
            parameters = KM.Parameters(parameter_file.read())

        self.model = KM.Model()

        # prepare the analysis
        self.analysis = StructuralMechanicsAnalysis(self.model, parameters)
        self.analysis.Initialize()

        self.main_model_part = self.model["IgaModelPart"]

        self.analysis.time = self.analysis._GetSolver().AdvanceInTime(self.analysis.time)
        self.analysis.InitializeSolutionStep()
        self.analysis._GetSolver().Predict()

        materials_filename = "optimization_materials.json"
        if materials_filename != "":
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KM.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            read_shit = KM.ReadMaterialsUtility(self.model)
            read_shit.ReadMaterialsToModelPart(material_settings)

        print()
        init_solution = []
        for model_part in self.main_model_part.SubModelParts:
            if self.IsThisModelPartPartOfTheOptimization(model_part.Name):
                print(f'ModelPart "{model_part.Name}" is part of the optimization')
                # add material id
                init_solution.append(1)

        self.optimization_iteration = 0
        self.min_norm = 1e10
        self.norms = []

        print(f"\nNumber of design variables: {len(init_solution)}\n\n{init_solution}\n")

    def ComputeMaxDisplacements(self):
        max_disp = 0.0
        for mp_name in model_parts_to_compare:
            for element in self.model[mp_name].Elements:
                displacement = element.CalculateOnIntegrationPoints(KM.DISPLACEMENT, self.model[mp_name].ProcessInfo)
                disp = math.sqrt(displacement[0][0]*displacement[0][0] + displacement[0][1]*displacement[0][1] + displacement[0][2]*displacement[0][2])
                max_disp = max([max_disp, disp])
        return max_disp

    def UpdateMaterials(self, updated_solution):
        print("Updated solution in UpdateMaterials:" + str(updated_solution))
        counter = 0
        for model_part in self.main_model_part.SubModelParts:
            if self.IsThisModelPartPartOfTheOptimization(model_part.Name):
                properties = self.main_model_part.GetProperties()[int(updated_solution[counter])]
                #model_part.AddProperties(properties)

                counter += 1
                for element in model_part.Elements:
                    element.Properties = properties

    def ResetModel(self, model_part, analysis):
        # set displacements back to zero
        KM.VariableUtils().SetVariable(KM.DISPLACEMENT, [0,0,0], model_part.Nodes)

        # restore initial configuration
        KM.VariableUtils().UpdateCurrentToInitialConfiguration(model_part.Nodes)

        # this resets the elements to the initial state
        analysis._GetSolver().Initialize()

    def PrintTime(self, label, start_time):
        print(label+": {0:.{1}f} [s]".format(time()-start_time,2))


    def OptimizableFunction(self, updated_solution):
        self.optimization_iteration += 1

        if verbose:
            print(f"\nOptimization iteration: {self.optimization_iteration}")

        start_time = time()

        self.UpdateMaterials(updated_solution)
        self.ResetModel(self.main_model_part, self.analysis)

        is_converged = self.analysis._GetSolver().SolveSolutionStep()

        if verbose:
            num_nl_iter = self.main_model_part.ProcessInfo[KM.NL_ITERATION_NUMBER]
            if is_converged:
                print(f"    Solver converged, iterations: {num_nl_iter}")
            else:
                print("    Solver did not converge")

        max_disp = self.ComputeMaxDisplacements()
        msg = f"    Max displacement: {round(max_disp, 5)}"
        new_min_found = False
        if max_disp < self.min_norm:
            self.min_norm = max_disp
            new_min_found = True
            msg += " " + "new minimum found!"
        print(msg)

        self.norms.append(max_disp)

        self.PrintTime("    Time needed", start_time)

        if self.optimization_iteration >= max_iterations:
            raise StopIteration("Maximum number of iterations reached!")

        KM.Logger.Flush()
        sys.stdout.flush()

        return max_disp