import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
import numpy as np
import json
from scipy.sparse import csr_matrix, load_npz
import time

def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return StressReconstructionProcess(model, parameters["Parameters"])

## All the processes python should be derived from "Process"

class StressReconstructionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, parameters):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "variable_to_reconstruct"     : "PK2_STRESS_VECTOR"
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        parameters.ValidateAndAssignDefaults(default_settings)
        # Alternative:
        # settings.RecursivelyValidateAndAssignDefaults(default_settings)
        if (parameters["variable_to_reconstruct"].GetString())!="PK2_STRESS_VECTOR":
            raise Exception("Variable to reconstruct not yet implemented")
        
        self.model_part = Model.GetModelPart("Structure")
        
    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.hr_elements = []
        self.hr_conditions = []
        ## Adding the weights to the corresponding elements
        with open('ElementsAndWeights.json') as f:
            hr_data = json.load(f)
            for key in hr_data["Elements"].keys():
                self.hr_elements.append(int(key))
            for key in hr_data["Conditions"].keys():
                self.hr_conditions.append(int(key))
            self.hr_elements.sort()
            self.hr_conditions.sort()
                
        #### Stresses ####
        ###################
        restricted_residual_parameters = KratosMultiphysics.Parameters("""
        {
            "nodal_unknowns" : ["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"],
            "number_of_dofs" : 1
        }""" )
        self.ResidualUtilityObject_Main_Part = romapp.RomModelPartUtility(self.model_part, restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.P = load_npz("shape.npz")
        with open('RomParameters.json') as s:
            rom_data = json.load(s) # Load SVD of Stresses
            pk2_dimension = 6 
            hr_modes_size = rom_data["stress_parameters"]["rom_settings"]["number_of_rom_modes"] # To initialize left singular values of Stresses
            FOM_gauss_point_size  = int(len(rom_data["stress_parameters"]["gauss_points_modes"])*pk2_dimension) # To initialize left singular values of Stresses
            self.u_stress = np.zeros((FOM_gauss_point_size ,hr_modes_size)) # Initialize left singular values of Stresses
            counter_in = 0
            for key in rom_data["stress_parameters"]["gauss_points_modes"].keys(): 
                counter_fin = counter_in + pk2_dimension
                self.u_stress[counter_in:counter_fin,:] = np.array(rom_data["stress_parameters"]["gauss_points_modes"][key]) # Create numpy array of left singular values
                counter_in = counter_fin
            hr_gauss_points = int(len(self.hr_elements)*pk2_dimension) 
            self.hr_u_stress = np.zeros((hr_gauss_points,hr_modes_size))
            counter_in = 0
            for i in self.hr_elements:
                counter_fin = counter_in + pk2_dimension
                gauss_point_in = (i)*pk2_dimension
                gauss_point_fin = gauss_point_in + pk2_dimension
                self.hr_u_stress[counter_in:counter_fin,:]=(self.u_stress[gauss_point_in:gauss_point_fin,:])
                counter_in = counter_fin
        ###################

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        model_part = self.model_part
        #### Stress #####
        ###################
        new_stress = []
        for elem in model_part.Elements:
            pk2_stress = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR,self.model_part.ProcessInfo)
            for i in range(pk2_stress[0].Size()):## Linear elements only have 1 Gauss-point
                new_stress.append(pk2_stress[0][i])
        q = np.linalg.pinv(self.hr_u_stress)@np.array(new_stress)
        stress = self.u_stress@q
        stress = self.ExtrapolateStressToNodes(stress)
        ###################

        #### Stress Assembling ####
        ###################
        self.ResidualUtilityObject_Main_Part.AssembleStresses(KratosMultiphysics.Matrix(stress))
        ###################

    def ExtrapolateStressToNodes(self,Vector):
        aux_shape_1 = int(np.shape(self.P)[1])
        aux_shape_2 = int(np.shape(Vector)[0])
        aux_shape_3 = int(aux_shape_2/aux_shape_1)
        extrapolated_matrix = np.zeros((aux_shape_1,aux_shape_3))
        for i in range(aux_shape_1):
            aux_index = i*aux_shape_3
            extrapolated_matrix[i,:] = Vector[aux_index:aux_index+aux_shape_3].T
        extrapolated_matrix = self.P@extrapolated_matrix
        return extrapolated_matrix