import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM
import numpy as np
import json

def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ReactionReconstructionProcess(model, parameters["Parameters"])

## All the processes python should be derived from "Process"

class ReactionReconstructionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, parameters):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "variable_to_reconstruct"     : "REACTION",
            "restricted_model_part" : "Name_of_restricted_model_part"
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        parameters.ValidateAndAssignDefaults(default_settings)
        # Alternative:
        # settings.RecursivelyValidateAndAssignDefaults(default_settings)
        if (parameters["variable_to_reconstruct"].GetString())!="REACTION":
            raise Exception("Variable to reconstruct not yet implemented")
        
        self.model_part = Model.GetModelPart("Structure")
        self.restricted_model_part_name = parameters["restricted_model_part"].GetString()
        
    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        computing_model_part = self.model_part
        try:
            restricted_model_part = computing_model_part.GetSubModelPart("COMPUTE_HROM").GetSubModelPart(self.restricted_model_part_name)
        except:
            restricted_model_part = computing_model_part.GetSubModelPart(self.restricted_model_part_name)
                
        #### Residuals ####
        ###################
        restricted_residual_parameters = KratosMultiphysics.Parameters("""
        {
            "nodal_unknowns" : ["REACTION_X","REACTION_Y","REACTION_Z"],
            "number_of_dofs" : 1
        }""" )
        self.ResidualUtilityObject_Restricted_Residuals = romapp.RomModelPartUtility(restricted_model_part, restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.ResidualUtilityObject_Main_Part = romapp.RomModelPartUtility(computing_model_part, restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.rom_restricted_residual_elements = self.ResidualUtilityObject_Restricted_Residuals.GetElementListFromNode(computing_model_part)
        with open('RomParameters.json') as s:
            rom_data = json.load(s) # Load SVD of Residuals
            rom_modes_size = rom_data["residual_parameters"]["rom_settings"]["number_of_rom_dofs"]
            rom_dofs_per_element = rom_data["residual_parameters"]["rom_settings"]["number_of_dofs_per_element"]
            rom_elements_size = len(rom_data["residual_parameters"]["elemental_modes"])
            total_restricted_dofs = int(rom_dofs_per_element*rom_elements_size)
            self.u_residual = np.zeros((total_restricted_dofs,rom_modes_size))
            restricted_residual_elements = []
            counter_in = 0
            for key in rom_data["residual_parameters"]["elemental_modes"].keys():
                counter_fin = counter_in+rom_dofs_per_element
                self.u_residual[counter_in:counter_fin,:] = np.array(rom_data["residual_parameters"]["elemental_modes"][key])
                restricted_residual_elements.append(int(key))
                counter_in = counter_fin
            restricted_residual_elements_index = self.FindDuplicateIndex(restricted_residual_elements) # Finds the restricted elements of the hyper reduced model part
            rom_total_restricted_dofs = int(rom_dofs_per_element*len(self.rom_restricted_residual_elements))
            self.hr_u_residual = np.zeros((rom_total_restricted_dofs,rom_modes_size))
            counter_in = 0
            for i in restricted_residual_elements_index:
                counter_fin = counter_in+rom_dofs_per_element
                aux_index_in = int(i*rom_dofs_per_element)
                aux_index_fin = aux_index_in + rom_dofs_per_element
                self.hr_u_residual[counter_in:counter_fin,:] = self.u_residual[aux_index_in:aux_index_fin,:]
                counter_in = counter_fin
        ###################

        #### Reaction Assembling ####
        ###################
        self.restricted_residual_nodes = rom_data["reaction_parameters"]["List_Nodes"]
        self.aux_list_nodes = rom_data["reaction_parameters"]["List_Nodal_Elements"]
        self.aux_list_elem = rom_data["reaction_parameters"]["List_Elemental_Nodes"]
        self.dof_node = len(rom_data["reaction_parameters"]["Nodal_Unknowns"])
        self.dof_elem = int(rom_data["reaction_parameters"]["Nodes_per_element"]*self.dof_node)
        ###################

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        #### Residuals #####
        ###################
        hr_restricted_residual = self.ResidualUtilityObject_Main_Part.GetNonProjectedResidualsFromElementList(self.rom_restricted_residual_elements) 
        hr_residuals = np.array(hr_restricted_residual).T
        B = np.linalg.pinv(self.hr_u_residual)@hr_residuals
        residual_reconstruction = self.u_residual@B
        ###################

        #### Assemble Reactions ####
        ###################
        self.ResidualUtilityObject_Main_Part.AssembleReactions(KratosMultiphysics.Vector(residual_reconstruction.tolist()), KratosMultiphysics.Vector(self.restricted_residual_nodes), self.aux_list_nodes, self.aux_list_elem, self.dof_elem, self.dof_node)
        ###################

    def FindDuplicateIndex(self,List_1):
        List_2 = self.rom_restricted_residual_elements
        Index_list = []
        for i in range(len(List_1)):
            for j in range(len(List_2)):
                if List_1[i]==List_2[j]:
                    Index_list.append(i)
                    continue
        return Index_list