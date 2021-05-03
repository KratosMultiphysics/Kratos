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
            "variable_to_reconstruct"     : "REACTION"
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        parameters.ValidateAndAssignDefaults(default_settings)
        # Alternative:
        # settings.RecursivelyValidateAndAssignDefaults(default_settings)
        if (default_settings["variable_to_reconstruct"].GetString())!="REACTION":
            raise Exception("Variable to reconstruct not yet implemented")
        
        self.model_part = Model.GetModelPart("Structure")
        
    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        computing_model_part = self.model_part
                
        #### Residuals ####
        ###################
        restricted_residual_parameters = KratosMultiphysics.Parameters("""
        {
            "nodal_unknowns" : ["REACTION_X","REACTION_Y","REACTION_Z"],
            "number_of_dofs" : 1
        }""" )
        self.ResidualUtilityObject_Restricted_Residuals = romapp.RomModelPartUtility(computing_model_part.GetSubModelPart("COMPUTE_HROM").GetSubModelPart("GENERIC_Interface_frame_1"), restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.ResidualUtilityObject_Main_Part = romapp.RomModelPartUtility(computing_model_part, restricted_residual_parameters, KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())
        self.hr_restricted_residual_elements = self.ResidualUtilityObject_Restricted_Residuals.GetElementListFromNode(computing_model_part)
        with open('RomParameters_Residuals.json') as s:
            hr_data_residuals = json.load(s) # Load SVD of Residuals
            hr_modes_size = hr_data_residuals["rom_settings"]["number_of_rom_dofs"]
            hr_dofs_per_element = hr_data_residuals["rom_settings"]["number_of_dofs_per_element"]
            hr_elements_size = len(hr_data_residuals["elemental_modes"])
            total_restricted_dofs = int(hr_dofs_per_element*hr_elements_size)
            self.u_residual = np.zeros((total_restricted_dofs,hr_modes_size))
            restricted_residual_elements = []
            counter_in = 0
            for key in hr_data_residuals["elemental_modes"].keys():
                counter_fin = counter_in+hr_dofs_per_element
                self.u_residual[counter_in:counter_fin,:] = np.array(hr_data_residuals["elemental_modes"][key])
                restricted_residual_elements.append(int(key))
                counter_in = counter_fin
            restricted_residual_elements_index = self.FindDuplicateIndex(restricted_residual_elements) #No need to be saved "self"
            hr_total_restricted_dofs = int(hr_dofs_per_element*len(self.hr_restricted_residual_elements))
            self.hr_u_residual = np.zeros((hr_total_restricted_dofs,hr_modes_size))
            counter_in = 0
            for i in restricted_residual_elements_index:
                counter_fin = counter_in+hr_dofs_per_element
                aux_index_in = int(i*hr_dofs_per_element)
                aux_index_fin = aux_index_in + hr_dofs_per_element
                self.hr_u_residual[counter_in:counter_fin,:] = self.u_residual[aux_index_in:aux_index_fin,:]
                counter_in = counter_fin
        ###################

        #### Reaction Assembling ####
        ###################
        with open('Reaction_parameters.json') as s:
            reaction_assembling = json.load(s)
        self.restricted_residual_nodes = reaction_assembling["List_Nodes"]
        self.aux_list_nodes = reaction_assembling["List_Nodal_Elements"]
        self.aux_list_elem = reaction_assembling["List_Elemental_Nodes"]
        self.dof_node = len(reaction_assembling["Nodal_Unknowns"])
        self.dof_elem = int(reaction_assembling["Nodes_per_element"]*self.dof_node)
        ###################

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        model_part = self.model_part
        #### Residuals #####
        ###################
        hr_restricted_residual = self.ResidualUtilityObject_Main_Part.GetNonProjectedResidualsFromElementList(self.hr_restricted_residual_elements) 
        hr_residuals = np.array(hr_restricted_residual).T
        B = np.linalg.pinv(self.hr_u_residual)@hr_residuals
        residual_reconstruction = self.u_residual@B
        ###################

        #### Assemble Reactions ####
        ###################
        self.ResidualUtilityObject_Main_Part.AssembleReactions(KratosMultiphysics.Vector(residual_reconstruction.tolist()), KratosMultiphysics.Vector(self.restricted_residual_nodes), self.aux_list_nodes, self.aux_list_elem, self.dof_elem, self.dof_node)
        ###################

    def FindDuplicateIndex(self,List_1):
        List_2 = self.hr_restricted_residual_elements
        Index_list = []
        for i in range(len(List_1)):
            for j in range(len(List_2)):
                if List_1[i]==List_2[j]:
                    Index_list.append(i)
                    continue
        return Index_list