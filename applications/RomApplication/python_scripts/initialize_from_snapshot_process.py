# Import Python modules
import json
import numpy
from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_database import RomDatabase

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return InitializeFromSnapshotProcess(model, settings["Parameters"])

class InitializeFromSnapshotProcess(KratosMultiphysics.Process):
    """A process to save the generalized coordinates (q) of the Reduced-Order Model (ROM) to disk."""

    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Validate input settings against defaults
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # # Get the rom parameters in order to pass them to the database initialization
        # rom_parameters_path = Path(settings["rom_data_folder"].GetString())+Path(settings["rom_parameters_file"].GetString())
        # with open(rom_parameters_path,'r') as parameter_file:
        #     rom_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # # Initialize the database
        # self.data_base = RomDatabase(rom_parameters) # Is the mu_names argument essential?

        # Get the model part from which the snapshots are to be retrieved
        if not settings["model_part_name"].GetString():
            raise Exception("\'model_part_name\' not provided. Please specify the model part to get the snapshots from.")
        self.model_part = model[settings["model_part_name"].GetString()]

        # # Get the mu values for the desired initial snapshot
        # self.initial_mu = settings["initial_mu"].GetVector()

        # initialize variables to be configured later
        # self.model_part = None
        # self.snapshot_variables_list = None
        # self.move_mesh = None

        # Get the values for the desired initial snapshot
        self.init_snapshot = numpy.array(settings["initial_snapshot"].GetVector())

        # Get list of variables that for the snapshot
        nodal_unknowns_list = settings["nodal_unknowns"].GetStringArray()
        if len(nodal_unknowns_list) == 0:
            err_msg = "The snapshots matrix variables need to be specified by the user in the \'nodal_unknowns\' string array."
            raise Exception(err_msg)
        if any(nodal_unknowns_list.count(var_name) > 1 for var_name in nodal_unknowns_list):
            err_msg = "There are repeated variables in the \'nodal_unknowns\' string array."
            raise Exception(err_msg)
        nodal_unknowns_list.sort()

        self.snapshot_variables_list = []
        for var_name in nodal_unknowns_list:
            if not KratosMultiphysics.KratosGlobals.HasVariable(var_name):
                err_msg = "\'{}\' variable in \'nodal_unknowns\' is not in KratosGlobals. Please check provided value.".format(var_name)
            if not KratosMultiphysics.KratosGlobals.GetVariableType(var_name):
                err_msg = "\'{}\' variable in \'nodal_unknowns\' is not double type. Please check provide double type variables (e.g. [\"DISPLACEMENT_X\",\"DISPLACEMENT_Y\"]).".format(var_name)
            self.snapshot_variables_list.append(KratosMultiphysics.KratosGlobals.GetVariable(var_name))

        # Set whether the mesh should be moved
        self.move_mesh = False
        if "DISPLACEMENT_X" in nodal_unknowns_list:
            self.move_mesh = True

    @classmethod
    def GetDefaultParameters(self):
        default_settings = KratosMultiphysics.Parameters("""{
        "help": "This process allows to retrieve a specific snapshot saved in the database and use it as initializatioin point for the new simulation",
        "initial_snapshot" : [],
        "nodal_unknowns": [],
        "model_part_name": ""
        }""")

        return default_settings
    
    # def ConfigureProcess(self, model_part, nodal_unknowns_list, move_mesh_flag):
    #     # Get the model part from which the snapshots are to be retrieved
    #     self.model_part = model_part

    #     # Set whether the mesh should be moved
    #     self.move_mesh = move_mesh_flag

    
    def ExecuteBeforeSolutionLoop(self):

        # We separate the full snapshot by the different variables included in it. Then update each one individually
        num_of_vars = len(self.snapshot_variables_list)
        nodes_array = self.model_part.Nodes
        for i, snapshot_var in enumerate(self.snapshot_variables_list):
            var_vector = self.init_snapshot[i::num_of_vars].copy()
            KratosMultiphysics.VariableUtils().SetSolutionStepValuesVector(nodes_array, snapshot_var, var_vector, 0)
            # If the simulation requires moving the mesh based on displacement, we need to keep the displacements in all
            # directions
            if self.move_mesh:
                displacement_z = None
                if snapshot_var == KratosMultiphysics.DISPLACEMENT_X:
                    displacement_x = var_vector.copy().reshape(-1,1)
                elif snapshot_var == KratosMultiphysics.DISPLACEMENT_Y:
                    displacement_y = var_vector.copy().reshape(-1,1)
                elif snapshot_var == KratosMultiphysics.DISPLACEMENT_Z:
                    displacement_z = var_vector.copy().reshape(-1,1)

        if self.move_mesh:
            if displacement_z is None:
                dim = 3
                displacements_vec = numpy.hstack([displacement_x, displacement_y]).reshape(-1)
            else:
                dim = 2
                displacements_vec = numpy.hstack([displacement_x, displacement_y, displacement_z]).reshape(-1)

            original_coords_vec = KratosMultiphysics.VariableUtils().GetInitialPositionsVector(nodes_array, dim)
            new_coords_vec=original_coords_vec+displacements_vec
            KratosMultiphysics.VariableUtils().SetCurrentPositionsVector(nodes_array, new_coords_vec)
