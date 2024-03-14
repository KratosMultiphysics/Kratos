# Import Python modules
import json
import numpy
from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return CalculateRomBasisOutputProcess(model, settings["Parameters"])

class CalculateRomBasisOutputProcess(KratosMultiphysics.OutputProcess):
    """A process to set the snapshots matrix and calculate the ROM basis from it."""

    def __init__(self, model, settings):
        KratosMultiphysics.OutputProcess.__init__(self)

        # Validate input settings against defaults
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # Get the model part from which the snapshots are to be retrieved
        if not settings["model_part_name"].GetString():
            raise Exception("\'model_part_name\' not provided. Please specify the model part to get the snapshots from.")
        self.model_part = model[settings["model_part_name"].GetString()]

        # Set the snapshots output control and interval
        snapshots_control_type = settings["snapshots_control_type"].GetString()
        if snapshots_control_type == "time":
            self.snapshots_control_is_time = True
        elif snapshots_control_type == "step":
            self.snapshots_control_is_time = False
        else:
            err_msg = "Unknown value \'{}\' for \'snapshots_control_type\'. Available options are \'time\' and \'step\'.".format(snapshots_control_type)
            raise Exception(err_msg)
        self.snapshots_interval = settings["snapshots_interval"].GetDouble()

        # Get the variables list to be used to get the snapshots matrix information
        # Note that we sort the snapshot variables list alphabetically
        # This is required in order to establish a consensum for the possible visualization model part projections
        nodal_unknowns = settings["nodal_unknowns"].GetStringArray()
        if len(nodal_unknowns) == 0:
            err_msg = "The snapshots matrix variables need to be specified by the user in the \'nodal_unknowns\' string array."
            raise Exception(err_msg)
        if any(nodal_unknowns.count(var_name) > 1 for var_name in nodal_unknowns):
            err_msg = "There are repeated variables in the \'nodal_unknowns\' string array."
            raise Exception(err_msg)
        nodal_unknowns.sort()

        self.snapshot_variables_list = []
        for var_name in nodal_unknowns:
            if not KratosMultiphysics.KratosGlobals.HasVariable(var_name):
                err_msg = "\'{}\' variable in \'nodal_unknowns\' is not in KratosGlobals. Please check provided value.".format(var_name)
            if not KratosMultiphysics.KratosGlobals.GetVariableType(var_name):
                err_msg = "\'{}\' variable in \'nodal_unknowns\' is not double type. Please check provide double type variables (e.g. [\"DISPLACEMENT_X\",\"DISPLACEMENT_Y\"]).".format(var_name)
            self.snapshot_variables_list.append(KratosMultiphysics.KratosGlobals.GetVariable(var_name))

        # Set the ROM basis output settings
        self.rom_basis_output_format = settings["rom_basis_output_format"].GetString()
        rom_basis_output_available_formats = ["json", "numpy"]
        if self.rom_basis_output_format not in rom_basis_output_available_formats:
            err_msg = "Provided \'rom_basis_output_format\' is {}. Available options are \'json\' and \'numpy\'.".format(self.rom_basis_output_format)
            raise Exception(err_msg)

        self.rom_basis_output_name = settings["rom_basis_output_name"].GetString()

        self.rom_basis_output_folder = Path(settings["rom_basis_output_folder"].GetString())

        # Get the SVD truncation tolerance
        self.svd_truncation_tolerance = settings["svd_truncation_tolerance"].GetDouble()

        # Initialize output interval data
        self.next_output = 0.0

        # Initialize the snapshots data list
        self.snapshots_data_list = []

        # Set the flag allowing to run multiple simulations using this process #TODO cope with arbitrarily large cases (parallelism)
        self.rom_manager = settings["rom_manager"].GetBool()


    @classmethod
    def GetDefaultParameters(self):
        default_settings = KratosMultiphysics.Parameters("""{
            "help": "A process to set the snapshots matrix and calculate the ROM basis from it.",
            "model_part_name": "",
            "rom_manager" : false,
            "snapshots_control_type": "step",
            "snapshots_interval": 1.0,
            "nodal_unknowns": [],
            "rom_basis_output_format": "numpy",
            "rom_basis_output_name": "RomParameters",
            "rom_basis_output_folder" : "rom_data",
            "svd_truncation_tolerance": 1.0e-6
        }""")

        return default_settings

    def IsOutputStep(self):
        if self.snapshots_control_is_time:
            time = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return time >= self.__GetPrettyFloat(self.next_output)
        else:
            step = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.STEP])
            return step >= self.next_output

    def PrintOutput(self):
        # Save the data in the snapshots data list
        aux_data_array = []
        for node in self.model_part.Nodes:
            for snapshot_var in self.snapshot_variables_list:
                aux_data_array.append(node.GetSolutionStepValue(snapshot_var))
        self.snapshots_data_list.append(aux_data_array)

        # Schedule next snapshot output
        if self.snapshots_interval > 0.0: # Note: if == 0, we'll just always print
            if self.snapshots_control_is_time:
                time = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
                while self.__GetPrettyFloat(self.next_output) <= time:
                    self.next_output += self.snapshots_interval
            else:
                step = self.__GetPrettyFloat(self.model_part.ProcessInfo[KratosMultiphysics.STEP])
                while self.next_output <= step:
                    self.next_output += self.snapshots_interval


    def _GetSnapshotsMatrix(self):
        snapshots_matrix = numpy.empty((self.n_nodal_unknowns*self.n_nodes,self.n_data_cols))
        for i_col in range(self.n_data_cols):
            aux_col = numpy.array(self.snapshots_data_list[i_col])
            snapshots_matrix[:,i_col] = aux_col.transpose()
        return snapshots_matrix


    def _PrintRomBasis(self, snapshots_matrix):
        # Initialize the Python dictionary with the default settings
        # Note that this order is kept if Python 3.6 onwards is used
        rom_basis_dict = {
            "rom_manager" : False,
            "train_hrom": False,
            "run_hrom": False,
            "projection_strategy": "galerkin",
            "assembling_strategy": "global",
            "rom_format": "numpy",
            "rom_settings": {
                "rom_bns_settings": {}
            },
            "hrom_settings": {},
            "nodal_modes": {},
            "elements_and_weights" : {}
        }
        #TODO: I'd rename elements_and_weights to hrom_weights

        if self.rom_manager:
            rom_basis_dict["rom_manager"] = True
        rom_basis_dict["hrom_settings"]["hrom_format"] = self.rom_basis_output_format
        n_nodal_unknowns = len(self.snapshot_variables_list)

        # Calculate the randomized SVD of the snapshots matrix
        u,_,_,_= RandomizedSingularValueDecomposition().Calculate(snapshots_matrix, self.svd_truncation_tolerance)

        # Save the nodal basis
        rom_basis_dict["rom_settings"]["nodal_unknowns"] = [var.Name() for var in self.snapshot_variables_list]
        rom_basis_dict["rom_settings"]["number_of_rom_dofs"] = numpy.shape(u)[1] #TODO: This is way misleading. I'd call it number_of_basis_modes or number_of_rom_modes
        rom_basis_dict["projection_strategy"] = "galerkin" # Galerkin: (Phi.T@K@Phi dq= Phi.T@b), LSPG = (K@Phi dq= b), Petrov-Galerkin = (Psi.T@K@Phi dq = Psi.T@b)
        rom_basis_dict["assembling_strategy"] = "global" # Assemble the ROM globally or element by element: "global" (Phi_g @ J_g @ Phi_g), "element by element" sum(Phi_e^T @ K_e @ Phi_e)
        rom_basis_dict["rom_format"] = self.rom_basis_output_format
        rom_basis_dict["rom_settings"]["petrov_galerkin_number_of_rom_dofs"] = 0
        #NOTE "petrov_galerkin_number_of_rom_dofs" is not used unless a Petrov-Galerkin simulation is called, in which case it shall be modified either manually or from the RomManager

        # Create the folder if it doesn't already exist
        if not self.rom_basis_output_folder.exists():
            self.rom_basis_output_folder.mkdir(parents=True)

        if self.rom_basis_output_format == "json":
            # Storing modes in JSON format
            i = 0
            for node in self.model_part.Nodes:
                rom_basis_dict["nodal_modes"][node.Id] = u[i:i+n_nodal_unknowns].tolist()
                i += n_nodal_unknowns

        elif self.rom_basis_output_format == "numpy":
            # Storing modes in Numpy format
            node_ids = []
            for node in self.model_part.Nodes:
                node_ids.append(node.Id)
            node_ids = numpy.array(node_ids)
            numpy.save(self.rom_basis_output_folder / "RightBasisMatrix.npy", u)
            numpy.save(self.rom_basis_output_folder / "NodeIds.npy", node_ids)
        else:
            err_msg = "Unsupported output format {}.".format(self.rom_basis_output_format)
            raise Exception(err_msg)

        # Creating the ROM JSON file containing or not the modes depending on "self.rom_basis_output_format"
        output_filename = self.rom_basis_output_folder / f"{self.rom_basis_output_name}.json"
        with output_filename.open('w') as f:
            json.dump(rom_basis_dict, f, indent = 4)



    def ExecuteFinalize(self):
        # Prepare a NumPy array with the snapshots data
        self.n_nodes = self.model_part.NumberOfNodes()
        self.n_data_cols = len(self.snapshots_data_list)
        self.n_nodal_unknowns = len(self.snapshot_variables_list)

        if not self.rom_manager:
            self._PrintRomBasis(self._GetSnapshotsMatrix())

    def __GetPrettyFloat(self, number):
        float_format = "{:.12f}"
        pretty_number = float(float_format.format(number))
        return pretty_number