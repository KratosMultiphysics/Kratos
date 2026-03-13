# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.kratos_utilities

VTK_QUAD = 9
VERSION = (2, 4)
GLOBAL_TYPE = 'PartitionedDataSetCollection'
SUB_TYPE = "UnstructuredGrid"

try:
    import h5py as h5
except ImportError as e:
    raise ImportError(
        "h5py is required for this functionality but it is not installed. "
        "Please install it with: pip install h5py"
    ) from e

def Factory(settings, model):
    if not (isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return IgaVTKOutputProcess(model, settings["Parameters"])

class IgaVTKOutputProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "output_file_name"          : "",
            "brep_surface_ids"          : [],                              
            "model_part_name"           : "",
            "file_label"                : "step",
            "output_control_type"       : "step",
            "output_frequency"          : 1.0
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        ## Get the model part
        self.model_part = model[self.params["model_part_name"].GetString()]

        ## Get the geometry
        self.brep_surface_ids = [
            self.params["brep_surface_ids"][i].GetInt()
            for i in range(self.params["brep_surface_ids"].size())
        ]

        self.output_file_name = self.params["output_file_name"].GetString()
        if not self.output_file_name.endswith(".vtkhdf"):
            self.output_file_name += ".vtkhdf"
        
        # self.nodal_results_scalar, self.nodal_results_vector = \
        #     CreateVariablesListFromInput(self.params["nodal_results"])

        # self.integration_point_results_scalar, self.integration_point_results_vector = \
        #     CreateVariablesListFromInput(self.params["integration_point_results"])

        # Set up output frequency and format
        output_file_label = self.params["file_label"].GetString()
        if output_file_label == "time":
            self.output_label_is_time = True
        elif output_file_label == "step":
            self.output_label_is_time = False
        else:
            msg = '{} Error: Unknown value "{}" read for parameter "file_label"'.format(self.__class__.__name__,output_file_label)
            raise Exception(msg)

        output_control_type = self.params["output_control_type"].GetString()
        if output_control_type == "time":
            self.output_control_is_time = True
        elif output_control_type == "step":
            self.output_control_is_time = False
        else:
            err_msg  = 'The requested "output_control_type" "' + output_control_type
            err_msg += '" is not available!\nAvailable options are: "time", "step"'
            raise Exception(err_msg)

        self.output_frequency = self.params["output_frequency"].GetDouble()

        self.step_count = 0
        self.printed_step_count = 0
        self.next_output = 0.0

    def ExecuteBeforeSolutionLoop(self):
        with open(self.output_file_name, 'w') as output_file:
            output_file.write("Post Results File 1.0\n")

    def PrintOutput(self):
        time = GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        self.printed_step_count += 1
        self.model_part.ProcessInfo[KratosMultiphysics.PRINTED_STEP] = self.printed_step_count
        if self.output_label_is_time:
            label = time
        else:
            label = self.printed_step_count

        #TODO: write vtkhdf output
        file = h5.File(self.output_file_name, "w")
        root = file.create_group("VTKHDF", track_order=True)
        root.attrs["Version"] = VERSION
        root.attrs["Type"] = GLOBAL_TYPE

        assembly = root.create_group("Assembly", track_order=True)

        # Output subgroup for each brep
        index = 0
        for brep_id in self.brep_surface_ids:
            brep_surface = self.model_part.GetGeometry(brep_id)
            group_name = f"brep_surface_{brep_id}"
            surface_group = root.create_group(group_name)
            surface_group.attrs["Version"] = VERSION
            surface_group.attrs["Type"] = SUB_TYPE
            self.__write_out_surface(surface_group, brep_surface)

            # ---- Assembly ----
            assem_group = assembly.create_group(f"brep_{brep_id}")
            #assem_group.attrs['Index'] = index
            assem_group[group_name] = h5.SoftLink(f"/VTKHDF/{group_name}")
            assem_group[group_name].attrs['Index'] = index

            index += 1

        # Schedule next output
        if self.output_frequency > 0.0: # Note: if == 0, we'll just always print
            if self.output_control_is_time:
                while GetPrettyTime(self.next_output) <= time:
                    self.next_output += self.output_frequency
            else:
                while self.next_output <= self.step_count:
                    self.next_output += self.output_frequency

    def IsOutputStep(self):
        if self.output_control_is_time:
            time = GetPrettyTime(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
            return (time >= GetPrettyTime(self.next_output))
        else:
            return ( self.step_count >= self.next_output )

    def __write_out_surface(self, surface_group:h5.Group, brep_surface:KM.BrepSurface):
        knots_u = brep_surface.KnotsU()
        knots_v = brep_surface.KnotsV()

        num_u = len(knots_u)
        num_v = len(knots_v)
        num_cells = (num_u-1) * (num_v-1)
        num_points = num_u * num_v

        # get the Points
        grid_points = KM.Matrix(num_points, 3)
        for j, v in enumerate(knots_v):
            for i, u in enumerate(knots_u):
                X = KM.Array3()
                X[0] = 0; X[1] = 0; X[2] = 0
                local_coord = KM.Array3()
                local_coord[0] = u; local_coord[1] = v; local_coord[2] = 0.0
                X = brep_surface.GlobalCoordinates(local_coord)
                print(X)

                # Column-major indexing
                idx = i + j * len(knots_u)
                grid_points[idx, 0] = local_coord[0]
                grid_points[idx, 1] = local_coord[1]
                grid_points[idx, 2] = local_coord[2]
        
        # get the Types
        cell_types = [VTK_QUAD] * num_cells

        # get the Connectivity & Offsets
        connectivity = KM.Vector(num_cells * 4)  # 4 nodes per quad
        offsets = KM.Vector(num_cells+1)
        cell_idx = 0
        conn_idx = 0
        for j in range(num_v - 1):
            for i in range(num_u - 1):
                # Column-major indices
                n0 = i     + j * num_u
                n1 = (i+1) + j * num_u
                n2 = (i+1) + (j+1) * num_u
                n3 = i     + (j+1) * num_u

                # Fill connectivity
                connectivity[conn_idx    ] = n0
                connectivity[conn_idx + 1] = n1
                connectivity[conn_idx + 2] = n2
                connectivity[conn_idx + 3] = n3

                offsets[cell_idx] = conn_idx  
                conn_idx += 4
                cell_idx += 1

        offsets[cell_idx] = conn_idx

        # write to hdf5 group
        surface_group.create_dataset("Points", data=grid_points, dtype="f")
        surface_group.create_dataset("Connectivity", data=connectivity, dtype="i8")
        surface_group.create_dataset("Offsets", data=offsets, dtype="i8")
        surface_group.create_dataset("Types", data=cell_types, dtype="uint8")
        surface_group.create_dataset("NumberOfPoints", data=(num_points,), dtype="i8")
        surface_group.create_dataset("NumberOfCells", data=(num_cells,), dtype="i8")
        surface_group.create_dataset("NumberOfConnectivityIds", data=(len(connectivity),), dtype="i8")

def GetPrettyTime(time):
    pretty_time = "{0:.12g}".format(time)
    pretty_time = float(pretty_time)
    return pretty_time

def CreateVariablesListFromInput(param):
    '''Parse a list of variables from input.'''
    scalar_variables = []
    vector_variables = []
    admissible_scalar_types = ["Bool", "Integer", "Unsigned Integer", "Double"]
    admissible_vector_types = ["Array", "Vector"]

    variable_list = KratosMultiphysics.kratos_utilities.GenerateVariableListFromInput(param)

    # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
    for variable in variable_list:
        if KratosMultiphysics.KratosGlobals.GetVariableType(variable.Name()) in admissible_scalar_types:
            scalar_variables.append(variable)
        elif KratosMultiphysics.KratosGlobals.GetVariableType(variable.Name()) in admissible_vector_types:
            vector_variables.append(variable)
        else:
            raise Exception("unsupported variables type: " + str(type(variable)))

    return scalar_variables, vector_variables
