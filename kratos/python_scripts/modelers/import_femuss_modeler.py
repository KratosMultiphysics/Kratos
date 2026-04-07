import KratosMultiphysics

class ImportFemussModeler(KratosMultiphysics.Modeler):

    def __init__(self, model, settings):
        super().__init__(model, settings)

        # Cannot validate as settings may differ among input types
        settings.AddMissingParameters(self.__GetDefaultSettings())

        # Declare required member variables
        self.model_part = None
        self.settings = settings

        # Create the import destination model part
        # It is mandatory to do this when the modeler is instantiated to have the model part created before the solvers add the variables
        model_part_name = self.settings["model_part_name"].GetString()
        if not model_part_name:
            err_msg = "Missing 'model_part_name' in input settings. This is where the imported model part is to be stored."
            raise Exception(err_msg)
        else:
            self.model_part = model.GetModelPart(model_part_name)

        # Save the problem dimension for later use
        if not self.model_part.ProcessInfo.Has(KratosMultiphysics.DOMAIN_SIZE):
            raise RuntimeError(f"Missing 'DOMAIN_SIZE' in model part ProcessInfo container of model part '{self.model_part.Name}'.")
        self.dimension = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.dimension not in [2, 3]:
            raise ValueError(f"Unsupported problem dimension: {self.dimension}. Only 2D and 3D are supported.")

    def SetupGeometryModel(self):
        super().SetupGeometryModel()
        self.__ReadFemussGeometry()
        self.__CreateFixitySubModelParts()

    def PrepareGeometryModel(self):
        super().PrepareGeometryModel()

    def SetupModelPart(self):
        super().SetupModelPart()

    @classmethod
    def __GetDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters('''{
            "echo_level" : 0,
            "geo_input_filename" : "",
            "fix_input_filename" : "",
            "model_part_name" : ""
        }''')
        return default_settings

    def __ReadFemussGeometry(self):
        # Read geometry from file
        element_ids, connectivity, node_ids, coordinates = self.__ParseGeoFile(self.settings["geo_input_filename"].GetString())

        # Consistency checks
        if len(node_ids) != len(coordinates):
            raise RuntimeError("Mismatch: node_ids vs coordinates")

        if len(element_ids) != len(connectivity):
            raise RuntimeError("Mismatch: element_ids vs connectivity")

        # Create nodes
        for i, node_id in enumerate(node_ids):
            x = coordinates[i][0]
            y = coordinates[i][1]
            z = coordinates[i][2] if len(coordinates[i]) == 3 else 0.0
            self.model_part.CreateNewNode(int(node_id), x, y, z)

        # Create geometries
        print(connectivity)
        geometry_name = self.__GetGeometryName(len(connectivity[0]))
        for i, geom_id in enumerate(element_ids):
            self.model_part.CreateNewGeometry(
                geometry_name,
                int(geom_id),
                [int(n) for n in connectivity[i]])

    def __GetGeometryName(self, num_nodes):
        if self.dimension == 2:
            if num_nodes == 3:
                return "Triangle2D3"
            elif num_nodes == 4:
                return "Quadrilateral2D4"
        elif self.dimension == 3:
            if num_nodes == 4:
                return "Tetrahedra3D4"
            elif num_nodes == 8:
                return "Hexahedra3D8"

        raise RuntimeError(f"Unsupported element with {num_nodes} nodes in {self.dimension}D")

    def __CreateFixitySubModelParts(self):
        groups = self.__ParseFixFile()
        for name, data in groups.items():
            sub_model_part = self.model_part.CreateSubModelPart(name)
            sub_model_part.AddNodes(data["nodes"])

    def __BuildAutoGroupName(self, fix, val):
        fix_str = ''.join('1' if f else '0' for f in fix)
        val_str = '_'.join(f"{v:.6g}" for v in val)
        return f"Fix_{fix_str}_{val_str}"

    # ---------------------------------------------------------------------
    # --- FEMUSS parser (no numpy, C++-ready logic) ------------------------
    # ---------------------------------------------------------------------

    def __ParseGeoFile(self, filepath):
        element_ids = []
        connectivity = []

        node_ids = []
        coordinates = []

        new_format = False
        in_elements = False
        in_coordinates = False

        with open(filepath, 'r') as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue

                # ---------------------------------
                # BLOCK DETECTION
                # ---------------------------------
                if line.startswith("ELEMENTS"):
                    in_elements = True
                    new_format = "NEWFORMAT" in line
                    continue

                elif line == "END_ELEMENTS":
                    in_elements = False
                    continue

                elif line == "COORDINATES":
                    in_coordinates = True
                    continue

                elif line == "END_COORDINATES":
                    in_coordinates = False
                    continue

                # ---------------------------------
                # ELEMENTS
                # ---------------------------------
                if in_elements:
                    parts = line.split()

                    elem_id = int(parts[0])

                    if new_format:
                        n_nodes = int(parts[1])
                        conn = [int(p) for p in parts[2:2 + n_nodes]]
                    else:
                        conn = [int(p) for p in parts[1:]]

                    element_ids.append(elem_id)
                    connectivity.append(conn)

                # ---------------------------------
                # COORDINATES
                # ---------------------------------
                elif in_coordinates:
                    parts = line.split()

                    node_ids.append(int(parts[0]))
                    coordinates.append([float(p) for p in parts[1:]])

        return element_ids, connectivity, node_ids, coordinates

    def __ParseFixFile(self):
        groups = {}
        in_nodes = False
        current_group = None

        with open(self.settings["fix_input_filename"].GetString(), 'r') as f:
            for raw_line in f:
                line = raw_line.strip()

                if not line:
                    continue

                if line == "ON_NODES":
                    in_nodes = True
                    continue

                if line == "END_ON_NODES":
                    break

                if not in_nodes:
                    continue

                # -----------------------------------------
                # GROUP START / END
                # -----------------------------------------
                if line.startswith("$"):
                    name = line[1:].strip()

                    if name != "": # START group
                        current_group = name
                        if current_group not in groups:
                            groups[current_group] = []
                    else: # END group
                        current_group = None

                    continue

                # -----------------------------------------
                # NODE LINE
                # -----------------------------------------
                parts = line.split()
                node_id = int(parts[0])
                fix = tuple(c == '1' for c in parts[1][:self.dimension])
                val = tuple(round(float(parts[2+i]), 12) for i in range(self.dimension))

                if not any(fix):
                    continue

                entry = {
                    "node_id": node_id,
                    "fixity": fix,
                    "value": val
                }

                # -----------------------------------------
                # CASE A: inside named group
                # -----------------------------------------
                if current_group is not None:
                    groups[current_group].append(entry)

                # -----------------------------------------
                # CASE B: outside groups → auto grouping
                # -----------------------------------------
                else:
                    name = self.__BuildAutoGroupName(fix, val)
                    if name not in groups:
                        groups[name] = []
                    groups[name].append(entry)

        return self.__SplitGroups(groups)

    def __SplitGroups(self, groups):
        split = {}

        for gname, entries in groups.items():
            subgroups = {}

            for entry in entries:
                key = (entry["fixity"], entry["value"])

                if key not in subgroups:
                    subgroups[key] = []

                subgroups[key].append(entry["node_id"])

            # Naming
            for i, ((fix, val), nodes) in enumerate(subgroups.items()):
                name = gname if len(subgroups) == 1 else f"{gname}_{i}"

                split[name] = {
                    "nodes": nodes,
                    "fixity": fix,
                    "value": val
                }

        return split

def Factory(model, settings):
    return ImportFemussModeler(model, settings)
