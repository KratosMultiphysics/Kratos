import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

## Distance import utility class definition
class DistanceImportUtility:
    def __init__(self, model_part, settings):
        self.settings = settings
        self.model_part = model_part

    def ImportDistance(self):
        import_mode = self.settings["import_mode"].GetString()

        if(import_mode == "from_GiD_file"):
            distance_file_name = self.settings["distance_file_name"].GetString()

            distance_file = open(distance_file_name, "r")
            distance_reading = False

            for line in distance_file:
                # Check if the distance reading is already finished
                if ("End Values" in line):
                    distance_reading = False

                # Read nodal distance
                if (distance_reading == True):
                    node_id = int(line[0:line.index(" ")])              # Get the nodal id. as integer
                    distance_value = float(line[line.index(" ")+1:])    # Get the distance value as float

                    if (node_id in self.model_part.Nodes):
                        # Note that the distance sign is swapped (Kratos distance criterion is opposite to GiD one)
                        self.model_part.Nodes[node_id].SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -distance_value)

                    # print("Node: ",node_id," distance value: ",self.model_part.Nodes[node_id].GetSolutionStepValue(KratosMultiphysics.DISTANCE))

                # Check if the distance reading has finished in the previous line
                if("Values" in line):
                    distance_reading = True

            # Recall to close the distance file to free memory
            distance_file.close()
