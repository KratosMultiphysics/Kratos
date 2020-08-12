# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO
from pathlib import Path
import os

def Create(model, settings, solver_name):
    return PingPongIO(model, settings, solver_name)

class PingPongIO(CoSimulationIO):
    """This is the IO wrapper for the PING-PONG example.
    """
    def ImportCouplingInterface(self, interface_config):
        pass

    def ExportCouplingInterface(self, interface_config):
        pass

    def ImportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            data_name = self.solver_name+"_"+interface_data.name
            data_file_name = "EMPIRE_datafield_" + data_name + ".dat"
            if( self.__FileExists(data_file_name) ):
                with open(data_file_name, 'r') as input_file:
                    length_of_data = int(input_file.readline() )
                    if(length_of_data != 1):
                        AttributeError("There is more than expected data in the file. PingPongIO takes only one scalar")
                    data = float(input_file.readline())
                    node = interface_data.GetModelPart().Nodes[1]
                    node.SetSolutionStepValue(interface_data.variable, data)
                    input_file.close()
                try:
                    cwd = os.getcwd()
                    file_to_rmv = os.path.join(cwd, data_file_name)
                    print("Removing file :: ", file_to_rmv)
                    os.remove(data_file_name)
                except OSError as e: # name the Exception `e`
                    print ("Failed with:", e.strerror) # look what it says
                    print ("Error code:", e.code)
        else:
            raise NotImplementedError('Importing interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            data_name = self.solver_name+"_"+interface_data.name
            data_file_name = "EMPIRE_datafield_" + data_name + ".dat"
            with open(data_file_name, 'w') as output_file:
                node = interface_data.GetModelPart().Nodes[1]
                data = node.GetSolutionStepValue(interface_data.variable)
                output_file.write(str(1))
                output_file.write("\n")
                output_file.write(str(data))
                output_file.close()
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def PrintInfo(self):
        print("This is the PingPong-IO")

    def Check(self):
        pass

    def __FileExists(self, file_name_including_path):
        given_file = Path(file_name_including_path)
        return given_file.is_file()

