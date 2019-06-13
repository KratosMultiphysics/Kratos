from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_io import CoSimulationBaseIO

def Create(model, custom_settings):
    return CoSimulationDummyIO(model, custom_settings)

class CoSimulationDummyIO(CoSimulationBaseIO):
    def ImportCouplingInterfaceData(self, data_object, from_solver=None):
        pass

    def ImportCouplingInterface(self, mesh_config, from_solver=None):
        pass

    def ExportCouplingInterfaceData(self, data_object, to_solver=None):
        pass

    def ExportCouplingInterface(self, mesh_config, to_solver=None):
        pass

    def PrintInfo(self):
        print("This is the dummy-IO")

    def Check(self):
        pass
