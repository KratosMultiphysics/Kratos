
class DummySolverWrapper(object):
    '''dummy object used for testing to emulate the behavior of the SolverWrapper'''
    def __init__(self, interface_data_dict):
        self.interface_data_dict = interface_data_dict

    def GetInterfaceData(self, data_name):
        return self.interface_data_dict[data_name]