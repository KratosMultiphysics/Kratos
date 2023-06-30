from abc import ABC
from abc import abstractclassmethod


class MeasurementDataGenerator(ABC):

    def __init__():
        pass

    @abstractclassmethod
    def write_measurement_data_file(self, file_name: str):
        pass
