from abc import ABC
from abc import abstractclassmethod

import KratosMultiphysics as Kratos


class MeasurementDataGenerator(ABC):

    """_summary_
    Interface for classes that are capable of writing a measurement.json-file based on the MeasurementFileStructure.py Dataclasses.
    """

    def __init__():
        pass

    @abstractclassmethod
    def write_measurement_data_file(self, simulation_parameters: Kratos.Parameters, mdpa_model_file_name:str ,output_file_name: str):
        pass
