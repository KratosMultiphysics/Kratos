
class ParameterFieldBase:
    """
    Base class of a user defined parameter field
    """

    def __init__(self):
        self.__input_dict = {}
        self.__output_dict = {}

    def get_input(self):
        return self.__input_dict

    def get_output(self):
        return self.__output_dict

    def validate_input(self, input_dict, output_dict):
        """
        Validates input for custom parameter field

        Parameters
        ----------
        input_dict dictionary with input values
        output_dict dictionary with output values

        Returns
        -------

        """

        self.__input_dict = input_dict
        self.__output_dict = output_dict

        if "values" not in self.__input_dict:
            raise KeyError("'values' is not a key in the input dictionary for the parameter field")

        if "coordinates" not in self.__input_dict:
            raise KeyError("'coordinates' is not a key in the input dictionary for the parameter field")

        if len(self.__input_dict["values"]) != len(self.__input_dict["coordinates"]):
            raise AssertionError("'values' and 'coordinates' do not have the same length in the input dictionary for "
                                 "the parameter field")

    def validate_output(self):
        """
        Validates output for custom parameter field


        """

        if "values" not in self.__output_dict:
            raise KeyError("'values' is not a key in the output dictionary for the parameter field")

        if len(self.__output_dict["values"]) != len(self.__input_dict["coordinates"]):
            raise AssertionError("'values' and 'coordinates' do not have the same length in the output dictionary for "
                                 "the parameter field")

    def generate_field(self):
        """
        Generates custom parameter field

        """
        pass
