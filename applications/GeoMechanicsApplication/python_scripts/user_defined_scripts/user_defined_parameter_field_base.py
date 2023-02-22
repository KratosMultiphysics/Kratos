
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

        assert "values" in self.__input_dict
        assert "coordinates" in self.__input_dict

        assert len(self.__input_dict["values"]) == len(self.__input_dict["coordinates"])

    def validate_output(self):
        """
        Validates output for custom parameter field


        """

        assert "values" in self.__output_dict
        assert len(self.__output_dict["values"]) == len(self.__input_dict["coordinates"])

    def run(self):
        """
        Runs custom parameter field

        """
        pass