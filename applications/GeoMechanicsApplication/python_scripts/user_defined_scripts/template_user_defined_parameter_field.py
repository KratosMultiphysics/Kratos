
from KratosMultiphysics.GeoMechanicsApplication.user_defined_scripts.user_defined_parameter_field_base \
    import ParameterFieldBase


class ParameterField(ParameterFieldBase):
    """
    Base class of a user defined parameter field
    """

    def __init__(self):
        super(ParameterField, self).__init__()

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
        super(ParameterField, self).validate_input(input_dict, output_dict)

        # add custom validation here

    def validate_output(self):
        """
        Validates input for custom parameter field

        Parameters
        ----------
        input_dict dictionary with input values
        output_dict dictionary with output values

        Returns
        -------

        """
        super(ParameterField, self).validate_output()

        # add custom validation here

    def generate_field(self):
        """
        Generates custom parameter field

        """
        super(ParameterField, self).run()

        input_dict = self.get_input()
        output_dict = self.get_output()

        # add custom run functionalities here