
from KratosMultiphysics.GeoMechanicsApplication.user_defined_scripts.user_defined_parameter_field_base \
    import ParameterFieldBase


class ParameterField(ParameterFieldBase):
    """
    Base class of a user defined parameter field
    """

    def generate_field(self):
        """
        Creates custom parameter field

        """
        super(ParameterField, self).generate_field()

        input_dict = self.get_input()
        output_dict = self.get_output()

        # add custom run functionalities here

        new_values = []
        for value, coord in zip(input_dict["values"], input_dict["coordinates"]):
            new_value = value * 2 * coord[0] + value * 3 * coord[1]
            new_values.append(new_value)

        output_dict["values"] = new_values



