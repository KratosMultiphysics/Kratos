
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

        default_umat_parameters = [0.3, 0.0, 30.0, 0.0, 0.0, 1.0, 0.0]

        new_values = []
        for value, coord in zip(input_dict["values"], input_dict["coordinates"]):
            new_value = value[0] * 2 * coord[0] + value[0] * 3 * coord[1]
            new_values.append([new_value * 10] + default_umat_parameters)

        output_dict["values"] = new_values



