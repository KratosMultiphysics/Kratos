from KratosMultiphysics import Parameters, THICKNESS
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

class ShellThicknessAdaptor():
    """ This simple process is designed to modify the thickness of a shell model
        after each shape change, in order to keep the mass constant.
        It assumes that the model one property and this property has THICKNESS
    """

    def __init__(self, model, settings):
        self.model_part = model[settings["model_part_name"].GetString()]

        dummy_parameters = Parameters("""{
            "response_type"          : "mass",
            "gradient_mode"          : "finite_differencing",
            "step_size"              : 1e-6
        }""")
        self.response_function_utility = StructuralMechanicsApplication.MassResponseFunctionUtility(self.model_part, dummy_parameters)
        self.initial_mass = None
        self.initial_thickness_dict = {}
        self.current_thickness_dict = {}

    def Initialize(self):
        """ This function relies on already read properties by the analyzer"""

        self.response_function_utility.Initialize()
        mass = self.response_function_utility.CalculateValue()
        self.initial_mass = mass

        for prop in self.model_part.GetProperties():
            print(prop)
            if not prop.Has(THICKNESS):
                print ("> WARNING: ShellThicknessAdaptor property {} does not have THICKNESS!".format(prop.Id))
                continue

            thickness = prop.GetValue(THICKNESS)
            self.initial_thickness_dict[prop.Id] = thickness
            self.current_thickness_dict[prop.Id] = thickness

        if len(self.initial_thickness_dict) == 0:
            raise RuntimeError("ShellThicknessAdaptor: The model part does not have any THICKNESS property!")


    def AfterMeshUpdate(self):

        mass = self.response_function_utility.CalculateValue()

        print("\n> ----------------------------------------")
        print("> Shell thickness adaptor:")
        print(">    initial mass:     ", self.initial_mass)
        print(">    current mass:     ", mass)

        for prop_id, initial_thickness in self.initial_thickness_dict.items():

            current_thickness = self.current_thickness_dict[prop_id]
            new_thickness = current_thickness * self.initial_mass/mass
            self.current_thickness_dict[prop_id] = new_thickness
            prop = self.model_part.GetProperties(prop_id, 0)
            prop.SetValue(THICKNESS, new_thickness)

            print("> -----------------------------------------")
            print(">  Property {}".format(prop_id))
            print(">    initial thickness:", initial_thickness)
            print(">    current thickness:", current_thickness)
            print(">    new thickness:    ", new_thickness)

        print("> -----------------------------------------")
        print(">    new mass:          ", self.response_function_utility.CalculateValue())
        print("> -----------------------------------------\n")
