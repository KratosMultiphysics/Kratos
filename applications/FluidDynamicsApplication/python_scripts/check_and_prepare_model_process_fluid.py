import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CheckAndPrepareModelProcessFluid(Model, settings["Parameters"])

class CheckAndPrepareModelProcessFluid(KratosMultiphysics.Process):
    def __init__(self, model, parameters):
        super().__init__()

        # Check input type to keep backwards compatibility
        if isinstance(model, KratosMultiphysics.Model):
            self.model = model
            __model_part_input = False
        else:
            __model_part_input = True
            self.model = model.GetModel() # Note that in here model is actually a model part
            root_model_part_name = model.Name # note that in here model is actually a model part
            IssueDeprecationWarning('CheckAndPrepareModelProcessFluid', "Construction from model part is deprecated. Please provide a \'Model\' container and use full model part names.")

        # Validate and assign settings
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        if not parameters["volume_model_part_name"].GetString():
            raise Exception("Please define the \"volume_model_part_name\" (string) argument.")

        self.skin_name_list = parameters["skin_parts"].GetStringArray()
        self.volume_model_part_name = parameters["volume_model_part_name"].GetString()
        self.assign_neighbour_elements = parameters["assign_neighbour_elements_to_conditions"].GetBool()

        # Check and issue partial model part names deprecation warnings
        # Note that we do not allow partial model part names when using the official model-based I/O
        if len(self.volume_model_part_name.split('.')) == 1:
            msg = f'Partial model part name \'{self.volume_model_part_name}\' found in \'volume_model_part_name\'. Please provide full model part names.'
            if __model_part_input:
                IssueDeprecationWarning('CheckAndPrepareModelProcessFluid', msg)
                # We assume that the volume model part is a submodelpart of the root
                # If not so, something that we allowed in the past, we take the root as volume model part
                self.volume_model_part_name = f"{root_model_part_name}.{self.volume_model_part_name}"
                if not self.model.HasModelPart(self.volume_model_part_name):
                    self.volume_model_part_name = f"{root_model_part_name}"
            else:
                raise Exception(msg)

        for i in range(len(self.skin_name_list)):
            name = self.skin_name_list[i]
            if len(name.split('.')) == 1:
                self.__deprecated_input = True
                msg = f'Partial model part name \'{name}\' found in \'skin_parts\'. Please provide full model part names.'
                if __model_part_input:
                    IssueDeprecationWarning('CheckAndPrepareModelProcessFluid', msg)
                    self.skin_name_list[i] = f"{root_model_part_name}.{name}"
                else:
                    raise Exception(msg)

    def GetDefaultParameters(self):
        default_parameters = KratosMultiphysics.Parameters(r'''{
            "volume_model_part_name" : "",
            "skin_parts" : [],
            "assign_neighbour_elements_to_conditions" : true
        }''')
        return default_parameters

    def Execute(self):
        # Get the model parts conforming the fluid mesh (skin and volume)

        print(self.volume_model_part_name)

        volume_model_part = self.model.GetModelPart(self.volume_model_part_name)
        skin_parts = [self.model.GetModelPart(name) for name in self.skin_name_list]

        # Check that volume model part and skin model part share the same root
        root_model_part = volume_model_part.GetRootModelPart()
        for skin_part in skin_parts:
            if skin_part.GetRootModelPart().Name != root_model_part.Name:
                raise Exception(f"Volume model part \'{volume_model_part.FullName()}\' and skin model part \'{skin_part.FullName()}\' have different root model parts.")

        # Construct a model part which contains both the skin and the volume
        computational_model_part_name = "fluid_computational_model_part"
        if root_model_part.HasSubModelPart(computational_model_part_name):
            fluid_computational_model_part = root_model_part.GetSubModelPart(computational_model_part_name)
        else:
            fluid_computational_model_part = root_model_part.CreateSubModelPart("fluid_computational_model_part")
            fluid_computational_model_part.AddNodes([node.Id for node in volume_model_part.Nodes])
            fluid_computational_model_part.AddElements([element.Id for element in volume_model_part.Elements])
            cond_ids = []
            for skin_part in skin_parts:
                cond_ids.extend([condition.Id for condition in skin_part.Conditions])
            fluid_computational_model_part.AddConditions(list(set(cond_ids)))

        # Check if the mesh is made of simplex elements
        # If the mesh contains non-simplex elements, we call the parent element assign utility, which works for any type of geometries
        # If the mesh is only made by simplex, the multipurpose tetrahedral mesh orientation check process can be used
        if KratosFluid.FluidMeshUtilities.AllElementsAreSimplex(fluid_computational_model_part):
            throw_errors = False # Skip the error thrown when swapping inverted geometries
            tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
            flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse()
            if self.assign_neighbour_elements:
                flags |= tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
            else:
                flags |= (tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS).AsFalse()
            KratosMultiphysics.TetrahedralMeshOrientationCheck(fluid_computational_model_part, throw_errors, flags).Execute()
        else:
            check_repeated_conditions = True
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Mesh contains non-simplex elements. Skin orientation cannot be checked.")
            KratosFluid.FluidMeshUtilities.AssignNeighbourElementsToConditions(fluid_computational_model_part, check_repeated_conditions)

class CheckAndPrepareModelProcess(CheckAndPrepareModelProcessFluid):
    def __init__(self, model, parameters):
        IssueDeprecationWarning('CheckAndPrepareModelProcess', '\'CheckAndPrepareModelProcess\' is deprecated. Use \'CheckAndPrepareModelProcessFluid\' instead.')
        super().__init__(model, parameters)