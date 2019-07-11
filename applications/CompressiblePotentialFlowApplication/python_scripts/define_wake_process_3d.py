import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
from stl import mesh #this requires numpy-stl

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess3D(Model, settings["Parameters"])

class DefineWakeProcess3D(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "wake_stl_file_name" : "",
            "epsilon": 1e-9
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

        trailing_edge_model_part_name = settings["model_part_name"].GetString()
        if trailing_edge_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the trailing edge nodes"
            raise Exception(err_msg)
        self.trailing_edge_model_part = Model[trailing_edge_model_part_name]

        self.wake_stl_file_name = settings["wake_stl_file_name"].GetString()
        if self.wake_stl_file_name == "":
            err_msg = "Empty wake_stl_file_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the stl file name that contains the wake surface nodes"
            raise Exception(err_msg)

        self.epsilon = settings["epsilon"].GetDouble()

        self.fluid_model_part = self.trailing_edge_model_part.GetRootModelPart()

    def ExecuteInitialize(self):

        # Mark trailing edge nodes
        for node in self.trailing_edge_model_part.Nodes:
            node.SetValue(CPFApp.TRAILING_EDGE, True)

        self.__CreateWakeModelPart()
        self.__MarkWakeElements()

    # This function imports the stl file containing the wake and creates the wake model part out of it.
    def __CreateWakeModelPart(self):

        wake_stl_mesh = mesh.Mesh.from_multi_file(self.wake_stl_file_name)
        self.wake_model_part = self.model.CreateModelPart("wake_model_part")

        #self.wake_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        #self.wake_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        dummy_property = self.wake_model_part.Properties[0]
        node_id = 1
        elem_id = 1

        # Looping over stl meshes
        for stl_mesh in wake_stl_mesh:
            for vertex in stl_mesh.points:
                node1 = self.wake_model_part.CreateNewNode(node_id, float(vertex[0]), float(vertex[1]), float(vertex[2]))
                node_id+=1
                node2 = self.wake_model_part.CreateNewNode(node_id, float(vertex[3]), float(vertex[4]), float(vertex[5]))
                node_id+=1
                node3 = self.wake_model_part.CreateNewNode(node_id, float(vertex[6]), float(vertex[7]), float(vertex[8]))
                node_id+=1

                self.wake_model_part.CreateNewElement("Element3D3N", elem_id,  [
                                              node1.Id, node2.Id, node3.Id], dummy_property)
                elem_id += 1

    # Check which elements are cut and mark them as wake
    def __MarkWakeElements(self):
        distance_calculator = KratosMultiphysics.CalculateSignedDistanceTo3DSkinProcess(
            self.wake_model_part, self.fluid_model_part)
        distance_calculator.Execute()

        for elem in self.fluid_model_part.Elements:
            if(elem.Is(KratosMultiphysics.TO_SPLIT)):
                elem.SetValue(CPFApp.WAKE, True)
                wake_elemental_distances = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                for i in range(len(wake_elemental_distances)):
                    if(abs(wake_elemental_distances[i]) < self.epsilon ):
                        if(wake_elemental_distances[i] < 0.0):
                            wake_elemental_distances[i] = -self.epsilon
                        else:
                            wake_elemental_distances[i] = self.epsilon
                elem.SetValue(CPFApp.WAKE_ELEMENTAL_DISTANCES,wake_elemental_distances)







