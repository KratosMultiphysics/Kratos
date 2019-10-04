import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import time
import math


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineEmbeddedWakeProcess(Model, settings["Parameters"])


class DefineEmbeddedWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        # Call the base Kratos process constructor
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "epsilon": 1e-9
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.main_model_part = Model[settings["model_part_name"].GetString()].GetRootModelPart()
        self.wake_model_part = Model.CreateModelPart("wake")
        self.plane_model_part = Model.CreateModelPart("plane")
        self.output_model_part = Model.CreateModelPart("output")

        self.epsilon = settings["epsilon"].GetDouble()

    def ExecuteInitialize(self):
        ini_time = time.time()

        self._DefineWakeModelPart()

                # Find nodal neigbours util call
        avg_elem_num = 10
        avg_node_num = 10
        KratosMultiphysics.FindNodalNeighboursProcess(
            self.main_model_part, avg_elem_num, avg_node_num).Execute()

        # self._MoveAndRotateWake()
        # Executing define wake process
        # KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess3D(self.main_model_part, self.wake_model_part).Execute()
        CPFApp.DefineEmbeddedWakeProcess3D(self.main_model_part, self.wake_model_part).Execute()


        deactivated_ids = [19074, 4038, 1379, 19725, 14670]
        kutta_elem_ids = [1379, 5265, 8721, 9629, 14669, 14670, 19073, 19074, 19724, 19725]
        structure_elem_ids = [4038, 5268, 15312, 15635, 17141, 17901]
        trailing_edge_nodes = [9, 10]

        wake_elem_ids = [170, 191, 314, 707, 708, 720, 785, 867, 910, 911, 929, 1035, 1054, 1094, 1213, 1292, 1348, 1532, 1694, 1698, 1704, 1785, 1799, 1814, 1892, 2163, 2165, 2422, 2560, 2638, 2639, 2789, 2807, 2848, 3147, 3196, 3200, 3212, 3247, 3248, 3249, 3263, 3299, 3307, 3311, 3429, 3763, 3895, 4038, 4047, 4281, 4291, 4296, 4375, 4398, 4481, 4755, 4878, 4886, 5026, 5255, 5268, 5445, 5454, 5456, 5560, 5719, 5736, 5828, 5842, 5955, 5960, 5962, 5975, 5978, 6012, 6057, 6080, 6081, 6082, 6156, 6273, 6415, 6460, 6569, 6740, 6904, 6990, 7272, 7523, 7853, 7879, 8009, 8010, 8011, 8090, 8235, 8282, 8352, 8422, 8445, 8476, 8508, 8646, 8652, 8769, 8934, 8983, 8984, 9016, 9167, 9218, 9402, 9448, 9451, 9508, 9623, 9631, 9665, 9769, 9771, 9772, 9932, 9934, 10085, 10357, 10367, 10590, 10594, 10643, 10720, 11008, 11131, 11171, 11188, 11352, 11632, 11677, 11705, 11706, 11710, 11712, 12245, 12424, 12636, 12742, 12754, 12847, 12874, 12883, 12885, 12919, 13076, 13137, 13144, 13211, 13339, 13553, 13560, 13714, 13801, 13802, 13820, 14409, 14431, 14500, 14778, 14781, 14782, 14896, 14979, 15002, 15312, 15377, 15447, 15544, 15635, 15637, 15644, 15771, 15773, 15940, 15968, 15977, 16029, 16144, 16160, 16265, 16390, 16410, 16417, 16420, 16523, 16570, 16581, 16673, 16854, 16987, 17033, 17037, 17089, 17092, 17141, 17212, 17213, 17362, 17562, 17563, 17564, 17565, 17578, 17580, 17722, 17819, 17820, 17824, 17879, 17901, 18074, 18140, 18296, 18298, 18299, 18306, 18307, 18780, 18907, 19442, 19443, 19967, ]

        wake_distances=[]
        with open('3d_wake_distances.dat') as dat_file:
            lines=dat_file.readlines()
            for line in lines:
                dist=[]
                dist.append(float(line.split(' ')[1]))
                dist.append(float(line.split(' ')[2]))
                dist.append(float(line.split(' ')[3]))
                dist.append(float(line.split(' ')[4]))
                wake_distances.append(dist)


        # for elem_id in deactivated_ids:
        #     self.main_model_part.GetElement(elem_id).Set(KratosMultiphysics.ACTIVE, True)
        # for ix, elem_id in enumerate(wake_elem_ids):
        #     self.main_model_part.GetElement(elem_id).SetValue(CPFApp.WAKE, True)
        #     # print(wake_distances[ix])
        #     self.main_model_part.GetElement(elem_id).SetValue(CPFApp.WAKE_ELEMENTAL_DISTANCES, wake_distances[ix])

        #     counter = 0
        #     for node in self.main_model_part.GetElement(elem_id).GetNodes():
        #         node.SetValue(CPFApp.WAKE_DISTANCE,wake_distances[ix][counter])
        #         counter += 1



        # for elem_id in kutta_elem_ids:
        #     self.main_model_part.GetElement(elem_id).SetValue(CPFApp.KUTTA, True)

        # for elem_id in structure_elem_ids:
        #     self.main_model_part.GetElement(elem_id).Set(KratosMultiphysics.STRUCTURE, True)

        # for node_id in trailing_edge_nodes:
        #     self.main_model_part.GetNode(node_id).SetValue(CPFApp.TRAILING_EDGE, True)


        self.wake_sub_model_part = self.main_model_part.CreateSubModelPart('aux')
        for elem in self.main_model_part.Elements:
            if elem.GetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WAKE) == 1:
                self.wake_sub_model_part.Elements.append(elem)
                # print (elem.Id, end = ' ', flush  = True)


        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Wake computation time: ',time.time()-ini_time)

    def _DefineWakeModelPart(self):
        ''' This function generates the modelpart of the wake. TODO: make end of the domain user-definable.
        '''
        self.wake_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.wake_model_part.CreateNewNode(2, 0.0, 0.5, 0.0)

        self.wake_model_part.CreateNewNode(3, 200.0*math.cos(math.radians(-5.0)), 0.5, 200.0*math.sin(math.radians(-5.0)))
        self.wake_model_part.CreateNewNode(4, 200.0*math.cos(math.radians(-5.0)), 0.0, 200.0*math.sin(math.radians(-5.0)))

        # self.wake_model_part.CreateNewNode(1, 0.4980874, -0.1, -0.04357787)
        # self.wake_model_part.CreateNewNode(2, 0.4980874, 0.6, -0.04357787)
        # self.wake_model_part.CreateNewNode(3, 200.0, 0.6, -0.04357787)
        # self.wake_model_part.CreateNewNode(4, 200.0, -0.1, -0.04357787)

        self.wake_model_part.CreateNewElement("Element3D3N", 1, [1,3,2], KratosMultiphysics.Properties(0))
        self.wake_model_part.CreateNewElement("Element3D3N", 2, [1,4,3], KratosMultiphysics.Properties(0))

    def _MoveAndRotateWake(self):
        ''' This function moves and rotates the wake with the same parameters as the geometry.
        '''
        self.moving_parameters = KratosMultiphysics.Parameters()
        self.moving_parameters.AddEmptyValue("origin")
        self.moving_parameters["origin"].SetVector(self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN))
        self.moving_parameters.AddEmptyValue("rotation_angle")
        # angle=math.radians(-self.main_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE))
        angle=math.radians(-5.0)
        self.moving_parameters["rotation_angle"].SetDouble(angle)
        CPFApp.MoveModelPartProcess(self.wake_model_part, self.moving_parameters).Execute()

    def ExecuteFinalizeSolutionStep(self):
        CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled3D(self.wake_sub_model_part, 1e-1, 0)