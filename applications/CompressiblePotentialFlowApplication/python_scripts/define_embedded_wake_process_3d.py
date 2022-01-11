import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import time
import KratosMultiphysics.MeshingApplication as MeshingApplication
import math
from KratosMultiphysics.gid_output_process import GiDOutputProcess


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
            "epsilon": 1e-9,
            "refinement_iterations" : 0,
            "target_h_wake" : 1.0
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

        self.main_model_part = Model[settings["model_part_name"].GetString()].GetRootModelPart()
        self.wake_model_part = Model.CreateModelPart("wake")
        # self.plane_model_part = Model.CreateModelPart("plane")
        # self.output_model_part = Model.CreateModelPart("output")
        self.wake_normal = KratosMultiphysics.Vector(3, 0.0)
        self.wake_normal[2] = 1.0
        self.main_model_part.ProcessInfo.SetValue(CPFApp.WAKE_NORMAL,self.wake_normal)
        self.refinement_iterations = settings["refinement_iterations"].GetInt()
        self.target_h_wake = settings["target_h_wake"].GetDouble()

        self.epsilon = settings["epsilon"].GetDouble()

    def ExecuteInitialize(self):
        ini_time = time.time()

        self._DefineWakeModelPart()

        # self._MoveAndRotateWake()
        # Executing define wake process
        # KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess3D(self.main_model_part, self.wake_model_part).Execute()
        # Find nodal neigbours util call
        KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.main_model_part).Execute()

        CPFApp.DefineEmbeddedWakeProcess3D(self.main_model_part, self.wake_model_part).Execute()

        # self.skin_model_part = self.model.GetModelPart("skin")

        for i in range(self.refinement_iterations):
            find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
            find_nodal_h.Execute()


            for node in self.main_model_part.Nodes:
                node.Set(KratosMultiphysics.BLOCKED)
                this_h = node.GetValue(KratosMultiphysics.NODAL_H)
                tensor = KratosMultiphysics.Vector(6, 0.0)
                tensor[0] = 1e-4
                tensor[1] = 1e-4
                tensor[2] = 1e-4
                node.SetValue(MeshingApplication.METRIC_TENSOR_3D, tensor)

            for elem in self.main_model_part.Elements:
                if elem.GetValue(CPFApp.WAKE):
                    for node in elem.GetNodes():
                        node.Set(KratosMultiphysics.BLOCKED, False)
                        tensor = KratosMultiphysics.Vector(6, 0.0)
                        tensor[0] = 1/self.target_h_wake/self.target_h_wake
                        tensor[1] = 1/self.target_h_wake/self.target_h_wake
                        tensor[2] = 1/self.target_h_wake/self.target_h_wake
                        old_tensor = node.GetValue(MeshingApplication.METRIC_TENSOR_3D)
                        if old_tensor.norm_2() < tensor.norm_2():
                            node.SetValue(MeshingApplication.METRIC_TENSOR_3D, tensor)

            mmg_parameters = KratosMultiphysics.Parameters("""
            {
                "discretization_type"              : "STANDARD",
                "save_external_files"              : false,
                "initialize_entities"              : false,
                "preserve_flags"                   : false,
                "interpolate_nodal_values"                   : false,
                "echo_level"                       : 0
            }
            """)

            mmg_process =MeshingApplication.MmgProcess3D(self.main_model_part, mmg_parameters)
            mmg_process.Execute()

            self._CalculateDistance()
            self._ModifyFinalDistance()

            KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.main_model_part).Execute()

            CPFApp.DefineEmbeddedWakeProcess3D(self.main_model_part, self.wake_model_part).Execute()

            # self.target_h_wake = self.target_h_wake*0.5


        KratosMultiphysics.Logger.PrintInfo('EmbeddedWake','Wake computation time: ',time.time()-ini_time)
        model_part_to_write=self.main_model_part
        KratosMultiphysics.ModelPartIO('./Meshes/after_wake_remeshed', KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(model_part_to_write)

        gid_output = GiDOutputProcess(self.main_model_part,
                                    "after_wake_process",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostAscii",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                            "nodal_results"       : ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL","GEOMETRY_DISTANCE"],
                                            "nodal_nonhistorical_results": ["TRAILING_EDGE","WAKE_DISTANCE","LOWER_SURFACE","UPPER_SURFACE","TEMPERATURE","KUTTA","WAKE"],
                                            "gauss_point_results" : ["PRESSURE_COEFFICIENT","VELOCITY","WAKE","KUTTA","DENSITY","VELOCITY_LOWER", "LOWER_WAKE"],
                                            "nodal_flags_results": ["MARKER"],
                                            "elemental_conditional_flags_results": ["TO_SPLIT","STRUCTURE","BOUNDARY","SOLID","MARKER"]
                                            }
                                        }
                                        """)
                                    )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()


        # for elem in self.main_model_part.Elements:
        #     for node in elem.GetNodes():
        #         if node.Id == 120883:
        #             print("ID WAKE", elem.Id, elem.GetValue(CPFApp.WAKE))

        # print("ELEM IDS")
        # for elem in self.main_model_part.Elements:
        #     for node in elem.GetNodes():
        #         if node.Id == 120883:
        #             print(elem.Id, end=" ")
        # print("FINISHED IDS")


    def _CalculateDistance(self):
        ''' This function calculate the distance to skin for every node in the main_model_part.'''
        ini_time=time.time()
        KratosMultiphysics.CalculateDistanceToSkinProcess3D(self.main_model_part, self.skin_model_part, 1e-16).Execute()
        KratosMultiphysics.Logger.PrintInfo('WakeProcess3D','CalculateDistance time: ',time.time()-ini_time)

    def _ModifyFinalDistance(self):
        ''' This function modifies the distance field to avoid ill defined cuts.
        '''
        distance_modification_parameters = KratosMultiphysics.Parameters("""{
                "distance_threshold"                          : 0.001,
                "check_at_each_time_step"                     : true,
                "avoid_almost_empty_elements"                 : true,
                "deactivate_full_negative_elements"           : true,
                "full_negative_elements_fixed_variables_list" : ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"]
            }
        """)

        ini_time = time.time()
        KratosMultiphysics.FluidDynamicsApplication.DistanceModificationProcess(self.main_model_part,distance_modification_parameters).Execute()
        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Modify distance time: ',time.time()-ini_time)
        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.DISTANCE,CPFApp.GEOMETRY_DISTANCE, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.DISTANCE, self.main_model_part.Nodes)

    def _DefineWakeModelPart(self):
        ''' This function generates the modelpart of the wake. TODO: make end of the domain user-definable.
        '''
        mid_chord_0 = [ 0.0, -0.01, 0.0]
        mid_chord_2 = [0.0, 2.0, 0.0]
        te_point_0 = [0.498097, -0.01, -0.0435779]
        te_point_2 = [0.498097, 2.0, -0.0435779]
        weight = 0.5
        self.wake_model_part.CreateNewNode(1, (1-weight)*mid_chord_0[0]+weight*te_point_0[0], -0.01, (1-weight)*mid_chord_0[2]+weight*te_point_0[2])
        self.wake_model_part.CreateNewNode(2, (1-weight)*mid_chord_2[0]+weight*te_point_2[0], 2.0, (1-weight)*mid_chord_2[2]+weight*te_point_2[2])


        # self.wake_model_part.CreateNewNode(3, 200.0*math.cos(math.radians(-5.0)), 2.0, 200.0*math.sin(math.radians(-5.0)))
        # self.wake_model_part.CreateNewNode(4, 200.0*math.cos(math.radians(-5.0)), 0.0, 200.0*math.sin(math.radians(-5.0)))
        # self.wake_model_part.CreateNewElement("Element3D3N", 1, [1,3,2], KratosMultiphysics.Properties(0))
        # self.wake_model_part.CreateNewElement("Element3D3N", 2, [1,4,3], KratosMultiphysics.Properties(0))
        self.wake_model_part.CreateNewNode(3, te_point_2[0], te_point_2[1], te_point_2[2])
        self.wake_model_part.CreateNewNode(4, te_point_0[0], te_point_0[1], te_point_0[2])
        self.wake_model_part.CreateNewNode(5, 200.0, 2.0, 0.0)
        self.wake_model_part.CreateNewNode(6, 200.0, 0.0, 0.0)
        print(self.wake_model_part.GetNode(1))
        print(self.wake_model_part.GetNode(2))
        print(self.wake_model_part.GetNode(3))
        print(self.wake_model_part.GetNode(4))

        self.wake_model_part.CreateNewElement("Element3D3N", 1, [1,3,2], KratosMultiphysics.Properties(0))
        self.wake_model_part.CreateNewElement("Element3D3N", 2, [1,4,3], KratosMultiphysics.Properties(0))
        self.wake_model_part.CreateNewElement("Element3D3N", 3, [4,5,3], KratosMultiphysics.Properties(0))
        self.wake_model_part.CreateNewElement("Element3D3N", 4, [4,6,5], KratosMultiphysics.Properties(0))


        # self.wake_model_part.CreateNewNode(1, 0.498097, 2.0, -0.0435779)
        # self.wake_model_part.CreateNewNode(2, 0.498097, -0.01, -0.0435779)
        # self.wake_model_part.CreateNewNode(3, 200.0, 2.0, 0.0)
        # self.wake_model_part.CreateNewNode(4, 200.0, 0.0, 0.0)

        # self.wake_model_part.CreateNewElement("Element3D3N", 1, [1,3,2], KratosMultiphysics.Properties(0))
        # self.wake_model_part.CreateNewElement("Element3D3N", 2, [1,4,3], KratosMultiphysics.Properties(0))

    def _MoveAndRotateWake(self):
        ''' This function moves and rotates the wake with the same parameters as the geometry.
        '''
        self.moving_parameters = KratosMultiphysics.Parameters()
        self.moving_parameters.AddEmptyValue("origin")
        self.moving_parameters["origin"].SetVector(self.main_model_part.ProcessInfo.GetValue(CPFApp.WAKE_ORIGIN))
        self.moving_parameters.AddEmptyValue("rotation_angle")
        angle=math.radians(-self.main_model_part.ProcessInfo.GetValue(CPFApp.ROTATION_ANGLE))
        self.moving_parameters["rotation_angle"].SetDouble(angle)
        CPFApp.MoveModelPartProcess(self.wake_model_part, self.moving_parameters).Execute()

    # def ExecuteFinalizeSolutionStep(self):
    #     CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled3D(self.wake_sub_model_part, 1e-1, 0)