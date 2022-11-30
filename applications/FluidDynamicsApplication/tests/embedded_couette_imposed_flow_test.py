import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class EmbeddedCouetteImposedTestFluidDynamicsAnalysis(FluidDynamicsAnalysis):
    def __init__(self, model, parameters, embedded_velocity, distance):
        super().__init__(model, parameters)
        self.distance = distance
        self.embedded_velocity = embedded_velocity

    def ModifyAfterSolverInitialize(self):
        super().ModifyAfterSolverInitialize()

        # Set up distance field in the DISTANCE nodal variable
        computing_model_part = self._GetSolver().GetComputingModelPart()
        domain_size = computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size == 2:
            distance = lambda node : self.distance - node.Y
        elif domain_size == 3:
            distance = lambda node : self.distance - node.Z
        else:
            raise Exception("DOMAIN_SIZE has not been set in the ProcessInfo container.")

        for node in computing_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance(node))

        # Set the ELEMENTAL_DISTANCES value
        # This is required for the discontinuous elements tests
        n_nodes = len(computing_model_part.Elements[1].GetNodes())
        n_edges = computing_model_part.Elements[1].GetGeometry().EdgesNumber()
        elem_edge_dist = KratosMultiphysics.Vector(n_edges, -1.0)
        for element in computing_model_part.Elements:
            elem_dist = KratosMultiphysics.Vector(n_nodes)
            elem_nodes = element.GetNodes()
            for i_node in range(n_nodes):
                elem_dist[i_node] = elem_nodes[i_node].GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            element.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, elem_dist)
            element.SetValue(KratosMultiphysics.ELEMENTAL_EDGE_DISTANCES, elem_edge_dist)

        # Set the EMBEDDED_VELOCITY in the intersected elements nodes
        for elem in computing_model_part.Elements:
            n_pos = 0
            n_neg = 0
            for node in elem.GetNodes():
                if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                    n_neg += 1
                else:
                    n_pos += 1

            if ((n_pos != 0) and (n_neg != 0)):
                for node in elem.GetNodes():
                    node.SetValue(KratosMultiphysics.EMBEDDED_VELOCITY, [self.embedded_velocity, 0.0, 0.0])

        # Set up the initial velocity condition
        # Note that we also initialize the buffer to avoid artificial inertial residual
        buffer_size = computing_model_part.GetBufferSize()
        if domain_size == 2:
            vel_x = lambda node : (self.embedded_velocity/self.distance)*node.Y
        else:
            vel_x = lambda node : (self.embedded_velocity/self.distance)*node.Z

        for node in computing_model_part.Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if dist > 0.0:
                init_v = KratosMultiphysics.Vector([vel_x(node),0.0,0.0])
                for step in range(buffer_size):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, step, init_v)

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        # Set velocity boundary condition
        main_model_part = self._GetSolver().GetComputingModelPart().GetRootModelPart()
        domain_size = main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size == 2:
            vel_x = lambda node : (self.embedded_velocity/self.distance)*node.Y
        else:
            vel_x = lambda node : (self.embedded_velocity/self.distance)*node.Z

        for node in main_model_part.GetSubModelPart("Inlet").Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if dist > 0.0:
                aux_vel = KratosMultiphysics.Vector([vel_x(node),0.0,0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)

@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
class EmbeddedCouetteImposedFlowTest(KratosUnittest.TestCase):

    # Embedded element tests
    def testEmbeddedCouetteImposedFlowEmbeddedWeaklyCompressibleNavierStokes2D(self):
        self.distance = 1.65
        self.embedded_velocity = 2.56
        self.settings_filename = "ProjectParameters2D.json"
        self.formulation_settings = KratosMultiphysics.Parameters("""{
            "element_type": "embedded_weakly_compressible_navier_stokes",
            "is_slip": false,
            "dynamic_tau": 1.0
        }""")
        self._ReadAndCustomizeTestSettings()
        self._RunEmbeddedCouetteFlowTest()

    def testEmbeddedCouetteImposedFlowEmbeddedWeaklyCompressibleNavierStokes3D(self):
        self.distance = 1.65
        self.embedded_velocity = 2.56
        self.settings_filename = "ProjectParameters3D.json"
        self.formulation_settings = KratosMultiphysics.Parameters("""{
            "element_type": "embedded_weakly_compressible_navier_stokes",
            "is_slip": false,
            "dynamic_tau": 1.0
        }""")
        self._ReadAndCustomizeTestSettings()
        self._RunEmbeddedCouetteFlowTest()

    def setUp(self):
        self.print_output = False
        self.check_tolerance = 1.0e-6
        self.check_relative_tolerance = 1.0e-8
        self.print_reference_values = False
        self.work_folder = "embedded_couette_imposed_flow_test"

    def _RunEmbeddedCouetteFlowTest(self):
        # Create the test simulation
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            simulation = EmbeddedCouetteImposedTestFluidDynamicsAnalysis(self.model, self.parameters, self.embedded_velocity, self.distance)
            simulation.Run()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('{0}.time'.format(self.parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()))

    def _ReadAndCustomizeTestSettings(self):
        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(self.settings_filename,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # Customize simulation settings
        self.parameters["solver_settings"]["formulation"] = self.formulation_settings

        # If required, add the output process to the test settings
        if self.print_output:
            self._AddOutput()

        # If required, add the reference values output process to the test settings
        if self.print_reference_values:
            self._AddReferenceValuesOutput()
        else:
            self._AddReferenceValuesCheck()

    def _AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "MainModelPart.fluid_computational_model_part",
                "output_name"            : "embedded_couette_imposed_flow_test",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "time",
                        "output_control_type"         : "step",
                        "output_interval"             : 1,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["DISTANCE","VELOCITY","PRESSURE"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        domain_size = self.parameters["solver_settings"]["domain_size"].GetInt()
        output_name = gid_output_settings["Parameters"]["output_name"].GetString()
        output_name += "_{0}_{1}d".format(self.formulation_settings["element_type"].GetString(), domain_size)
        gid_output_settings["Parameters"]["output_name"].SetString(output_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

    def _AddReferenceValuesOutput(self):
        json_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "json_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "JsonOutputProcess",
            "Parameters"    : {
                "output_variables" : ["VELOCITY","PRESSURE"],
                "output_file_name" : "reference_embedded_couette_imposed_flow_test",
                "model_part_name"  : "MainModelPart.Parts_Fluid",
                "time_frequency"   : 199.0
            }
        }""")
        domain_size = self.parameters["solver_settings"]["domain_size"].GetInt()
        output_file_name = json_output_settings["Parameters"]["output_file_name"].GetString()
        output_file_name += "_{0}_{1}d".format(self.formulation_settings["element_type"].GetString(), domain_size)
        json_output_settings["Parameters"]["output_file_name"].SetString(output_file_name)
        self.parameters["processes"]["json_check_process_list"].Append(json_output_settings)

    def _AddReferenceValuesCheck(self):
        json_check_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "from_json_check_result_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "FromJsonCheckResultProcess",
            "Parameters"    : {
                "check_variables"      : ["VELOCITY_X","VELOCITY_Y","PRESSURE"],
                "input_file_name"      : "reference_embedded_couette_imposed_flow_test",
                "model_part_name"      : "MainModelPart.Parts_Fluid",
                "tolerance"            : 0.0,
                "relative_tolerance"   : 0.0,
                "time_frequency"       : 199.0
            }
        }""")
        domain_size = self.parameters["solver_settings"]["domain_size"].GetInt()
        input_file_name = json_check_settings["Parameters"]["input_file_name"].GetString()
        input_file_name += "_{0}_{1}d".format(self.formulation_settings["element_type"].GetString(), domain_size)
        json_check_settings["Parameters"]["input_file_name"].SetString(input_file_name)
        json_check_settings["Parameters"]["tolerance"].SetDouble(self.check_tolerance)
        json_check_settings["Parameters"]["relative_tolerance"].SetDouble(self.check_relative_tolerance)
        self.parameters["processes"]["json_check_process_list"].Append(json_check_settings)

if __name__ == '__main__':
    KratosUnittest.main()
