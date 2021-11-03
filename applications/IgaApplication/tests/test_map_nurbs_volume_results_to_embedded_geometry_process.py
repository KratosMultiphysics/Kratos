import KratosMultiphysics
import KratosMultiphysics.IgaApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.IgaApplication.map_nurbs_volume_results_to_embedded_geometry_process import MapNurbsVolumeResultsToEmbeddedGeometryProcess

def run_modelers(current_model, modelers_list):
    from KratosMultiphysics.modeler_factory import KratosModelerFactory
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()

def run_processes(current_model, parameters):
    process = MapNurbsVolumeResultsToEmbeddedGeometryProcess(current_model, parameters)
    process.ExecuteBeforeOutputStep()

class TestMapNurbsVolumeResultsToEmbeddedGeometryProcess(KratosUnittest.TestCase):
    def test_MapNurbsVolumeResultsToEmbeddedGeometryProcess(self):
        current_model = KratosMultiphysics.Model()
        modelers_list = KratosMultiphysics.Parameters(
        """ [{
            "modeler_name": "NurbsGeometryModeler",
                "Parameters": {
                    "model_part_name" : "NurbsMesh",
                    "geometry_name"   : "NurbsVolume",
                    "lower_point": [-2.1, -2.1,  -2.1],
                    "upper_point": [2.1, 2.1, 12.1],
                    "polynomial_order" : [1, 2, 2],
                    "number_of_knot_spans" : [1,2,2]}
            },{
            "modeler_name":"CadIoModeler",
                "Parameters" : {"echo_level":0,
                                "cad_model_part_name": "IgaModelPart",
                                "geometry_file_name" : "solid_cantilever_circular_cross_section/cylinder.cad.json"}
            }] """ )

        model_part = current_model.CreateModelPart("NurbsMesh")
        iga_model_part = current_model.CreateModelPart("IgaModelPart")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        iga_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # Create nurbs volume and read Cad geometry
        run_modelers(current_model, modelers_list)

        # Deform nurbs volume
        model_part = current_model.GetModelPart("NurbsMesh")
        nodal_disp = self._GetNodalDisplacement()
        for i, node in enumerate(model_part.Nodes):
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, nodal_disp[i])

        # Map deformation to embedded geometry
        process_params = KratosMultiphysics.Parameters(
        """ {
                "main_model_part_name"                    : "NurbsMesh",
                "nurbs_volume_name"                       : "NurbsVolume",
                "embedded_model_part_name"                : "IgaModelPart",
                "nodal_results": ["DISPLACEMENT"]
        } """ )
        process = MapNurbsVolumeResultsToEmbeddedGeometryProcess(current_model, process_params)
        process.ExecuteBeforeOutputStep()

        self._CheckResults(current_model)

        # Create Reference solution
        # from KratosMultiphysics.json_output_process import JsonOutputProcess

        # out_parameters = KratosMultiphysics.Parameters("""
        # {
        #     "output_variables": ["DISPLACEMENT"],
        #     "output_file_name"     : "solid_cantilever_circular_cross_section/cylinder_results.json",
        #     "model_part_name"      : "IgaModelPart"
        # }
        # """)

        # out = JsonOutputProcess(current_model, out_parameters)
        # out.ExecuteInitialize()
        # out.ExecuteBeforeSolutionLoop()
        # out.ExecuteFinalizeSolutionStep()

    def _CheckResults(self, current_model):
        from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

        iga_model_part = current_model.GetModelPart("IgaModelPart")
        iga_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] =  0.2

        json_check_process_param = KratosMultiphysics.Parameters(
        """ {
                "check_variables": ["DISPLACEMENT"],
                "gauss_points_check_variables": [],
                "input_file_name": "solid_cantilever_circular_cross_section/cylinder_results.json",
                "model_part_name": "IgaModelPart",
                "time_frequency" : 0.0,
                "tolerance" : 1e-12
        } """)
        check_process = FromJsonCheckResultProcess(current_model, json_check_process_param)
        check_process.ExecuteInitialize()
        check_process.ExecuteBeforeSolutionLoop()
        check_process.ExecuteInitializeSolutionStep()
        check_process.ExecuteFinalizeSolutionStep()
        check_process.ExecuteFinalize()

    def _GetNodalDisplacement(self):
        displacement = [[0.00328209,0.00368871,-0.0245826],
                        [-0.00328382,0.00366463,-0.0245683],
                        [-0.00421475,0.00630134,-0.0113094],
                        [0.00421482,0.00630691,-0.0113134],
                        [0.00421426,0.00630818,0.0113092],
                        [-0.00421456,0.00630229,0.011313],
                        [-0.00327469,0.00365256,0.0245871],
                        [0.00327906,0.00368821,0.0245685],
                        [0.00581626,-0.0138629,0.0214322],
                        [-0.00581616,-0.0138561,0.0214254],
                        [0.00437176,-0.0195617,0.00928335],
                        [-0.00437125,-0.0195631,0.0092857],
                        [-0.00437154,-0.0195635,-0.00928335],
                        [0.00437127,-0.019562,-0.00928529],
                        [-0.00581791,-0.0138532,-0.0214323],
                        [0.00581721,-0.0138625,-0.0214261],
                        [0.000768276,0.132468,0.0636883],
                        [-0.000768087,0.132464,0.0636914],
                        [-0.000537143,0.132773,0.0325486],
                        [0.000537043,0.132775,0.032548],
                        [0.000536351,0.132775,-0.0325486],
                        [-0.000537709,0.132775,-0.0325476],
                        [-0.000767579,0.132461,-0.0636884],
                        [0.000766588,0.132469,-0.0636901],
                        [0.000684277,0.254191,0.0632124],
                        [-0.000682272,0.2542,0.0632054],
                        [0.000982665,0.25295,0.0286405],
                        [-0.000982577,0.252951,0.0286432],
                        [-0.000983551,0.252948,-0.0286405],
                        [0.000981584,0.252952,-0.0286429],
                        [-0.000686987,0.2542,-0.0632127],
                        [0.000682113,0.254195,-0.0632063]]

        return displacement

if __name__ == '__main__':
    KratosUnittest.main()
