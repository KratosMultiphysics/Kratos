import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
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
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        iga_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        iga_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        volume_properties = model_part.GetProperties()[1]
        volume_properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, 100)
        volume_properties.SetValue(KratosMultiphysics.POISSON_RATIO, 0.2)
        volume_properties.SetValue(KratosMultiphysics.DENSITY, 100)
        volume_properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, cl)

        # Create nurbs volume and read Cad geometry
        run_modelers(current_model, modelers_list)

        nurbs_volume = model_part.GetGeometry("NurbsVolume")
        quadrature_point_geometries = KratosMultiphysics.GeometriesVector()
        nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2)
        for i in range(len(quadrature_point_geometries)):
            model_part.CreateNewElement('SmallDisplacementElement3D8N', i+1, quadrature_point_geometries[i], volume_properties)

        # Deform nurbs volume
        model_part = current_model.GetModelPart("NurbsMesh")
        nodal_disp = self._GetNodalDisplacement()
        nodal_temp = self._GetNodalTemperature()
        for i, node in enumerate(model_part.Nodes):
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, nodal_disp[i])
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, nodal_temp[i])

        # Map deformation to embedded geometry
        process_params = KratosMultiphysics.Parameters(
        """ {
                "main_model_part_name"                    : "NurbsMesh",
                "nurbs_volume_name"                       : "NurbsVolume",
                "embedded_model_part_name"                : "IgaModelPart",
                "nodal_results": ["DISPLACEMENT", "TEMPERATURE"],
                "gauss_point_results": ["CAUCHY_STRESS_VECTOR", "STRAIN_ENERGY", "CAUCHY_STRESS_TENSOR", "INTEGRATION_COORDINATES"]
        } """ )
        process = MapNurbsVolumeResultsToEmbeddedGeometryProcess(current_model, process_params)
        process.ExecuteBeforeOutputStep()

        lower_point = [-2.1, -2.1, -2.1]
        upper_point = [2.1, 2.1, 12.1]
        iga_model_part = current_model.GetModelPart("IgaModelPart")
        el_count = model_part.NumberOfElements()
        for i, node in enumerate(iga_model_part.Nodes):
            u = (node.X0 - lower_point[0]) / ( upper_point[0]- lower_point[0] )
            v = (node.Y0 - lower_point[1]) / ( upper_point[1]- lower_point[1] )
            w = (node.Z0 - lower_point[2]) / ( upper_point[2]- lower_point[2] )

            integration_points = [[u, v, w, 1]]
            quadrature_point_geometries = KratosMultiphysics.GeometriesVector()
            nurbs_volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 2, integration_points)

            el = model_part.CreateNewElement('SmallDisplacementElement3D8N', el_count+i+1, quadrature_point_geometries[0], volume_properties)
            el.Initialize(model_part.ProcessInfo)

            # Scalar
            energy_ref = el.CalculateOnIntegrationPoints(KratosMultiphysics.STRAIN_ENERGY, model_part.ProcessInfo)[0]
            energy = node.GetValue(KratosMultiphysics.STRAIN_ENERGY)
            self.assertAlmostEqual(energy, energy_ref)

            # Array
            coord_ref = el.CalculateOnIntegrationPoints(KratosMultiphysics.INTEGRATION_COORDINATES, model_part.ProcessInfo)[0]
            coord = node.GetValue(KratosMultiphysics.INTEGRATION_COORDINATES)
            self.assertVectorAlmostEqual(coord, coord_ref)

            # Vector
            stress_vector_ref = el.CalculateOnIntegrationPoints(KratosMultiphysics.CAUCHY_STRESS_VECTOR, model_part.ProcessInfo)[0]
            stress_vector = node.GetValue(KratosMultiphysics.CAUCHY_STRESS_VECTOR)
            self.assertVectorAlmostEqual(stress_vector, stress_vector_ref)

            # Matrix
            stress_tensor_ref = el.CalculateOnIntegrationPoints(KratosMultiphysics.CAUCHY_STRESS_TENSOR, model_part.ProcessInfo)[0]
            stress_tensor = node.GetValue(KratosMultiphysics.CAUCHY_STRESS_TENSOR)
            self.assertMatrixAlmostEqual(stress_tensor, stress_tensor_ref)


        self._CheckNodalResults(current_model)

        # #Create Reference solution
        # from KratosMultiphysics.json_output_process import JsonOutputProcess
        # out_parameters = KratosMultiphysics.Parameters("""
        # {
        #     "output_variables"     : ["DISPLACEMENT", "TEMPERATURE"],
        #     "output_file_name"     : "solid_cantilever_circular_cross_section/cylinder_results.json",
        #     "model_part_name"      : "IgaModelPart",
        #     "time_frequency"                : -0.1
        # }
        # """)

        # out = JsonOutputProcess(current_model, out_parameters)
        # out.ExecuteInitialize()
        # out.ExecuteBeforeSolutionLoop()
        # out.ExecuteFinalizeSolutionStep()

    def _CheckNodalResults(self, current_model):
        from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

        iga_model_part = current_model.GetModelPart("IgaModelPart")
        iga_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] =  0.2

        json_check_process_param = KratosMultiphysics.Parameters(
        """ {
                "check_variables": ["DISPLACEMENT", "TEMPERATURE"],
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

    def _GetNodalTemperature(self):
        temperature = [0.00328209, -0.00328382, -0.00421475, 0.00421482,  0.00421426, -0.00421456, -0.00327469,
                       0.00327906,  0.00581626, -0.00581616, 0.00437176, -0.00437125, -0.00437154,  0.00437127,
                      -0.00581791,  0.00581721, 0.000768276, -0.000768087, -0.000537143, 0.000537043, 0.000536351,
                      -0.000537709,-0.000767579,0.000766588, 0.000684277, -0.000682272, 0.000982665, -0.000982577,
                      -0.000983551, 0.000981584, -0.000686987, 0.000682113]
        return temperature

if __name__ == '__main__':
    KratosUnittest.main()
