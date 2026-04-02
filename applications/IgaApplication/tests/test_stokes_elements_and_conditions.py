import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as DFA
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.modeler_factory import KratosModelerFactory
try:
    from .test_creation_utility import TestCreationUtility
except ImportError:
    from test_creation_utility import TestCreationUtility


def _run_modelers(current_model, modelers_list):
    factory = KratosModelerFactory()
    list_of_modelers = factory.ConstructListOfModelers(current_model, modelers_list)

    for modeler in list_of_modelers:
        modeler.SetupGeometryModel()

    for modeler in list_of_modelers:
        modeler.PrepareGeometryModel()

    for modeler in list_of_modelers:
        modeler.SetupModelPart()


def _create_cube_outer_skin(model_part, min_coord=0.25, max_coord=1.75):
    model_part.CreateNewProperties(1)
    prop = model_part.GetProperties()[1]

    coords = {
        1: (min_coord, min_coord, min_coord),
        2: (max_coord, min_coord, min_coord),
        3: (max_coord, max_coord, min_coord),
        4: (min_coord, max_coord, min_coord),
        5: (min_coord, min_coord, max_coord),
        6: (max_coord, min_coord, max_coord),
        7: (max_coord, max_coord, max_coord),
        8: (min_coord, max_coord, max_coord),
    }

    for node_id, (x, y, z) in coords.items():
        model_part.CreateNewNode(node_id, x, y, z)

    triangles = {
        1: [1, 3, 2],
        2: [1, 4, 3],
        3: [5, 6, 7],
        4: [5, 7, 8],
        5: [1, 2, 6],
        6: [1, 6, 5],
        7: [4, 8, 7],
        8: [4, 7, 3],
        9: [1, 5, 8],
        10: [1, 8, 4],
        11: [2, 7, 6],
        12: [2, 3, 7],
    }

    for cond_id, node_ids in triangles.items():
        model_part.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)


def _create_support_condition_model_part_3d(condition_name):
    current_model = KM.Model()
    iga_model_part = current_model.CreateModelPart("IgaModelPart")
    iga_model_part.SetBufferSize(2)
    iga_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
    iga_model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
    iga_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

    skin_model_part_outer_initial = current_model.CreateModelPart("skin_model_part_outer_initial")
    skin_model_part_outer_initial.AddNodalSolutionStepVariable(KM.VELOCITY)
    _create_cube_outer_skin(skin_model_part_outer_initial)

    modeler_settings = KM.Parameters(
        f"""
        [
            {{
                "modeler_name": "NurbsGeometryModelerSbm",
                "Parameters": {{
                    "model_part_name" : "IgaModelPart",
                    "lower_point_xyz": [0.0, 0.0, 0.0],
                    "upper_point_xyz": [2.0, 2.0, 2.0],
                    "lower_point_uvw": [0.0, 0.0, 0.0],
                    "upper_point_uvw": [2.0, 2.0, 2.0],
                    "polynomial_order" : [1, 1, 1],
                    "number_of_knot_spans" : [4, 4, 4],
                    "lambda_outer": 0.5,
                    "number_of_inner_loops": 0,
                    "skin_model_part_outer_initial_name": "skin_model_part_outer_initial",
                    "skin_model_part_name": "skin_model_part",
                    "echo_level": 0
                }}
            }},
            {{
                "modeler_name": "IgaModelerSbm",
                "Parameters": {{
                    "echo_level": 0,
                    "skin_model_part_name": "skin_model_part",
                    "analysis_model_part_name": "IgaModelPart",
                    "element_condition_list": [
                        {{
                            "geometry_type": "SurfaceEdge",
                            "iga_model_part": "Support_outer",
                            "type": "condition",
                            "name": "{condition_name}",
                            "shape_function_derivatives_order": 3,
                            "sbm_parameters": {{
                                "is_inner" : false
                            }}
                        }}
                    ]
                }}
            }}
        ]
        """
    )

    _run_modelers(current_model, modeler_settings)
    support_model_part = current_model.GetModelPart("IgaModelPart.Support_outer")
    return current_model, support_model_part, next(iter(support_model_part.Conditions))


class FluidTests(KratosUnittest.TestCase):

    @staticmethod
    def create_element(model_part, integration_point):
        props = model_part.GetProperties()[1]
        props.SetValue(KM.DYNAMIC_VISCOSITY, 1.0)

        law = DFA.Newtonian2DLaw()
        props.SetValue(KM.CONSTITUTIVE_LAW, law)
        geometry = TestCreationUtility.GetQuadraturePointGeometryP2(model_part, integration_point)
        element = model_part.CreateNewElement("StokesElement", 1, geometry, props)

        bf = KM.Vector(3)
        bf[0] = 30.0
        bf[1] = 0.0
        bf[2] = 0.0
        element.SetValue(KM.BODY_FORCE, bf)
        divergence_stress = KM.Vector(2)
        divergence_stress[0] = 0.1
        divergence_stress[1] = 100.0
        element.SetValue(IGA.DIVERGENCE_STRESS, divergence_stress)

        return element

    @staticmethod
    def create_navier_stokes_element_3d(model_part, integration_point):
        props = model_part.CreateNewProperties(1)
        props.SetValue(KM.DENSITY, 2.5)
        props.SetValue(KM.DYNAMIC_VISCOSITY, 0.8)
        props.SetValue(KM.CONSTITUTIVE_LAW, DFA.Newtonian3DLaw())

        volume = TestCreationUtility.GenerateNurbsVolumeP2Rectangular(model_part)
        volume.SetId(1)
        quadrature_point_geometries = KM.GeometriesVector()
        volume.CreateQuadraturePointGeometries(quadrature_point_geometries, 3)
        model_part.AddGeometry(volume)

        geometry = quadrature_point_geometries[0]
        element = model_part.CreateNewElement("NavierStokesElement", 1, geometry, props)

        bf = KM.Vector(3)
        bf[0] = 1.1
        bf[1] = -0.7
        bf[2] = 0.25
        element.SetValue(KM.BODY_FORCE, bf)

        divergence_stress = FluidTests._compute_divergence_stress_3d(
            volume,
            quadrature_point_geometries)
        element.SetValue(IGA.DIVERGENCE_STRESS, divergence_stress)

        return element

    @staticmethod
    def _compute_divergence_stress_3d(volume, quadrature_point_geometries):
        num_gauss_points = len(quadrature_point_geometries)
        sigma_values = KM.Matrix(num_gauss_points, 6)
        shape_function_values = KM.Matrix(num_gauss_points, num_gauss_points)
        shape_function_values_dx = KM.Matrix(num_gauss_points, num_gauss_points)
        shape_function_values_dy = KM.Matrix(num_gauss_points, num_gauss_points)
        shape_function_values_dz = KM.Matrix(num_gauss_points, num_gauss_points)

        sigma_coefficients_xx = [0.1 * node.X for node in volume]
        sigma_coefficients_yy = [-0.15 * node.Y for node in volume]
        sigma_coefficients_zz = [0.1 * node.Z for node in volume]

        for i, quadrature_point_geometry in enumerate(quadrature_point_geometries):
            shape_functions = quadrature_point_geometry.ShapeFunctionsValues()
            shape_function_derivatives = quadrature_point_geometry.ShapeFunctionDerivatives(1, 0)

            sigma_xx = 0.0
            sigma_yy = 0.0
            sigma_zz = 0.0

            for j in range(num_gauss_points):
                shape_function_value = shape_functions[0, j]
                shape_function_values[i, j] = shape_function_value
                shape_function_values_dx[i, j] = shape_function_derivatives[j, 0]
                shape_function_values_dy[i, j] = shape_function_derivatives[j, 1]
                shape_function_values_dz[i, j] = shape_function_derivatives[j, 2]

                sigma_xx += shape_function_value * sigma_coefficients_xx[j]
                sigma_yy += shape_function_value * sigma_coefficients_yy[j]
                sigma_zz += shape_function_value * sigma_coefficients_zz[j]

            sigma_values[i, 0] = sigma_xx
            sigma_values[i, 1] = sigma_yy
            sigma_values[i, 2] = sigma_zz
            sigma_values[i, 3] = 0.0
            sigma_values[i, 4] = 0.0
            sigma_values[i, 5] = 0.0

        div_sigma_utility = DFA.ComputeDivSigmaUtility()
        divergence_matrix = div_sigma_utility.ComputeDivergence(
            sigma_values,
            shape_function_values,
            shape_function_values_dx,
            shape_function_values_dy,
            shape_function_values_dz)

        divergence_stress = KM.Vector(3)
        divergence_stress[0] = divergence_matrix[0, 0]
        divergence_stress[1] = divergence_matrix[0, 1]
        divergence_stress[2] = divergence_matrix[0, 2]

        return divergence_stress

    @staticmethod
    def create_stokes_element_3d(model_part, integration_point):
        props = model_part.CreateNewProperties(1)
        props.SetValue(KM.DYNAMIC_VISCOSITY, 0.8)
        props.SetValue(KM.CONSTITUTIVE_LAW, DFA.Newtonian3DLaw())

        geometry = TestCreationUtility.GetQuadraturePointGeometryFromRectangularVolumeP2(model_part, integration_point)
        element = model_part.CreateNewElement("StokesElement", 1, geometry, props)

        bf = KM.Vector(3)
        bf[0] = 1.1
        bf[1] = -0.7
        bf[2] = 0.25
        element.SetValue(KM.BODY_FORCE, bf)

        divergence_stress = KM.Vector(3)
        divergence_stress[0] = 0.2
        divergence_stress[1] = -0.15
        divergence_stress[2] = 0.05
        element.SetValue(IGA.DIVERGENCE_STRESS, divergence_stress)

        return element

    def create_condition(self, model_part, integration_point, name_condition):
        # Define properties
        props = model_part.CreateNewProperties(0)
        props.SetValue(IGA.PENALTY_FACTOR, 100.0)

        # Assign constitutive law
        law = DFA.Newtonian2DLaw()
        props.SetValue(KM.CONSTITUTIVE_LAW, law)

        # Create quadrature point geometry on a curve
        geometry = TestCreationUtility.GetQuadraturePointGeometryOnCurveP2(model_part, integration_point)

        # Create the condition and add it to the model part
        condition = model_part.CreateNewCondition(name_condition, 1, geometry, props)

        return condition

    @staticmethod
    def create_condition_3d(model_part, integration_point, name_condition, properties):
        geometry = TestCreationUtility.GetQuadraturePointGeometryOnVolumeSurfaceP2(model_part, integration_point)
        condition = model_part.CreateNewCondition(name_condition, 1, geometry, properties)
        return condition

    @staticmethod
    def add_fluid_dofs_and_set_3d_state(model_part):
        for node in model_part.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.VELOCITY_Z)
            node.AddDof(KM.PRESSURE)

            velocity = node.GetSolutionStepValue(KM.VELOCITY)
            velocity[0] = 0.1 + 0.01 * node.X - 0.02 * node.Y + 0.03 * node.Z
            velocity[1] = -0.05 - 0.02 * node.X + 0.03 * node.Y + 0.01 * node.Z
            velocity[2] = 0.02 + 0.02 * node.X + 0.01 * node.Y - 0.01 * node.Z
            node.SetSolutionStepValue(KM.PRESSURE, 1.0 + 0.2 * node.X - 0.1 * node.Y + 0.05 * node.Z)

    @staticmethod
    def add_fluid_dofs_and_set_3d_state_on_geometry(geometry):
        for i in range(geometry.PointsNumber()):
            node = geometry[i]
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.VELOCITY_Z)
            node.AddDof(KM.PRESSURE)

            velocity = node.GetSolutionStepValue(KM.VELOCITY)
            velocity[0] = 0.1 + 0.01 * node.X - 0.02 * node.Y + 0.03 * node.Z
            velocity[1] = -0.05 - 0.02 * node.X + 0.03 * node.Y + 0.01 * node.Z
            velocity[2] = 0.02 + 0.02 * node.X + 0.01 * node.Y - 0.01 * node.Z
            node.SetSolutionStepValue(KM.PRESSURE, 1.0 + 0.2 * node.X - 0.1 * node.Y + 0.05 * node.Z)

    @staticmethod
    def add_fluid_dofs_on_geometry(geometry):
        for i in range(geometry.PointsNumber()):
            node = geometry[i]
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.VELOCITY_Z)
            node.AddDof(KM.PRESSURE)

    # test for stokes element
    def test_StokesElementP3(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)

        # integration point
        ipt = [0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364]

        element = self.create_element(model_part, ipt)

        # Add DOFs
        for node in model_part.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.PRESSURE)

        process_info = model_part.ProcessInfo
        element.Initialize(process_info)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, process_info)

        # expected values
        expected_LHS = [4.3418508236068837e-01, 1.8334784071599580e-01, 5.4225706107883624e-02, -2.2743271837751716e-01, -1.2041087202037005e-01, 8.0918109110538431e-03, -1.9386388862260143e-02, -1.0004830866198284e-02, 3.0187437158485690e-04, 2.6953192999768616e-02, 4.0136875065223079e-02, 2.9059468321229607e-02, -1.5258019639070167e-01, -7.3200768631259711e-02, 4.3363883978444467e-03, -1.1534410007260125e-02, -5.6851206783428559e-03, 1.6177398816363004e-04, -2.3951000078331196e-02, -2.4091394389335677e-03, 3.8932305345756304e-03, -2.4554805894705031e-02, -1.0968975189311853e-02, 5.8096588463505425e-04, -1.6987557496816854e-03, -8.0500895680252855e-04, 2.1673604742403321e-05]
        expected_RHS = [1.4052448134068332e+00,0.0000000000000000e+00,-1.6507842487568241e+00,2.0969713683772687e-01,0.0000000000000000e+00,1.1422357135115979e-01,7.8229943954274059e-03,0.0000000000000000e+00,1.7712396889478853e-02,7.5306842584076872e-01,0.0000000000000000e+00,7.8302699758045113e-01,1.1237635694157502e-01,0.0000000000000000e+00,3.1007109960153995e-01,4.1923300612959425e-03,0.0000000000000000e+00,1.8776013144483572e-02,1.0089203827470532e-01,0.0000000000000000e+00,3.2833242822625847e-01,1.5055577045423296e-02,0.0000000000000000e+00,7.4882417285054639e-02,5.6166572716448307e-04,0.0000000000000000e+00,3.7593246783978154e-03]

        tolerance = 1e-7

        for i in range(rhs.Size()):
                self.assertAlmostEqual(lhs[0, i], expected_LHS[i], delta=tolerance)

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

    def test_StokesElement3DRectangularP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)

        ipt = [0.23, 0.61, 0.37, 0.42]
        element = self.create_stokes_element_3d(model_part, ipt)

        for node in model_part.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.VELOCITY_Z)
            node.AddDof(KM.PRESSURE)

            velocity = node.GetSolutionStepValue(KM.VELOCITY)
            velocity[0] = 0.1 + 0.01 * node.X - 0.02 * node.Y
            velocity[1] = -0.05 + 0.03 * node.Y + 0.01 * node.Z
            velocity[2] = 0.02 + 0.02 * node.X - 0.01 * node.Z
            node.SetSolutionStepValue(KM.PRESSURE, 1.0 + 0.2 * node.X - 0.1 * node.Y + 0.05 * node.Z)

        process_info = model_part.ProcessInfo
        element.Initialize(process_info)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, process_info)

        geometry = element.GetGeometry()
        self.assertEqual(geometry.WorkingSpaceDimension(), 3)
        self.assertEqual(geometry.LocalSpaceDimension(), 3)
        self.assertEqual(geometry.PointsNumber(), 27)
        self.assertGreater(geometry.DomainSize(), 0.0)
        self.assertEqual(lhs.Size1(), 108)
        self.assertEqual(lhs.Size2(), 108)
        self.assertEqual(rhs.Size(), 108)

        tolerance = 1e-10
        expected_lhs = [
            8.9900319974901654e-02, 2.7661636915354355e-02, 2.7661636915354365e-02, 1.1505058297442602e-02,
           -3.1718793944655044e-02,-1.6354366987248646e-02,-1.6354366987248650e-02, 2.9226680132595053e-03,
           -5.4791976390410994e-03,-2.5235485295714630e-03,-2.5235485295714635e-03, 1.8561375559541957e-04]
        expected_rhs = [
           -1.2934143575736417e-02,-3.1761225011348601e-02,-2.1824709809219953e-02,-3.7191100322717335e-04,
            2.4258112912556050e-02,-8.0684090425895944e-03,-5.5442032199225657e-03, 7.9313375873266733e-04,
            3.2898538470421090e-03,-5.1241115900977463e-04,-3.5210307047031452e-04, 1.0674132192544911e-04]

        for i, expected_value in enumerate(expected_lhs):
            self.assertAlmostEqual(lhs[0, i], expected_value, delta=tolerance)

        for i, expected_value in enumerate(expected_rhs):
            self.assertAlmostEqual(rhs[i], expected_value, delta=tolerance)

    # test for support fluid condition (body-fitted)
    def test_SupportFluidConditionP3(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        # Add required variables
        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        process_info = model_part.ProcessInfo

        # Define integration point in parametric space
        ipt = [0.333333333333333, 0.05, 0.0, 0.086963711284364]

        # Create the support condition
        condition = self.create_condition(model_part, ipt, "SupportFluidCondition")

        # Set prescribed displacement value
        u_D = KM.Vector(3)
        u_D[0] = 0.1
        u_D[1] = -0.5
        u_D[2] = 0.0
        condition.SetValue(KM.VELOCITY_X, u_D[0])
        condition.SetValue(KM.VELOCITY_Y, u_D[1])
        mesh_size = [0.1, 0.1]
        condition.SetValue(IGA.KNOT_SPAN_SIZES, mesh_size)

        # Add DOFs
        for node in model_part.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.PRESSURE)

        # Initialize and step the condition
        condition.Initialize(process_info)
        
        # Compute local system
        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)
        
        # expected values
        expected_LHS = [1.1193281796e+02,0.0000000000e+00,2.2329836010e-18,1.1193281796e+02,0.0000000000e+00,2.2329836010e-18,2.7983204490e+01,0.0000000000e+00,5.5824590025e-19,1.1782401890e+01,0.0000000000e+00,2.3505090537e-19,1.1782401890e+01,0.0000000000e+00,2.3505090537e-19,2.9456004726e+00,0.0000000000e+00,5.8762726342e-20,3.1006320764e-01,0.0000000000e+00,6.1855501413e-21,3.1006320764e-01,0.0000000000e+00,6.1855501413e-21,7.7515801910e-02,0.0000000000e+00,1.5463875353e-21]
        expected_RHS = [2.7905688688e+01,-1.3952844344e+02,-3.4882110860e-02,2.7905688688e+01,-1.3952844344e+02,-3.4882110860e-02,6.9764221719e+00,-3.4882110860e+01,-8.7205277149e-03,2.9374409145e+00,-1.4687204572e+01,-3.6718011431e-03,2.9374409145e+00,-1.4687204572e+01,-3.6718011431e-03,7.3436022862e-01,-3.6718011431e+00,-9.1795028578e-04,7.7301076697e-02,-3.8650538349e-01,-9.6626345872e-05,7.7301076697e-02,-3.8650538349e-01,-9.6626345872e-05,1.9325269174e-02,-9.6626345872e-02,-2.4156586468e-05]

        tolerance = 1e-7

        for i in range(rhs.Size()):
                self.assertAlmostEqual(lhs[0, i], expected_LHS[i], delta=tolerance)

        # Confronto RHS
        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

    def test_SupportFluidCondition3DRectangularP2(self):
        current_model, support_model_part, condition = _create_support_condition_model_part_3d("SupportFluidCondition")
        iga_model_part = current_model.GetModelPart("IgaModelPart")

        props = support_model_part.CreateNewProperties(1)
        props.SetValue(IGA.PENALTY_FACTOR, 100.0)
        props.SetValue(KM.DYNAMIC_VISCOSITY, 0.8)
        props.SetValue(KM.CONSTITUTIVE_LAW, DFA.Newtonian3DLaw())
        for existing_condition in support_model_part.Conditions:
            existing_condition.Properties = props

        geometry = condition.GetGeometry()
        normal = geometry.Calculate(KM.NORMAL)
        self.assertAlmostEqual(normal[0], -1.0)
        self.assertAlmostEqual(normal[1], 0.0)
        self.assertAlmostEqual(normal[2], 0.0)

        condition.SetValue(KM.VELOCITY_X, 0.1)
        condition.SetValue(KM.VELOCITY_Y, -0.5)
        condition.SetValue(KM.VELOCITY_Z, 0.2)
        condition.SetValue(IGA.KNOT_SPAN_SIZES, [0.1, 0.1, 0.1])

        self.add_fluid_dofs_and_set_3d_state_on_geometry(geometry)

        process_info = iga_model_part.ProcessInfo
        condition.Initialize(process_info)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)

        self.assertEqual(lhs.Size1(), 32)
        self.assertEqual(lhs.Size2(), 32)
        self.assertEqual(rhs.Size(), 32)

        tolerance = 1e-10
        expected_lhs = [
            9.6723633543579935e+01, 3.2704174144156270e-01, 3.2704174144156270e-01, -9.6723633543579930e-02,
           -4.8310230833900063e+01,-1.6352087072078136e-01,-1.6352087072078136e-01, 4.8361816771789970e-02,
            2.5917019497006095e+01,-7.8238354270304230e-02, 8.7630570510534767e-02,-2.5917019497006095e-02,
           -1.2944687338104643e+01, 3.9119177135152120e-02,-4.3815285255267387e-02, 1.2958509748503048e-02]
        expected_rhs = [
            7.8134190774569320e+00,-3.8917145337366820e+01, 1.5552136199100953e+01, 7.7751058491018281e-03,
           -3.9108562618479867e+00, 1.9474122880381610e+01,-7.7822881842297582e+00,-3.8875529245509141e-03,
            2.1935993319305127e+00,-1.0414484331538736e+01, 4.1671823351279317e+00, 2.0833333333333337e-03,
           -1.0979107770763674e+00, 5.2114088324360348e+00,-2.0852578342306332e+00,-1.0416666666666669e-03]

        for i, expected_value in enumerate(expected_lhs):
            self.assertAlmostEqual(lhs[0, i], expected_value, delta=tolerance)

        for i, expected_value in enumerate(expected_rhs):
            self.assertAlmostEqual(rhs[i], expected_value, delta=tolerance)

    def test_SupportPressureCondition3DRectangularP2(self):
        current_model, support_model_part, condition = _create_support_condition_model_part_3d("SupportPressureCondition")
        iga_model_part = current_model.GetModelPart("IgaModelPart")

        geometry = condition.GetGeometry()
        normal = geometry.Calculate(KM.NORMAL)
        self.assertAlmostEqual(normal[0], -1.0)
        self.assertAlmostEqual(normal[1], 0.0)
        self.assertAlmostEqual(normal[2], 0.0)

        condition.SetValue(KM.PRESSURE, 2.5)
        condition.SetValue(IGA.KNOT_SPAN_SIZES, [0.1, 0.1, 0.1])

        normal_stress = KM.Vector(3)
        normal_stress[0] = 0.3
        normal_stress[1] = -0.4
        normal_stress[2] = 0.2
        condition.SetValue(KM.NORMAL_STRESS, normal_stress)

        self.add_fluid_dofs_on_geometry(geometry)

        process_info = iga_model_part.ProcessInfo
        condition.Initialize(process_info)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)

        self.assertEqual(lhs.Size1(), 32)
        self.assertEqual(lhs.Size2(), 32)
        self.assertEqual(rhs.Size(), 32)

        tolerance = 1e-12
        expected_lhs = [
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0]
        expected_rhs = [
            2.1770296377485118e-01,-3.1100423396407312e-02, 1.5550211698203656e-02, 0.0,
           -1.0885148188742559e-01, 1.5550211698203656e-02,-7.7751058491018281e-03, 0.0,
            5.8333333333333337e-02,-8.3333333333333350e-03, 4.1666666666666675e-03, 0.0,
           -2.9166666666666670e-02, 4.1666666666666675e-03,-2.0833333333333337e-03, 0.0]

        for i, expected_value in enumerate(expected_lhs):
            self.assertAlmostEqual(lhs[0, i], expected_value, delta=tolerance)

        for i, expected_value in enumerate(expected_rhs):
            self.assertAlmostEqual(rhs[i], expected_value, delta=tolerance)

    def test_NavierStokesElement3DRectangularP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        model_part.ProcessInfo.SetValue(KM.DYNAMIC_TAU, 1.0)
        model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0.1)

        ipt = [0.23, 0.61, 0.37, 0.42]
        element = self.create_navier_stokes_element_3d(model_part, ipt)
        divergence_stress = element.GetValue(IGA.DIVERGENCE_STRESS)

        self.assertAlmostEqual(divergence_stress[0], 0.2, delta=1e-10)
        self.assertAlmostEqual(divergence_stress[1], -0.15, delta=1e-10)
        self.assertAlmostEqual(divergence_stress[2], 0.05, delta=1e-10)

        for node in model_part.Nodes:
            node.AddDof(KM.VELOCITY_X)
            node.AddDof(KM.VELOCITY_Y)
            node.AddDof(KM.VELOCITY_Z)
            node.AddDof(KM.PRESSURE)

            velocity = node.GetSolutionStepValue(KM.VELOCITY)
            velocity[0] = 0.1 + 0.01 * node.X - 0.02 * node.Y
            velocity[1] = -0.05 + 0.03 * node.Y + 0.01 * node.Z
            velocity[2] = 0.02 + 0.02 * node.X - 0.01 * node.Z
            node.SetSolutionStepValue(KM.PRESSURE, 1.0 + 0.2 * node.X - 0.1 * node.Y + 0.05 * node.Z)

        process_info = model_part.ProcessInfo
        element.Initialize(process_info)
        lhs = KM.Matrix()
        rhs = KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, process_info)

        geometry = element.GetGeometry()
        self.assertEqual(geometry.WorkingSpaceDimension(), 3)
        self.assertEqual(geometry.LocalSpaceDimension(), 3)
        self.assertEqual(geometry.PointsNumber(), 27)
        self.assertGreater(geometry.DomainSize(), 0.0)
        self.assertEqual(lhs.Size1(), 108)
        self.assertEqual(lhs.Size2(), 108)
        self.assertEqual(rhs.Size(), 108)

        tolerance = 1e-10
        expected_lhs = [
            1.4914900216177548e-01, 8.691031910222818e-02, 8.691031910222821e-02, 1.1505058297442602e-02,
           -8.344190677853194e-02,-1.3032282812547660e-03,-1.3032282812547695e-03, 2.9226680132595053e-03,
           -1.3004766992038040e-02,-1.5676758924697773e-03,-1.5676758924697777e-03, 1.8561375559541957e-04]
        expected_rhs = [
           -1.2934143575736417e-02,-3.1761225011348601e-02,-2.1824709809219953e-02,-1.9290603275341875e-04,
            2.4258112912556050e-02,-8.0684090425895944e-03,-5.5442032199225657e-03, 4.1138951392215021e-04,
            3.2898538470421090e-03,-5.1241115900977463e-04,-3.5210307047031452e-04, 5.5365516924263562e-05]

        for i, expected_value in enumerate(expected_lhs):
            self.assertAlmostEqual(lhs[0, i], expected_value, delta=tolerance)

        for i, expected_value in enumerate(expected_rhs):
            self.assertAlmostEqual(rhs[i], expected_value, delta=tolerance)



if __name__ == '__main__':
    KratosUnittest.main()
