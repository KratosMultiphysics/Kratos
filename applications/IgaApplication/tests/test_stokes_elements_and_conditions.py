import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as DFA
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.KratosUnittest as KratosUnittest
try:
    from .test_creation_utility import TestCreationUtility
except ImportError:
    from test_creation_utility import TestCreationUtility

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

        geometry = TestCreationUtility.GetQuadraturePointGeometryFromRectangularVolumeP2(model_part, integration_point)
        element = model_part.CreateNewElement("NavierStokesElement", 1, geometry, props)

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
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)

        props = model_part.CreateNewProperties(1)
        props.SetValue(IGA.PENALTY_FACTOR, 100.0)
        props.SetValue(KM.DYNAMIC_VISCOSITY, 0.8)
        props.SetValue(KM.CONSTITUTIVE_LAW, DFA.Newtonian3DLaw())

        ipt = [0.23, 0.61, 0.0, 0.42]
        condition = self.create_condition_3d(model_part, ipt, "SupportFluidCondition", props)
        condition.SetValue(KM.VELOCITY_X, 0.1)
        condition.SetValue(KM.VELOCITY_Y, -0.5)
        condition.SetValue(KM.VELOCITY_Z, 0.2)
        condition.SetValue(IGA.KNOT_SPAN_SIZES, [0.1, 0.1, 0.1])

        self.add_fluid_dofs_and_set_3d_state(model_part)

        process_info = model_part.ProcessInfo
        condition.Initialize(process_info)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)

        self.assertEqual(lhs.Size1(), 108)
        self.assertEqual(lhs.Size2(), 108)
        self.assertEqual(rhs.Size(), 108)

        tolerance = 1e-10
        expected_lhs = [
            2.6684597075089113e-02, 0.0, -4.6207094502318810e-05, 0.0,
            8.5252089356778207e-02, 0.0, -7.5611609185612600e-05, 0.0,
            6.8090954486257920e-02, 0.0, -2.8757662120687565e-06, 0.0,
            2.2031898097894076e-01, 0.0, -3.8150472896786286e-04, 0.0]
        expected_rhs = [
            4.7297899548900030e-01, -2.3694691470450016e+00, 9.4441149564949760e-01, 2.3672273625000015e-04,
            1.5127792718220012e+00, -7.5699923399100050e+00, 3.0159331123347590e+00, 7.5628302750000050e-04,
            1.2096193726890008e+00, -6.0461627130450030e+00, 2.4078092408907485e+00, 6.0404423625000030e-04,
            3.9051086294220010e+00, -1.9556668887910003e+01, 7.8223505589522590e+00, 1.9544800275000003e-03]

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
