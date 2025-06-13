import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_creation_utility import TestCreationUtility

class TestSupportSimmetrySolidIGACondition(KratosUnittest.TestCase):

    def create_condition(self, model_part, degree, integration_point):
        # Define properties
        props = model_part.CreateNewProperties(0)
        props.SetValue(IGA.PENALTY_FACTOR, -1.0)  # penalty-free formulation
        props.SetValue(KM.THICKNESS, 1.0)
        props.SetValue(KM.YOUNG_MODULUS, 100)
        props.SetValue(KM.POISSON_RATIO, 0.3)

        # Assign constitutive law
        law = SMA.LinearElasticPlaneStrain2DLaw()
        props.SetValue(KM.CONSTITUTIVE_LAW, law)

        # Create quadrature point geometry on a curve
        geometry = TestCreationUtility.GetQuadraturePointGeometryOnCurveP2(model_part, integration_point)

        # Create the condition and add it to the model part
        condition = model_part.CreateNewCondition("SupportSolidCondition", 1, geometry, props)
        model_part.AddCondition(condition)

        return condition

    def test_SupportSimmetrySolidIGAConditionP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        # Add required variables
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        process_info = model_part.ProcessInfo

        # Define integration point in parametric space
        ipt = [0.333333333333333, 0.05, 0.0, 0.086963711284364]

        # Create the support condition
        condition = self.create_condition(model_part, 2, ipt)

        # Set prescribed displacement value
        u_D = KM.Vector(3)
        u_D[0] = 0.1
        u_D[1] = -0.5
        u_D[2] = 0.0
        condition.SetValue(KM.DISPLACEMENT, u_D)

        direction = KM.Vector(3)
        direction[0] = 0.36
        direction[1] = 0.64
        direction[2] = 0.0
        condition.SetValue(KM.DIRECTION, direction)

        # Add degrees of freedom
        for node in model_part.Nodes:
            node.AddDof(KM.DISPLACEMENT_X)
            node.AddDof(KM.DISPLACEMENT_Y)

        # Initialize and step the condition
        condition.Initialize(process_info)
        
        # Compute local system
        lhs = KM.Matrix()
        rhs = KM.Vector()
        condition.CalculateLocalSystem(lhs, rhs, process_info)

        # Check the RHS
        expected_RHS = [
            2.48866,3.96487,-0.466624,2.85663,-0.855478,0.4371,-0.283795,-2.97848,-0.594878,-3.09514,-0.22649,-0.802949,-0.0218304,-0.167745,-0.0300168,-0.170815,-0.00955079,-0.0434713
        ]
        tolerance = 1e-4

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

        expected_LHS = [
            (0,0.130105,0.836913,0.443947,0.418457,0.189447,0.154555,0.975368,0.242651,1.0084,0.0826867,0.26036,0.00813445,0.0509748,0.0104528,0.0518441,0.00319277,0.0131784),(-0.130105,0,1.35774,0.557942,0.711397,0.278971,0.261068,1.70964,0.417684,1.76837,0.143575,0.456776,0.0141008,0.0899811,0.0182223,0.0915266,0.00558594,0.023268),(-0.836913,-1.35774,0,-1.0439,0.209228,-0.182514,0.0664585,0.818753,0.154555,0.851789,0.0606627,0.221206,0.00581613,0.0468533,0.00813445,0.0477227,0.00261319,0.012148),(-0.443947,-0.557942,1.0439,0,0.632936,0.139486,0.228032,1.65091,0.384648,1.70964,0.135316,0.442093,0.0132315,0.0884356,0.0173529,0.0899811,0.0053686,0.0228817),(-0.418457,-0.711397,-0.209228,-0.632936,0,-0.138619,-0.00540941,0.165534,0.0166146,0.173793,0.00965966,0.0455131,0.000874454,0.010683,0.00145403,0.0109003,0.000508403,0.00277941),(-0.189447,-0.278971,0.182514,-0.139486,0.138619,0,0.0487491,0.398045,0.0879029,0.412728,0.0317642,0.106853,0.00309053,0.0217225,0.00412089,0.0221089,0.00128781,0.00562382),(-0.154555,-0.261068,-0.0664585,-0.228032,0.00540941,-0.0487491,0,0.0737478,0.00927328,0.0772253,0.00463664,0.0201757,0.000428129,0.00460465,0.000672163,0.00469616,0.000229049,0.00119692),(-0.975368,-1.70964,-0.818753,-1.65091,-0.165534,-0.398045,-0.0737478,0,-0.057262,0.00618218,-0.010194,0.00309109,-0.00117961,0.00473585,-0.000745777,0.00489854,-7.79849e-05,0.00126531),(-0.242651,-0.417684,-0.154555,-0.384648,-0.0166146,-0.0879029,-0.00927328,0.057262,0,0.0607395,0.00231832,0.0160542,0.000184095,0.00417081,0.000428129,0.00426232,0.000168041,0.00108846),(-1.0084,-1.76837,-0.851789,-1.70964,-0.173793,-0.412728,-0.0772253,-0.00618218,-0.0607395,0,-0.0110634,0.00154555,-0.00127113,0.00457316,-0.00083729,0.00473585,-0.000100863,0.00122463),(-0.0826867,-0.143575,-0.0606627,-0.135316,-0.00965966,-0.0317642,-0.00463664,0.010194,-0.00231832,0.0110634,0,0.0029832,-1.49845e-05,0.000934243,4.60239e-05,0.000957121,2.67581e-05,0.000245),(-0.26036,-0.456776,-0.221206,-0.442093,-0.0455131,-0.106853,-0.0201757,-0.00309109,-0.0160542,-0.00154555,-0.0029832,0,-0.00034066,0.00110262,-0.000232201,0.00114329,-3.09353e-05,0.00029599),(-0.00813445,-0.0141008,-0.00581613,-0.0132315,-0.000874454,-0.00309053,-0.000428129,0.00117961,-0.000184095,0.00127113,1.49845e-05,0.00034066,0,0.000101145,6.42194e-06,0.000103554,3.21097e-06,2.64905e-05),(-0.0509748,-0.0899811,-0.0468533,-0.0884356,-0.010683,-0.0217225,-0.00460465,-0.00473585,-0.00417081,-0.00457316,-0.000934243,-0.00110262,-0.000101145,0,-8.97287e-05,4.28129e-06,-1.9578e-05,2.14065e-06),(-0.0104528,-0.0182223,-0.00813445,-0.0173529,-0.00145403,-0.00412089,-0.000672163,0.000745777,-0.000428129,0.00083729,-4.60239e-05,0.000232201,-6.42194e-06,8.97287e-05,0,9.21369e-05,1.60548e-06,2.36363e-05),(-0.0518441,-0.0915266,-0.0477227,-0.0899811,-0.0109003,-0.0221089,-0.00469616,-0.00489854,-0.00426232,-0.00473585,-0.000957121,-0.00114329,-0.000103554,-4.28129e-06,-9.21369e-05,0,-2.018e-05,1.07032e-06),(-0.00319277,-0.00558594,-0.00261319,-0.0053686,-0.000508403,-0.00128781,-0.000229049,7.79849e-05,-0.000168041,0.000100863,-2.67581e-05,3.09353e-05,-3.21097e-06,1.9578e-05,-1.60548e-06,2.018e-05,0,5.19552e-06),(-0.0131784,-0.023268,-0.012148,-0.0228817,-0.00277941,-0.00562382,-0.00119692,-0.00126531,-0.00108846,-0.00122463,-0.000245,-0.00029599,-2.64905e-05,-2.14065e-06,-2.36363e-05,-1.07032e-06,-5.19552e-06,0)
        ]
        for i in range(lhs.Size1()):
            for j in range(lhs.Size2()):
                self.assertAlmostEqual(lhs[i, j], expected_LHS[i][j], delta=tolerance)

if __name__ == '__main__':
    KratosUnittest.main()