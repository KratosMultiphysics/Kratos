import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_creation_utility import TestCreationUtility

class TestSupportSolidIGACondition(KratosUnittest.TestCase):

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

    def test_SupportSolidIGAConditionP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        # Add required variables
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        process_info = model_part.ProcessInfo

        # Define integration point in parametric space
        ipt = [0.333333333333333, 0.05, 0.0, 0.086963711284364]

        # Create the support condition
        condition = self.create_condition(model_part, 3, ipt)

        # Set prescribed displacement value
        u_D = KM.Vector(3)
        u_D[0] = 0.1
        u_D[1] = -0.5
        u_D[2] = 0.0
        condition.SetValue(KM.DISPLACEMENT, u_D)

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
            -2.7362,-4.54032,1.79177,-5.14405,1.57993,-1.43695,-0.585333,4.72503,-0.108705,4.66148,0.0919808,1.14948,-0.0232275,0.261263,-0.0106846,0.25959,0.00046455,0.0644795
        ]
        tolerance = 1e-4

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

        expected_LHS = [
            (0,0.807208,2.22045e-16,3.22883,1.66533e-16,1.41261,1.19255,0.0849692,1.19255,0.339877,0.298138,0.148696,0.0627658,0.00223603,0.0627658,0.00894413,0.0156915,0.00391306),(-0.807208,0,2.82523,8.88178e-16,1.61442,2.22045e-16,-0.0849692,4.17393,0.297392,4.17393,0.169938,1.04348,-0.00223603,0.21968,0.00782611,0.21968,0.00447207,0.0549201),(-2.22045e-16,-2.82523,0,-0.403604,5.55112e-17,0.504505,1.19255,-0.297392,1.19255,-0.0424846,0.298138,0.0531058,0.0627658,-0.00782611,0.0627658,-0.00111802,0.0156915,0.00139752),(-3.22883,-8.88178e-16,0.403604,0,1.00901,1.11022e-16,-0.339877,4.17393,0.0424846,4.17393,0.106212,1.04348,-0.00894413,0.21968,0.00111802,0.21968,0.00279504,0.0549201),(-1.66533e-16,-1.61442,-5.55112e-17,-1.00901,0,-0.100901,0.298138,-0.169938,0.298138,-0.106212,0.0745344,-0.0106212,0.0156915,-0.00447207,0.0156915,-0.00279504,0.00392286,-0.000279504),(-1.41261,-2.22045e-16,-0.504505,-1.11022e-16,0.100901,0,-0.148696,1.04348,-0.0531058,1.04348,0.0106212,0.26087,-0.00391306,0.0549201,-0.00139752,0.0549201,0.000279504,0.01373),(-1.19255,0.0849692,-1.19255,0.339877,-0.298138,0.148696,0,0.00894413,0,0.0357765,-3.46945e-18,0.0156522,0.00330346,0.000235372,0.00330346,0.000941487,0.000825866,0.000411901),(-0.0849692,-4.17393,0.297392,-4.17393,0.169938,-1.04348,-0.00894413,0,0.0313045,-5.55112e-17,0.0178883,0,-0.000235372,0.0115621,0.000823802,0.0115621,0.000470744,0.00289053),(-1.19255,-0.297392,-1.19255,-0.0424846,-0.298138,0.0531058,0,-0.0313045,0,-0.00447207,0,0.00559008,0.00330346,-0.000823802,0.00330346,-0.000117686,0.000825866,0.000147107),(-0.339877,-4.17393,0.0424846,-4.17393,0.106212,-1.04348,-0.0357765,5.55112e-17,0.00447207,0,0.0111802,0,-0.000941487,0.0115621,0.000117686,0.0115621,0.000294215,0.00289053),(-0.298138,-0.169938,-0.298138,-0.106212,-0.0745344,-0.0106212,3.46945e-18,-0.0178883,0,-0.0111802,0,-0.00111802,0.000825866,-0.000470744,0.000825866,-0.000294215,0.000206467,-2.94215e-05),(-0.148696,-1.04348,-0.0531058,-1.04348,0.0106212,-0.26087,-0.0156522,0,-0.00559008,0,0.00111802,0,-0.000411901,0.00289053,-0.000147107,0.00289053,2.94215e-05,0.000722633),(-0.0627658,0.00223603,-0.0627658,0.00894413,-0.0156915,0.00391306,-0.00330346,0.000235372,-0.00330346,0.000941487,-0.000825866,0.000411901,0,6.194e-06,2.71051e-20,2.4776e-05,0,1.08395e-05),(-0.00223603,-0.21968,0.00782611,-0.21968,0.00447207,-0.0549201,-0.000235372,-0.0115621,0.000823802,-0.0115621,0.000470744,-0.00289053,-6.194e-06,0,2.1679e-05,1.0842e-19,1.2388e-05,0),(-0.0627658,-0.00782611,-0.0627658,-0.00111802,-0.0156915,0.00139752,-0.00330346,-0.000823802,-0.00330346,-0.000117686,-0.000825866,0.000147107,-2.71051e-20,-2.1679e-05,0,-3.097e-06,-6.77626e-21,3.87125e-06),(-0.00894413,-0.21968,0.00111802,-0.21968,0.00279504,-0.0549201,-0.000941487,-0.0115621,0.000117686,-0.0115621,0.000294215,-0.00289053,-2.4776e-05,-1.0842e-19,3.097e-06,0,7.7425e-06,-2.71051e-20),(-0.0156915,-0.00447207,-0.0156915,-0.00279504,-0.00392286,-0.000279504,-0.000825866,-0.000470744,-0.000825866,-0.000294215,-0.000206467,-2.94215e-05,0,-1.2388e-05,6.77626e-21,-7.7425e-06,0,-7.7425e-07),(-0.00391306,-0.0549201,-0.00139752,-0.0549201,0.000279504,-0.01373,-0.000411901,-0.00289053,-0.000147107,-0.00289053,2.94215e-05,-0.000722633,-1.08395e-05,0,-3.87125e-06,2.71051e-20,7.7425e-07,0)
        ]
        for i in range(lhs.Size1()):
            for j in range(lhs.Size2()):
                self.assertAlmostEqual(lhs[i, j], expected_LHS[i][j], delta=tolerance)

if __name__ == '__main__':
    KratosUnittest.main()