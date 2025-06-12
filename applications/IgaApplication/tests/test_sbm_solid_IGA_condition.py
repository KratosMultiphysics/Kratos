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
        geometry = TestCreationUtility.GetQuadraturePointGeometryOnCurve(model_part, degree, integration_point)

        # Create the condition and add it to the model part
        condition = model_part.CreateNewCondition("SbmSolidIGACondition", 1, geometry, props)

        node_1 = model_part.CreateNewNode(88, 0.1, 0.04, 0.0)
        
        u_D = KM.Vector(3)
        u_D[0] = 3.1
        u_D[1] = -0.5
        u_D[2] = 0.0

        node_1.SetValue(KM.DISPLACEMENT, u_D)
        condition.SetValue(IGA.PROJECTION_NODE, node_1)
        model_part.AddCondition(condition)

        return condition

    def test_SupportSolidIGAConditionP3(self):
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
            -2.18648,
            -16.9195,
            1.48656,
            -26.0148,
            3.12642,
            -13.3251,
            0.91826,
            -2.27382,
            -1.15828,
            17.3655,
            -1.48656,
            26.0148,
            -0.617851,
            12.9907,
            -0.0820705,
            2.16232
        ]
        tolerance = 1e-4

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

        expected_LHS = [
            [0, 0.596275, 0, 2.68324, 0, 2.23603, 0, 0.521741, 2.93641, 0.0313829, 4.40462, 0.141223, 2.20231, 0.117686, 0.367052, 0.0274601],
            [-0.596275, 0, 1.78883, 0, 2.23603, 0, 0.596275, 0, -0.0313829, 10.2774, 0.0941487, 15.4162, 0.117686, 7.70808, 0.0313829, 1.28468],
            [0, -1.78883, 0, -1.34055e-15, 0, 1.34162, 1.11022e-16, 0.447207, 4.40462, -0.0941487, 6.60693, -2.57848e-17, 3.30346, 0.0706116, 0.550577, 0.0235372],
            [-2.68324, 0, 1.34055e-15, 0, 2.01243, 1.77636e-15, 0.67081, 0, -0.141223, 15.4162, 1.15325e-16, 23.1243, 0.105917, 11.5621, 0.0353058, 1.92702],
            [0, -2.23603, 0, -2.01243, 0, -0.335405, 5.55112e-17, 0.0559008, 2.20231, -0.117686, 3.30346, -0.105917, 1.65173, -0.0176529, 0.275289, 0.00294215],
            [-2.23603, 0, -1.34162, -1.77636e-15, 0.335405, 0, 0.223603, 0, -0.117686, 7.70808, -0.0706116, 11.5621, 0.0176529, 5.78106, 0.0117686, 0.963511],
            [0, -0.596275, -1.11022e-16, -0.67081, -5.55112e-17, -0.223603, 0, -0.0186336, 0.367052, -0.0313829, 0.550577, -0.0353058, 0.275289, -0.0117686, 0.0458815, -0.000980716],
            [-0.521741, 0, -0.447207, 0, -0.0559008, 0, 0.0186336, 0, -0.0274601, 1.28468, -0.0235372, 1.92702, -0.00294215, 0.963511, 0.000980716, 0.160585],
            [-2.93641, 0.0313829, -4.40462, 0.141223, -2.20231, 0.117686, -0.367052, 0.0274601, 0, 0.00165173, 0, 0.0074328, 0, 0.006194, 0, 0.00144527],
            [-0.0313829, -10.2774, 0.0941487, -15.4162, 0.117686, -7.70808, 0.0313829, -1.28468, -0.00165173, 0, 0.0049552, 0, 0.006194, 0, 0.00165173, -1.38778e-17],
            [-4.40462, -0.0941487, -6.60693, -1.15325e-16, -3.30346, 0.0706116, -0.550577, 0.0235372, 0, -0.0049552, 0, -3.71343e-18, 0, 0.0037164, 0, 0.0012388],
            [-0.141223, -15.4162, 2.57848e-17, -23.1243, 0.105917, -11.5621, 0.0353058, -1.92702, -0.0074328, 0, 3.71343e-18, 0, 0.0055746, 1.11022e-16, 0.0018582, 0],
            [-2.20231, -0.117686, -3.30346, -0.105917, -1.65173, -0.0176529, -0.275289, 0.00294215, 0, -0.006194, 0, -0.0055746, 0, -0.000929099, 0, 0.00015485],
            [-0.117686, -7.70808, -0.0706116, -11.5621, 0.0176529, -5.78106, 0.0117686, -0.963511, -0.006194, 0, -0.0037164, -1.11022e-16, 0.000929099, 0, 0.0006194, -6.93889e-18],
            [-0.367052, -0.0313829, -0.550577, -0.0353058, -0.275289, -0.0117686, -0.0458815, -0.000980716, 0, -0.00165173, 0, -0.0018582, 0, -0.0006194, 0, -5.16166e-05],
            [-0.0274601, -1.28468, -0.0235372, -1.92702, -0.00294215, -0.963511, 0.000980716, -0.160585, -0.00144527, 1.38778e-17, -0.0012388, 0, -0.00015485, 6.93889e-18, 5.16166e-05, 0]
        ]
        for i in range(lhs.Size1()):
            for j in range(lhs.Size2()):
                self.assertAlmostEqual(lhs[i, j], expected_LHS[i][j], delta=tolerance)

if __name__ == '__main__':
    KratosUnittest.main()