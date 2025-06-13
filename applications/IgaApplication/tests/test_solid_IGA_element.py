import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest


from test_creation_utility import TestCreationUtility

class SolidIGAElementTests(KratosUnittest.TestCase):

    @staticmethod
    def create_element(model_part, polynomial_degree, integration_point):
        # Proprietà
        props = model_part.GetProperties()[1]
        props.SetValue(KM.THICKNESS, 1.0)
        props.SetValue(KM.YOUNG_MODULUS, 100)
        props.SetValue(KM.POISSON_RATIO, 0.3)

        law = SMA.LinearElasticPlaneStrain2DLaw()
        props.SetValue(KM.CONSTITUTIVE_LAW, law)
        # Geometria del punto di quadratura
        geometry = TestCreationUtility.GetQuadraturePointGeometry(model_part, polynomial_degree, integration_point)
        element = model_part.CreateNewElement("SolidElement", 1, geometry, props)

        bf = KM.Vector(3)
        bf[0] = 30.0
        bf[1] = 0.0
        bf[2] = 0.0
        element.SetValue(KM.BODY_FORCE, bf)

        model_part.AddElement(element)
        return element


    def test_SolidIGAElementP3(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        # integration point
        ipt = [0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364]

        element = self.create_element(model_part, 3, ipt)

        # Add DOFs
        for node in model_part.Nodes:
            node.AddDof(KM.DISPLACEMENT_X)
            node.AddDof(KM.DISPLACEMENT_Y)

        process_info = model_part.ProcessInfo
        element.Initialize(process_info)

        lhs = KM.Matrix()
        rhs = KM.Vector()
        element.CalculateLocalSystem(lhs, rhs, process_info)

        # expected values
        expected_LHS = [[26.634, 13.8059, 0.680681, -2.84413, -0.34324, -0.654986, -0.0183364, -0.0273023, -20.4029, -6.80383, -5.98196, -3.11307, -0.551908, -0.350916, -0.0163522, -0.0116783],
        [13.8059, 77.423, -5.81133, 15.8213, -1.09777, 1.06788, -0.0438207, 0.0237591, -3.30279, -75.6427, -3.12446, -17.3359, -0.411087, -1.32364, -0.0146501, -0.0336702],
        [0.680681, -5.81133, 4.64531, -2.62913, 0.681826, -0.295276, 0.0252949, -0.00980869, -5.98196, 7.37865, -0.135104, 1.29569, 0.0797439, 0.0701177, 0.00421729, 0.0010834],
        [-2.84413, 15.8213, -2.62913, 4.82509, -0.344831, 0.455791, -0.0122736, 0.013718, 5.19511, -17.3359, 0.628968, -3.53646, 0.00709391, -0.238199, -0.000814292, -0.00528591],
        [-0.34324, -1.09777, 0.681826, -0.344831, 0.107478, -0.0331235, 0.00408086, -0.00100772, -0.551908, 1.21471, 0.0797439, 0.24534, 0.0211172, 0.0163239, 0.000902422, 0.000356709],
        [-0.654986, 1.06788, -0.295276, 0.455791, -0.0331235, 0.0501805, -0.00109968, 0.00165026, 0.830399, -1.32364, 0.146039, -0.238199, 0.00792415, -0.0134391, 0.00012316, -0.000226463],
        [-0.0183364, -0.0438207, 0.0252949, -0.0122736, 0.00408086, -0.00109968, 0.00015605, -3.1924e-05, -0.0163522, 0.0467291, 0.00421729, 0.00979922, 0.000902422, 0.000681862, 3.7062e-05, 1.57328e-05],
        [-0.0273023, 0.0237591, -0.00980869, 0.013718, -0.00100772, 0.00165026, -3.1924e-05, 5.66304e-05, 0.0316649, -0.0336702, 0.00609706, -0.00528591, 0.000380995, -0.000226463, 7.6372e-06, -1.45571e-06],
        [-20.4029, -3.30279, -5.98196, 5.19511, -0.551908, 0.830399, -0.0163522, 0.0316649, 22.0725, -3.69928, 4.56149, 0.762082, 0.312053, 0.175503, 0.00705742, 0.00731562],
        [-6.80383, -75.6427, 7.37865, -17.3359, 1.21471, -1.32364, 0.0467291, -0.0336702, -3.69928, 76.1198, 1.55714, 16.9301, 0.294145, 1.25511, 0.0117417, 0.0310145],
        [-5.98196, -3.12446, -0.135104, 0.628968, 0.0797439, 0.146039, 0.00421729, 0.00609706, 4.56149, 1.55714, 1.34361, 0.704473, 0.124318, 0.0791189, 0.00369049, 0.00262823],
        [-3.11307, -17.3359, 1.29569, -3.53646, 0.24534, -0.238199, 0.00979922, -0.00528591, 0.762082, 16.9301, 0.704473, 3.88174, 0.0923972, 0.296503, 0.00328871, 0.00754528],
        [-0.551908, -0.411087, 0.0797439, 0.00709391, 0.0211172, 0.00792415, 0.000902422, 0.000380995, 0.312053, 0.294145, 0.124318, 0.0923972, 0.0133397, 0.00887541, 0.000432846, 0.000270018],
        [-0.350916, -1.32364, 0.0701177, -0.238199, 0.0163239, -0.0134391, 0.000681862, -0.000226463, 0.175503, 1.25511, 0.0791189, 0.296503, 0.00887541, 0.0232839, 0.000294658, 0.000607968],
        [-0.0163522, -0.0146501, 0.00421729, -0.000814292, 0.000902422, 0.00012316, 3.7062e-05, 7.6372e-06, 0.00705742, 0.0117417, 0.00369049, 0.00328871, 0.000432846, 0.000294658, 1.46821e-05, 8.55402e-06],
        [-0.0116783, -0.0336702, 0.0010834, -0.00528591, 0.000356709, -0.000226463, 1.57328e-05, -1.45571e-06, 0.00731562, 0.0310145, 0.00262823, 0.00754528, 0.000270018, 0.000607968, 8.55402e-06, 1.62397e-05]]  # ← inserisci qui la matrice
        expected_RHS = [0.165807,
        0,
        0.0371137,
        0,
        0.00276914,
        0,
        6.88706e-05,
        0,
        0.0444278,
        0,
        0.00994458,
        0,
        0.000741988,
        0,
        1.84538e-05,
        0] 

        tolerance = 1e-4

        # Compare LHS
        self.assertEqual(lhs.Size1(), len(expected_LHS))
        self.assertEqual(lhs.Size2(), len(expected_LHS[0]))

        for i in range(lhs.Size1()):
            for j in range(lhs.Size2()):
                self.assertAlmostEqual(lhs[i, j], expected_LHS[i][j], delta=tolerance)

        # Compare RHS
        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)
                

if __name__ == '__main__':
    KratosUnittest.main()