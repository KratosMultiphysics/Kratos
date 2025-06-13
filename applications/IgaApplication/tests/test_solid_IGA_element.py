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
        geometry = TestCreationUtility.GetQuadraturePointGeometryP2(model_part, integration_point)
        element = model_part.CreateNewElement("SolidElement", 1, geometry, props)

        bf = KM.Vector(3)
        bf[0] = 30.0
        bf[1] = 0.0
        bf[2] = 0.0
        element.SetValue(KM.BODY_FORCE, bf)

        model_part.AddElement(element)
        return element


    def test_SolidIGAElementP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        # integration point
        ipt = [0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364]

        element = self.create_element(model_part, 2, ipt)

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
        expected_LHS = [
        (21.9289,13.2222,-13.5867,-3.71042,-1.13582,-0.350451,3.83914,-2.97329,-8.46185,-3.48946,-0.652731,-0.243805,-0.545731,-1.746,-1.29186,-0.668603,-0.0933508,-0.040166),(13.2222,26.3239,-6.55217,-0.888702,-0.56248,-0.212853,0.379725,-13.5869,-4.512,-4.60886,-0.338765,-0.268239,-0.847564,-5.53058,-0.738562,-1.17113,-0.0503875,-0.0565924),(-13.5867,-6.55217,13.5737,-1.82586,1.0884,-0.0997559,-8.46185,5.79722,7.09793,0.410584,0.5767,-0.00163853,-1.29186,2.02379,0.927336,0.241107,0.0763825,0.00672311),(-3.71042,-0.888702,-1.82586,4.32486,-0.115576,0.327635,4.21727,-4.60886,-0.0524364,1.701,-0.0273899,0.152573,1.39641,-1.17113,0.117041,0.14527,0.000958864,0.0173587),(-1.13582,-0.56248,1.0884,-0.115576,0.0875313,-0.00549206,-0.652731,0.449097,0.5767,0.0500604,0.0466627,0.00123501,-0.0933508,0.160719,0.0763825,0.0217116,0.00621876,0.000725231),(-0.350451,-0.212853,-0.0997559,0.327635,-0.00549206,0.0256306,0.312546,-0.268239,0.0212059,0.152573,-0.000157725,0.0128771,0.108908,-0.0565924,0.0128443,0.0173587,0.00035205,0.00161022),(3.83914,0.379725,-8.46185,4.21727,-0.652731,0.312546,7.84977,-5.18713,-3.67033,1.45561,-0.317551,0.137484,1.8277,-1.41715,-0.375928,0.0872443,-0.0382236,0.0143988),(-2.97329,-13.5869,5.79722,-4.60886,0.449097,-0.268239,-5.18713,12.9921,2.57045,0.555397,0.220664,-0.0308875,-1.17642,4.45672,0.272527,0.47972,0.026883,0.0109824),(-8.46185,-4.512,7.09793,-0.0524364,0.5767,0.0212059,-3.67033,2.57045,3.93275,0.716295,0.313864,0.0391347,-0.375928,1.0127,0.54417,0.195695,0.0426945,0.0089636),(-3.48946,-4.60886,0.410584,1.701,0.0500604,0.152573,1.45561,0.555397,0.716295,1.36301,0.045341,0.0986054,0.640563,0.47972,0.162452,0.243091,0.0085549,0.015467),(-0.652731,-0.338765,0.5767,-0.0273899,0.0466627,-0.000157725,-0.317551,0.220664,0.313864,0.045341,0.025186,0.00215456,-0.0382236,0.0834489,0.0426945,0.0141156,0.00339833,0.000588637),(-0.243805,-0.268239,-0.00163853,0.152573,0.00123501,0.0128771,0.137484,-0.0308875,0.0391347,0.0986054,0.00215456,0.00752913,0.054343,0.0109824,0.0106038,0.015467,0.000488644,0.00109289),(-0.545731,-0.847564,-1.29186,1.39641,-0.0933508,0.108908,1.8277,-1.17642,-0.375928,0.640563,-0.0382236,0.054343,0.528913,-0.254367,-0.00797792,0.0713806,-0.00353971,0.00674193),(-1.746,-5.53058,2.02379,-1.17113,0.160719,-0.0565924,-1.41715,4.45672,1.0127,0.47972,0.0834489,0.0109824,-0.254367,1.59125,0.12605,0.212624,0.0108209,0.00700588),(-1.29186,-0.738562,0.927336,0.117041,0.0763825,0.0128443,-0.375928,0.272527,0.54417,0.162452,0.0426945,0.0106038,-0.00797792,0.12605,0.0792301,0.0351257,0.00595595,0.00191909),(-0.668603,-1.17113,0.241107,0.14527,0.0217116,0.0173587,0.0872443,0.47972,0.195695,0.243091,0.0141156,0.015467,0.0713806,0.212624,0.0351257,0.0547061,0.00222344,0.00289807),(-0.0933508,-0.0503875,0.0763825,0.000958864,0.00621876,0.00035205,-0.0382236,0.026883,0.0426945,0.0085549,0.00339833,0.000488644,-0.00353971,0.0108209,0.00595595,0.00222344,0.000464093,0.000105656),(-0.040166,-0.0565924,0.00672311,0.0173587,0.000725231,0.00161022,0.0143988,0.0109824,0.0089636,0.015467,0.000588637,0.00109289,0.00674193,0.00700588,0.00191909,0.00289807,0.000105656,0.00017723)
        ]

        expected_RHS = [
        5.62098,0,0.838789,0,0.031292,0,3.01227,0,0.449505,0,0.0167693,0,0.403568,0,0.0602223,0,0.00224666,0] 

        tolerance = 1e-4

        # Confront LHS
        self.assertEqual(lhs.Size1(), len(expected_LHS))
        self.assertEqual(lhs.Size2(), len(expected_LHS[0]))

        for i in range(lhs.Size1()):
            for j in range(lhs.Size2()):
                self.assertAlmostEqual(lhs[i, j], expected_LHS[i][j], delta=tolerance)

        # Confronto RHS
        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)
                

if __name__ == '__main__':
    KratosUnittest.main()