import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest


from test_creation_utility import TestCreationUtility

class SolidIGAElementTests(KratosUnittest.TestCase):

    @staticmethod
    def create_element(model_part, integration_point):
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
    
    def create_support_condition(self, model_part, degree, integration_point):
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
    
    def create_load_condition(self, model_part, degree, integration_point):
        # Define properties
        props = model_part.CreateNewProperties(0)
        props.SetValue(KM.THICKNESS, 1.0)
        props.SetValue(KM.YOUNG_MODULUS, 100)
        props.SetValue(KM.POISSON_RATIO, 0.3)

        # Assign constitutive law
        law = SMA.LinearElasticPlaneStrain2DLaw()
        props.SetValue(KM.CONSTITUTIVE_LAW, law)

        # Create quadrature point geometry on a curve
        geometry = TestCreationUtility.GetQuadraturePointGeometryOnCurveP2(model_part, integration_point)

        # Create the condition and add it to the model part
        condition = model_part.CreateNewCondition("LoadSolidCondition", 1, geometry, props)

        force_vector = KM.Vector(3)
        force_vector[0] = -12.0
        force_vector[1] = 3.7
        condition.SetValue(KM.FORCE, force_vector)  # penalty-free formulation
        model_part.AddCondition(condition)

        return condition


    def testSolidIGAElementP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        # integration point
        ipt = [0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364]

        element = self.create_element(model_part, ipt)

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
    
    def testSupportSimmetrySolidIGAConditionP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        # Add required variables
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        process_info = model_part.ProcessInfo

        # Define integration point in parametric space
        ipt = [0.333333333333333, 0.05, 0.0, 0.086963711284364]

        # Create the support condition
        condition = self.create_support_condition(model_part, 2, ipt)

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
            -1.38611,-2.20831,0.259896,-1.59106,0.476476,-0.243452,0.158065,1.65893,0.331329,1.7239,0.126148,0.447218,0.0121589,0.0934291,0.0167184,0.0951389,0.0053195,0.0242122   ]
        tolerance = 1e-4

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

        expected_LHS = [
            (0,0.130105,0.836913,0.443947,0.418457,0.189447,0.154555,0.975368,0.242651,1.0084,0.0826867,0.26036,0.00813445,0.0509748,0.0104528,0.0518441,0.00319277,0.0131784),(-0.130105,0,1.35774,0.557942,0.711397,0.278971,0.261068,1.70964,0.417684,1.76837,0.143575,0.456776,0.0141008,0.0899811,0.0182223,0.0915266,0.00558594,0.023268),(-0.836913,-1.35774,0,-1.0439,0.209228,-0.182514,0.0664585,0.818753,0.154555,0.851789,0.0606627,0.221206,0.00581613,0.0468533,0.00813445,0.0477227,0.00261319,0.012148),(-0.443947,-0.557942,1.0439,0,0.632936,0.139486,0.228032,1.65091,0.384648,1.70964,0.135316,0.442093,0.0132315,0.0884356,0.0173529,0.0899811,0.0053686,0.0228817),(-0.418457,-0.711397,-0.209228,-0.632936,0,-0.138619,-0.00540941,0.165534,0.0166146,0.173793,0.00965966,0.0455131,0.000874454,0.010683,0.00145403,0.0109003,0.000508403,0.00277941),(-0.189447,-0.278971,0.182514,-0.139486,0.138619,0,0.0487491,0.398045,0.0879029,0.412728,0.0317642,0.106853,0.00309053,0.0217225,0.00412089,0.0221089,0.00128781,0.00562382),(-0.154555,-0.261068,-0.0664585,-0.228032,0.00540941,-0.0487491,0,0.0737478,0.00927328,0.0772253,0.00463664,0.0201757,0.000428129,0.00460465,0.000672163,0.00469616,0.000229049,0.00119692),(-0.975368,-1.70964,-0.818753,-1.65091,-0.165534,-0.398045,-0.0737478,0,-0.057262,0.00618218,-0.010194,0.00309109,-0.00117961,0.00473585,-0.000745777,0.00489854,-7.79849e-05,0.00126531),(-0.242651,-0.417684,-0.154555,-0.384648,-0.0166146,-0.0879029,-0.00927328,0.057262,0,0.0607395,0.00231832,0.0160542,0.000184095,0.00417081,0.000428129,0.00426232,0.000168041,0.00108846),(-1.0084,-1.76837,-0.851789,-1.70964,-0.173793,-0.412728,-0.0772253,-0.00618218,-0.0607395,0,-0.0110634,0.00154555,-0.00127113,0.00457316,-0.00083729,0.00473585,-0.000100863,0.00122463),(-0.0826867,-0.143575,-0.0606627,-0.135316,-0.00965966,-0.0317642,-0.00463664,0.010194,-0.00231832,0.0110634,0,0.0029832,-1.49845e-05,0.000934243,4.60239e-05,0.000957121,2.67581e-05,0.000245),(-0.26036,-0.456776,-0.221206,-0.442093,-0.0455131,-0.106853,-0.0201757,-0.00309109,-0.0160542,-0.00154555,-0.0029832,0,-0.00034066,0.00110262,-0.000232201,0.00114329,-3.09353e-05,0.00029599),(-0.00813445,-0.0141008,-0.00581613,-0.0132315,-0.000874454,-0.00309053,-0.000428129,0.00117961,-0.000184095,0.00127113,1.49845e-05,0.00034066,0,0.000101145,6.42194e-06,0.000103554,3.21097e-06,2.64905e-05),(-0.0509748,-0.0899811,-0.0468533,-0.0884356,-0.010683,-0.0217225,-0.00460465,-0.00473585,-0.00417081,-0.00457316,-0.000934243,-0.00110262,-0.000101145,0,-8.97287e-05,4.28129e-06,-1.9578e-05,2.14065e-06),(-0.0104528,-0.0182223,-0.00813445,-0.0173529,-0.00145403,-0.00412089,-0.000672163,0.000745777,-0.000428129,0.00083729,-4.60239e-05,0.000232201,-6.42194e-06,8.97287e-05,0,9.21369e-05,1.60548e-06,2.36363e-05),(-0.0518441,-0.0915266,-0.0477227,-0.0899811,-0.0109003,-0.0221089,-0.00469616,-0.00489854,-0.00426232,-0.00473585,-0.000957121,-0.00114329,-0.000103554,-4.28129e-06,-9.21369e-05,0,-2.018e-05,1.07032e-06),(-0.00319277,-0.00558594,-0.00261319,-0.0053686,-0.000508403,-0.00128781,-0.000229049,7.79849e-05,-0.000168041,0.000100863,-2.67581e-05,3.09353e-05,-3.21097e-06,1.9578e-05,-1.60548e-06,2.018e-05,0,5.19552e-06),(-0.0131784,-0.023268,-0.012148,-0.0228817,-0.00277941,-0.00562382,-0.00119692,-0.00126531,-0.00108846,-0.00122463,-0.000245,-0.00029599,-2.64905e-05,-2.14065e-06,-2.36363e-05,-1.07032e-06,-5.19552e-06,0)        ]
        for i in range(lhs.Size1()):
            for j in range(lhs.Size2()):
                self.assertAlmostEqual(lhs[i, j], expected_LHS[i][j], delta=tolerance)
    
    def testSupportSolidIGAConditionP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        # Add required variables
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        process_info = model_part.ProcessInfo

        # Define integration point in parametric space
        ipt = [0.333333333333333, 0.05, 0.0, 0.086963711284364]

        # Create the support condition
        condition = self.create_support_condition(model_part, 3, ipt)

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

    def testLoadSolidIGAConditionP2(self):
        model = KM.Model()
        model_part = model.CreateModelPart("ModelPart")
        model_part.SetBufferSize(2)

        # Add required variables
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        process_info = model_part.ProcessInfo

        # Define integration point in parametric space
        ipt = [0.333333333333333, 0.05, 0.0, 0.086963711284364]

        # Create the support condition
        condition = self.create_load_condition(model_part, 3, ipt)

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
            -0.837171,0.258128,-0.837171,0.258128,-0.209293,0.0645319,-0.0881232,0.0271713,-0.0881232,0.0271713,-0.0220308,0.00679283,-0.00231903,0.000715035,-0.00231903,0.000715035,-0.000579758,0.000178759
        ]
        tolerance = 1e-4

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)
    
if __name__ == '__main__':
    KratosUnittest.main()
                
