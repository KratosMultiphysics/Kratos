import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as DFA
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.KratosUnittest as KratosUnittest
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
        divergence_stress = [0.1, 100.0, -0.43]
        element.SetValue(IGA.DIVERGENCE_STRESS, divergence_stress)

        model_part.AddElement(element)
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
        model_part.AddCondition(condition)

        return condition

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
        condition.SetValue(KM.VELOCITY, u_D)
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


if __name__ == '__main__':
    KratosUnittest.main()