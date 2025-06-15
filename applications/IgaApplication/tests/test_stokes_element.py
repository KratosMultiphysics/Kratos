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
                

if __name__ == '__main__':
    KratosUnittest.main()