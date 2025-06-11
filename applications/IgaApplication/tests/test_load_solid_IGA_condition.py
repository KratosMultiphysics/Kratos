import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

from test_creation_utility import TestCreationUtility

class TestLoadSolidIGACondition(KratosUnittest.TestCase):

    def create_condition(self, model_part, degree, integration_point):
        # Define properties
        props = model_part.CreateNewProperties(0)
        props.SetValue(KM.THICKNESS, 1.0)
        props.SetValue(KM.YOUNG_MODULUS, 100)
        props.SetValue(KM.POISSON_RATIO, 0.3)

        # Assign constitutive law
        law = SMA.LinearElasticPlaneStrain2DLaw()
        props.SetValue(KM.CONSTITUTIVE_LAW, law)

        # Create quadrature point geometry on a curve
        geometry = TestCreationUtility.GetQuadraturePointGeometryOnCurve(model_part, degree, integration_point)

        # Create the condition and add it to the model part
        condition = model_part.CreateNewCondition("LoadSolidIGACondition", 1, geometry, props)

        force_vector = KM.Vector(3)
        force_vector[0] = -12.0
        force_vector[1] = 3.7
        condition.SetValue(KM.FORCE, force_vector)  # penalty-free formulation
        model_part.AddCondition(condition)

        return condition

    def test_LoadSolidIGAConditionP3(self):
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
            -0.293744,0.0905711,-0.440616,0.135857,-0.220308,0.0679283,-0.036718,0.0113214,-0.0154602,0.0047669,-0.0231903,0.00715035,-0.0115952,0.00357517,-0.00193253,0.000595862
        ]
        tolerance = 1e-4

        self.assertEqual(rhs.Size(), len(expected_RHS))
        for i in range(rhs.Size()):
            self.assertAlmostEqual(rhs[i], expected_RHS[i], delta=tolerance)

if __name__ == '__main__':
    KratosUnittest.main()