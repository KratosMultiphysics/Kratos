from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils

dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
if dependencies_are_available:
    import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

class TestSpecificationsUtilities(KratosUnittest.TestCase):
    def test_specifications_utilities_elements(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        node1 = model_part.CreateNewNode(1, 0.0,0.0,0.0)
        node2 = model_part.CreateNewNode(2, 1.0,0.0,0.0)
        node3 = model_part.CreateNewNode(3, 1.0,1.0,0.0)
        prop1 = KratosMultiphysics.Properties(1)
        model_part.AddProperties(prop1)
        elem1 = model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])
        elem1.Initialize(model_part.ProcessInfo)

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        KratosMultiphysics.SpecificationsUtilities.AddMissingVariables(model_part)
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)

        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_X), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Y), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z), False)
        KratosMultiphysics.SpecificationsUtilities.AddMissingDofs(model_part)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_X), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Y), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z), False)

        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineTimeIntegration(model_part), [])
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineFramework(model_part), "lagrangian")
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineSymmetricLHS(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DeterminePositiveDefiniteLHS(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineIfCompatibleGeometries(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineIfRequiresTimeIntegration(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckCompatibleConstitutiveLaws(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckGeometricalPolynomialDegree(model_part), -1)
        docu = KratosMultiphysics.Parameters("""{"Element2D3N"   : "This is a pure geometric element, no computation"}""")
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.GetDocumention(model_part).IsEquivalentTo(docu), True)

    def test_specifications_utilities_conditions(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        node1 = model_part.CreateNewNode(1, 0.0,0.0,0.0)
        node2 = model_part.CreateNewNode(2, 1.0,0.0,0.0)
        node3 = model_part.CreateNewNode(3, 1.0,1.0,0.0)
        prop1 = KratosMultiphysics.Properties(1)
        model_part.AddProperties(prop1)
        cond1 = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,3], model_part.GetProperties()[1])
        cond1.Initialize(model_part.ProcessInfo)

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        KratosMultiphysics.SpecificationsUtilities.AddMissingVariables(model_part)
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)

        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_X), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Y), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z), False)
        KratosMultiphysics.SpecificationsUtilities.AddMissingDofs(model_part)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_X), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Y), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z), False)

        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineTimeIntegration(model_part), [])
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineFramework(model_part), "lagrangian")
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineSymmetricLHS(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DeterminePositiveDefiniteLHS(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineIfCompatibleGeometries(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineIfRequiresTimeIntegration(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckCompatibleConstitutiveLaws(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckGeometricalPolynomialDegree(model_part), -1)
        docu = KratosMultiphysics.Parameters("""{"SurfaceCondition3D3N"   : "This is a pure geometric condition, no computation"}""")
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.GetDocumention(model_part).IsEquivalentTo(docu), True)

    @KratosUnittest.skipUnless(dependencies_are_available,"StructuralMechanicsApplication is not available")
    def test_specifications_utilities_elements_dependencies(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        node1 = model_part.CreateNewNode(1, 0.0,0.0,0.0)
        node2 = model_part.CreateNewNode(2, 1.0,0.0,0.0)
        node3 = model_part.CreateNewNode(3, 1.0,1.0,0.0)
        prop1 = KratosMultiphysics.Properties(1)
        model_part.AddProperties(prop1)
        prop1.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, StructuralMechanicsApplication.LinearElasticPlaneStrain2DLaw())
        elem1 = model_part.CreateNewElement("SmallDisplacementElement2D3N", 1, [1,2,3], model_part.GetProperties()[1])
        elem1.Initialize(model_part.ProcessInfo)

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        KratosMultiphysics.SpecificationsUtilities.AddMissingVariables(model_part)
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 3)

        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_X), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Y), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z), False)
        KratosMultiphysics.SpecificationsUtilities.AddMissingDofs(model_part)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_X), True)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Y), True)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z), True)

        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineTimeIntegration(model_part).sort(), ['explicit', 'static', 'implicit'].sort())
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineFramework(model_part), "lagrangian")
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineSymmetricLHS(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DeterminePositiveDefiniteLHS(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineIfCompatibleGeometries(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineIfRequiresTimeIntegration(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckCompatibleConstitutiveLaws(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckGeometricalPolynomialDegree(model_part), -1)
        docu = KratosMultiphysics.Parameters("""{"SmallDisplacementElement2D3N"   : "This is a pure displacement element"}""")
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.GetDocumention(model_part).IsEquivalentTo(docu), True)

        # Changing the law to get a False check
        prop1.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw())
        elem1.Initialize(model_part.ProcessInfo)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckCompatibleConstitutiveLaws(model_part), False)

    @KratosUnittest.skipUnless(dependencies_are_available,"StructuralMechanicsApplication is not available")
    def test_specifications_utilities_conditions_dependencies(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        node1 = model_part.CreateNewNode(1, 0.0,0.0,0.0)
        node2 = model_part.CreateNewNode(2, 1.0,0.0,0.0)
        node3 = model_part.CreateNewNode(3, 1.0,1.0,0.0)
        prop1 = KratosMultiphysics.Properties(1)
        model_part.AddProperties(prop1)
        cond1 = model_part.CreateNewCondition("SurfaceLoadCondition3D3N", 1, [1,2,3], model_part.GetProperties()[1])
        cond1.Initialize(model_part.ProcessInfo)

        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 0)
        KratosMultiphysics.SpecificationsUtilities.AddMissingVariables(model_part)
        self.assertEqual(model_part.GetNodalSolutionStepDataSize(), 3)

        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_X), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Y), False)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z), False)
        KratosMultiphysics.SpecificationsUtilities.AddMissingDofs(model_part)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_X), True)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Y), True)
        self.assertEqual(node1.HasDofFor(KratosMultiphysics.DISPLACEMENT_Z), True)

        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineTimeIntegration(model_part).sort(), ['explicit', 'static', 'implicit'].sort())
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineFramework(model_part), "lagrangian")
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineSymmetricLHS(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DeterminePositiveDefiniteLHS(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineIfCompatibleGeometries(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.DetermineIfRequiresTimeIntegration(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckCompatibleConstitutiveLaws(model_part), True)
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.CheckGeometricalPolynomialDegree(model_part), -1)
        docu = KratosMultiphysics.Parameters("""{"SurfaceLoadCondition3D3N"   : "This is a pure displacement condition"}""")
        self.assertEqual(KratosMultiphysics.SpecificationsUtilities.GetDocumention(model_part).IsEquivalentTo(docu), True)

if __name__ == '__main__':
    KratosUnittest.main()
