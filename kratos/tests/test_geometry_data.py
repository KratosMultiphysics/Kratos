import KratosMultiphysics as kratos
from KratosMultiphysics import KratosUnittest

class TestGeometryData(KratosUnittest.TestCase):

    def setUp(self):
        self.model = kratos.Model()
        self.model_part = self.model.CreateModelPart("model_part")

        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        self.model_part.CreateNewNode(4, 0.5, 0.0, 0.0)
        self.model_part.CreateNewNode(5, 0.5, 0.5, 0.0)
        self.model_part.CreateNewNode(6, 0.0, 0.5, 0.0)

        self.model_part.CreateNewNode(7, 0.0, 0.0, 1.0)

        prop = self.model_part.GetProperties()[1]

        self.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
        self.model_part.CreateNewElement("Element2D6N", 2, [1, 2, 3, 4, 5, 6], prop)
        self.model_part.CreateNewElement("Element3D4N", 3, [1, 2, 3, 7], prop)

        self.model_part.CreateNewCondition("LineCondition3D2N", 1, [1, 7], prop)


    def test_geometry_data(self):

        triangle_2d_3 = self.model_part.Elements[1].GetGeometry()
        self.assertEqual(
            triangle_2d_3.GetDefaultIntegrationMethod(),
            kratos.GeometryData.IntegrationMethod.GI_GAUSS_1)
        self.assertEqual(
            triangle_2d_3.GetGeometryFamily(),
            kratos.GeometryData.KratosGeometryFamily.Kratos_Triangle)
        self.assertEqual(
            triangle_2d_3.GetGeometryType(),
            kratos.GeometryData.KratosGeometryType.Kratos_Triangle2D3)

        triangle_2d_6 = self.model_part.Elements[2].GetGeometry()
        self.assertEqual(
            triangle_2d_6.GetDefaultIntegrationMethod(),
            kratos.GeometryData.IntegrationMethod.GI_GAUSS_2)
        self.assertEqual(
            triangle_2d_6.GetGeometryFamily(),
            kratos.GeometryData.KratosGeometryFamily.Kratos_Triangle)
        self.assertEqual(
            triangle_2d_6.GetGeometryType(),
            kratos.GeometryData.KratosGeometryType.Kratos_Triangle2D6)

        tetrahedron_3d_4 = self.model_part.Elements[3].GetGeometry()
        self.assertEqual(
            tetrahedron_3d_4.GetDefaultIntegrationMethod(),
            kratos.GeometryData.IntegrationMethod.GI_GAUSS_1)
        self.assertEqual(
            tetrahedron_3d_4.GetGeometryFamily(),
            kratos.GeometryData.KratosGeometryFamily.Kratos_Tetrahedra)
        self.assertEqual(
            tetrahedron_3d_4.GetGeometryType(),
            kratos.GeometryData.KratosGeometryType.Kratos_Tetrahedra3D4)

        line_3d_2 = self.model_part.Conditions[1].GetGeometry()
        self.assertEqual(
            line_3d_2.GetDefaultIntegrationMethod(),
            kratos.GeometryData.IntegrationMethod.GI_GAUSS_1)
        self.assertEqual(
            line_3d_2.GetGeometryFamily(),
            kratos.GeometryData.KratosGeometryFamily.Kratos_Linear)
        self.assertEqual(
            line_3d_2.GetGeometryType(),
            kratos.GeometryData.KratosGeometryType.Kratos_Line3D2)

if __name__ == "__main__":
    KratosUnittest.main()