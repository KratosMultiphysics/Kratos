import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetCosineDistances
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSubDistances

class TestSensorUtils(UnitTest.TestCase):
    def test_GetSubDistances(self):
        distances_matrix = [
            [0,  2,  3,  4,  5,  6],
            [2,  0,  7,  8,  9, 10],
            [3,  7,  0, 11, 12, 13],
            [4,  8, 11,  0, 14, 15],
            [5,  9, 12, 14,  0, 16],
            [6, 10, 13, 15, 16,  0]
        ]


        def get_distances_vector(matrix_of_distances):
            distances_vector: 'list[int]' = []
            for i, row in enumerate(matrix_of_distances):
                for col in row[i+1:]:
                    distances_vector.append(col)
            return distances_vector

        sub_distances = GetSubDistances(get_distances_vector(distances_matrix), [1,3,4,5])
        sub_matrix = []
        for i in [1,3,4,5]:
            row = []
            for j in [1,3,4,5]:
                row.append(distances_matrix[i][j])
            sub_matrix.append(row)

        self.assertEqual(sub_distances, get_distances_vector(sub_matrix))

if __name__ == '__main__':
    UnitTest.main()
