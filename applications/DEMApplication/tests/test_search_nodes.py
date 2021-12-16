import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication as DEM

class TestSearchNodes(KratosUnittest.TestCase):

    def setUp(self):
        self.model = Kratos.Model()
        self.base_model_part = self.model.CreateModelPart('Base')
        self.target_model_part = self.model.CreateModelPart('Target')
        # creating a search tool with a periodic bounding box
        dx = self.unit_length = 0.1
        self.epsilon = 0.0001
        self.domain_box = [0, 0, 0, 10*dx, 10*dx, 10*dx]
        dimension = int(len(self.domain_box) / 2)
        self.domain_box_periods = [self.domain_box[i + dimension] - self.domain_box[i] for i in range(dimension)]
        self.search_strategy = DEM.OMP_DEMSearch(*self.domain_box)
        self.CreateNodes()

    def test_SearchNodesInTargetModelPart(self):
        dx = self.unit_length
        epsilon = self.epsilon

        neighbours = []
        distances = []
        # setting search radius just about large enough to catch some of the neighbours
        search_radius = 2*dx + epsilon
        radii = [search_radius for node in self.base_model_part.Nodes]
        self.search_strategy.SearchNodesInRadiusExclusive(self.base_model_part.Nodes,
                                                          self.target_model_part.Nodes,
                                                          radii,
                                                          DEM.VectorResultNodesContainer(),
                                                          DEM.VectorDistances(),
                                                          neighbours,
                                                          distances)

        neighbours = [[n for n in l] for l in neighbours]
        distances = [[d for d in l] for l in distances]
        print('neighbours', neighbours)
        print('distances', distances)

        obtained_distances_dictionary = dict()
        for i, l in enumerate(neighbours):
            for j, n in enumerate(l):
                obtained_distances_dictionary[frozenset([self.base_ids[i], self.target_ids[j]])] = distances[i][j]

        for k in obtained_distances_dictionary.keys():
            print('obtained pair', obtained_distances_dictionary[k], self.expected_distances[k])

        self.AssertCorrectnessOfNeighbours(neighbours, search_radius)
        self.AssertCorrectnessOfDistances(distances, neighbours)

        neighbours = []
        search_radius = 2*dx - epsilon

        # setting search radius a tad too short to catch some of the neighbours
        radii = [search_radius for node in self.base_model_part.Nodes]
        self.search_strategy.SearchNodesInRadiusExclusive(self.base_model_part.Nodes,
                                                          self.target_model_part.Nodes,
                                                          radii,
                                                          DEM.VectorResultNodesContainer(),
                                                          DEM.VectorDistances(),
                                                          neighbours,
                                                          distances)

        self.AssertCorrectnessOfNeighbours(neighbours, search_radius)
        self.AssertCorrectnessOfDistances(distances, neighbours)

    def GetNeighborPairs(self, neighbours):
        pairs = []
        for i, l in enumerate(neighbours):
            for j, n in enumerate(l):
                pairs.append([self.base_ids[i], n])
        return pairs

    def AssertCorrectnessOfNeighbours(self, neighbours, search_radius):
        neighbours_sets = [frozenset(l) for l in neighbours]
        self.assertEqual(neighbours_sets, self.GetExpectedListsOfNeighbors(search_radius))

    def AssertCorrectnessOfDistances(self, distances, neighbours):
        obtained_distances_dictionary = dict()
        for i, l in enumerate(neighbours):
            for j, n in enumerate(l):
                obtained_distances_dictionary[frozenset([self.base_ids[i], self.target_ids[j]])] = distances[i][j]

        for k in obtained_distances_dictionary.keys():
            self.assertAlmostEqual(obtained_distances_dictionary[k], self.expected_distances[k])


    def CreateNodes(self):
        dx = self.unit_length
        # Base nodes, for which neighbours must be searched for
        # ------------------------------------------
        B1 = [dx, dx, dx]
        B2 = [0.9*dx, 0.9*dx, 0.9*dx]
        base_points = [B1, B2]
        # Target nodes, which are search candidates
        # ------------------------------------------
        # the first one is simply displaced two length units in the x-direction with respect to
        # the first base node
        T1 = [3*dx, dx, dx]
        # the second candidate is displaced along the diagonal two length units towards the origin,
        # so that it comes in through the opposite end of the diagonal
        T2 = [(11 - 2.0/3**0.5) * dx for x in T1]
        # the third one is displaced 3 units in the y direction from the first base node
        T3 = [4*dx, dx, dx]
        # the fourth candidate is also displaced 3 units in the y direction, but in the opposite
        # direction (through the boundary)
        T4 = [1 - 2*dx, dx, dx]

        disp = self.GetPeriodicDisplacement(B1, T2)

        target_points = [T1, T2, T3, T4]

        self.base_ids = [i + 1 for i in range(len(base_points))]
        base_points = dict(zip(self.base_ids, base_points))

        for id in base_points.keys():
             self.CreateNode('Base', id, base_points[id])

        self.target_ids = [i + len(self.base_ids) + 1 for i in range(len(target_points))]
        target_points = dict(zip(self.target_ids, target_points))

        for id in target_points.keys():
             self.CreateNode('Target', id, target_points[id])

        self.expected_distances = dict()
        for base_id in self.base_ids:
            for target_id in self.target_ids:
                self.expected_distances[frozenset([base_id, target_id])] = sum(d**2 for d in self.GetPeriodicDisplacement(base_points[base_id], target_points[target_id]))**0.5


    def CreateNode(self, model_part_name, id, coordinates):
        model_part = self.model.GetModelPart(model_part_name)
        node = model_part.CreateNewNode(id, coordinates[0], coordinates[1], coordinates[2])

    def GetPeriodicDisplacement(self, X, Y):
        difference = [x - y for (x, y) in zip(X, Y)]
        sign = lambda x: -1 if x < 0 else 1
        for i, d in enumerate(difference):
            thorugh_boundary_difference = d - sign(d) * self.domain_box_periods[i]
            if abs(d) > abs(thorugh_boundary_difference):
                difference[i] = thorugh_boundary_difference
        return difference

    def GetExpectedListsOfNeighbors(self, search_radius):
        expected_neighbours = [frozenset() for id in self.base_ids]

        for i, base_id in enumerate(self.base_ids):
            neighs = set()
            for target_id in self.target_ids:
                expected_distance = self.expected_distances[frozenset([base_id, target_id])]
                if expected_distance < search_radius:
                    neighs.add(target_id)
                expected_neighbours[i] = frozenset(neighs)

        return expected_neighbours

if __name__ == '__main__':
    KratosUnittest.main()
