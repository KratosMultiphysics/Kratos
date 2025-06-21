import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication as DEM

# This test consists in applying OMP_DEMSearch to assign, to each node in
# a 'Base' model part, the corresponding neighbors from a 'Target' model part.

# The results calculated by Kratos are compared to the expected results. We
# assume the domain to be a periodic parallepiped, so that the closest neighbors
# are to be found across the boundaries in some cases.

# In particular, it is checked that the search is sharp, in the sense that if we
# vary the search radius by a small amount (epsilon), neighbors placed strategically
# at a distance equal to the base search radius are either found or not.

class TestSearchNodes(KratosUnittest.TestCase):

    def setUp(self):
        # input parameters:
        dx = self.unit_length = 0.1
        self.epsilon = 0.0001 # small distance by which we modify the sensitivity to search radius
        self.domain_box = [0, 0, 0, 10*dx, 10*dx, 10*dx]

        # creating a search tool with a periodic bounding box
        dimension = int(len(self.domain_box) / 2)
        self.domain_box_periods = [self.domain_box[i + dimension] - self.domain_box[i] for i in range(dimension)]
        self.search_strategy = DEM.OMP_DEMSearch(*self.domain_box)

        # creating data containers and initializing global variables
        self.next_node_id = 1
        self.all_points_db = dict()
        self.expected_distances_db = dict()
        self.base_ids = []
        self.target_ids = []

        # creating model parts
        self.model = Kratos.Model()
        self.base_model_part = self.model.CreateModelPart('Base')
        self.target_model_part = self.model.CreateModelPart('Target')

        # creating nodes
        self.CreateNodes()

    def test_SearchNodesInTargetModelPart(self):
        dx = self.unit_length
        epsilon = self.epsilon

        neighbors, distances = [], []
        # setting search radius just about large enough to catch some of the potential neighbors
        search_radius = 2*dx + epsilon
        radii = [search_radius for node in self.base_model_part.Nodes]
        self.search_strategy.SearchNodesInRadiusExclusive(self.base_model_part.Nodes,
                                                          self.target_model_part.Nodes,
                                                          radii,
                                                          DEM.VectorResultNodesContainer(),
                                                          DEM.VectorDistances(),
                                                          neighbors,
                                                          distances)

        self.AssertCorrectnessOfNeighbors(neighbors, search_radius)
        self.AssertCorrectnessOfDistances(distances, neighbors)

        neighbors, distances = [], []
        search_radius = 2*dx - epsilon
        # setting search radius a tad too short to catch some of the potential neighbors
        radii = [search_radius for node in self.base_model_part.Nodes]
        self.search_strategy.SearchNodesInRadiusExclusive(self.base_model_part.Nodes,
                                                          self.target_model_part.Nodes,
                                                          radii,
                                                          DEM.VectorResultNodesContainer(),
                                                          DEM.VectorDistances(),
                                                          neighbors,
                                                          distances)

        self.AssertCorrectnessOfNeighbors(neighbors, search_radius)
        self.AssertCorrectnessOfDistances(distances, neighbors)

    def AssertCorrectnessOfNeighbors(self, neighbors, search_radius):
        obtained_neighbors = [frozenset(l) for l in neighbors]
        expected_neighbors = self.GetExpectedListsOfNeighbors(search_radius)
        self.assertEqual(obtained_neighbors, expected_neighbors)

    def AssertCorrectnessOfDistances(self, distances, neighbor_lists):
        for i, i_neighbors in enumerate(neighbor_lists):
            for j, neigh_id in enumerate(i_neighbors):
                expected_distance = self.GetDistanceById(self.base_ids[i], neigh_id)
                obtained_distance = distances[i][j]
                self.assertAlmostEqual(expected_distance, obtained_distance)

    def CreateNodes(self):
        dx = self.unit_length
        # Base nodes, for which neighbors must be searched for
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
        target_points = [T1, T2, T3, T4]

        for point in base_points:
             self.CreateNode('Base', point)

        for point in target_points:
             self.CreateNode('Target', point)

    @staticmethod
    def Norm(V):
        return sum(v**2 for v in V)**0.5

    def CreateNode(self, model_part_name, coordinates):
        model_part = self.model.GetModelPart(model_part_name)
        id = self.next_node_id
        model_part.CreateNewNode(id, *coordinates)
        self.all_points_db[id] = coordinates

        if model_part_name == 'Base':
            self.base_ids.append(id)
        else:
            self.target_ids.append(id)

        self.next_node_id += 1

    def GetExpectedListsOfNeighbors(self, search_radius):
        expected_neighbors = [None for __ in self.base_ids]

        for i, base_id in enumerate(self.base_ids):
            neighs = set()
            for target_id in self.target_ids:
                if self.GetDistanceById(base_id, target_id) < search_radius:
                    neighs.add(target_id)
            expected_neighbors[i] = frozenset(neighs)

        return expected_neighbors

    def GetDistanceById(self, id1, id2):
        X = self.all_points_db[id1]
        Y = self.all_points_db[id2]
        distance = type(self).Norm(self.GetPeriodicDisplacement(X, Y))
        return distance

    def GetPeriodicDisplacement(self, X, Y):
        sign = lambda x: -1 if x < 0 else 1
        difference = [x - y for (x, y) in zip(X, Y)]
        for i, d in enumerate(difference):
            thorugh_boundary_difference = d - sign(d) * self.domain_box_periods[i]
            if abs(d) > abs(thorugh_boundary_difference):
                difference[i] = thorugh_boundary_difference
        return difference

if __name__ == '__main__':
    KratosUnittest.main()
