#/////////////////////////////////////////////////
#// Main author: Chengshun Shang (CIMNE)
#// Email: cshang@cimne.upc.edu // chengshun.shang1996@gmail.com
#// Date: April 2024
#/////////////////////////////////////////////////

import numpy as np

'''
how to build a k-d tree?

example:

    points_list = []

    for node in self.spheres_model_part.Nodes:
        
        r = node.GetSolutionStepValue(RADIUS)
        x = node.X
        y = node.Y
        z = node.Z

        points_list.append(((x, y, z), r))

    self.kdtree = Kdtree()

    self.kdtree.BuildKdtree(points_list)

    target_point = (center_x, center_y, center_z)
    resulted_nodes = []

    resulted_nodes = self.kdtree.SearchKdtree(self.kdtree.Kdtree, target_point, radius, resulted_nodes)
'''

class KdtreeNode:
    def __init__(self, point, radius, left=None, right=None, axis=None):
        self.point = point
        self.radius = radius
        self.left = left
        self.right = right
        self.axis = axis

class Kdtree:
    def __init__(self):
        self.Kdtree = None

    def BuildKdtree(self, points, depth=0):
        if len(points) == 0:
            return None
        
        axis = depth % len(points[0])
        
        points.sort(key=lambda x: x[0][axis])
        
        median_idx = len(points) // 2
        median_point = points[median_idx][0]
        median_radius = points[median_idx][1]
        
        self.Kdtree = KdtreeNode(
                                point = median_point,
                                radius = median_radius,
                                left = self.BuildKdtree(points[:median_idx], depth + 1),
                                right = self.BuildKdtree(points[median_idx + 1:], depth + 1),
                                axis = axis
                            )

        return self.Kdtree

    def SearchKdtree(self, node, target_point, radius, result=[]):
        if node is None:
            return

        dist = np.linalg.norm(np.array(target_point) - np.array(node.point))
        
        if dist <= (radius + node.radius):
            result.append((node.point, node.radius))
        
        if target_point[node.axis] - (radius + node.radius) <= node.point[node.axis]:
            self.SearchKdtree(node.left, target_point, radius, result)
        if target_point[node.axis] + (radius + node.radius) >= node.point[node.axis]:
            self.SearchKdtree(node.right, target_point, radius, result)

        return result
