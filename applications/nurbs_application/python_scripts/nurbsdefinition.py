from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# NurbsDefinition is a class which contains all necessary information for reading in a NURBS-surface expect of
# the Control Points, which itself a read in by a seperate mdpa-file.
from KratosMultiphysics import *
class NurbsDefinition:
    def __init__(self,connectivities_start,connectivities_end, polynomial_degree_xi, polynomial_degree_eta, knots_xi, knots_eta, number_of_cp_xi, number_of_cp_eta, weights):
#       if connectivities will be connected it should be added in the constructor aswell
        self.connectivities_start = connectivities_start
        self.connectivities_end = connectivities_end
        self.polynomial_degree_xi = polynomial_degree_xi
        self.polynomial_degree_eta = polynomial_degree_eta
        self.knots_xi = Vector(len(knots_xi))
        self.knots_eta = Vector(len(knots_eta))
        self.number_of_cp_xi = number_of_cp_xi
        self.number_of_cp_eta = number_of_cp_eta
        self.weights = Vector(len(weights))

        for i in range(len(weights)):
            self.weights[i] = weights[i]

        for i in range(len(knots_xi)):
            self.knots_xi[i] = knots_xi[i]

        for i in range(len(knots_eta)):
            self.knots_eta[i] = knots_eta[i]




        pass
