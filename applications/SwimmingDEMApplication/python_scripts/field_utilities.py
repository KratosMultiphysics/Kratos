from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *

class PorosityField:

    def __init__(self, real_field_type, fixed_mesh_option):
        self.real_field_type = real_field_type
        self.fixed_mesh_option = fixed_mesh_option

        if (self.real_field_type == 0):
            self.force_formula = TimeDependantForceField(2)
            self.porosity_default_value = 1.0
            self.force_default_vector = Array3()
            self.force_default_vector[0] = 0.0
            self.force_default_vector[1] = 0.0
            self.force_default_vector[2] = 0.0
            self.porosity_formula = self.force_formula.GetPorosityField()
            self.b_box_rule = BoundingBoxRule()
            self.b_box_rule.SetTimeBoundingInterval(0,3)
            self.b_box_rule.SetYBoundingInterval(-1,1)
            self.domain = SpaceTimeSet()
            self.domain.AddAndRule(self.b_box_rule)
            self.porosity_field_utility = FieldUtility(self.domain)
            self.force_field_utility = FieldUtility(self.domain)

    def EvaluatePorosityAtPoint(self, time, coor):
        return self.porosity_field_utility.EvaluateFieldAtPoint(time, coor, self.porosity_formula)

    def ImposePorosityField(self, model_part):
        self.porosity_field_utility.MarkNodesInside(model_part, model_part.ProcessInfo)
        self.porosity_field_utility.ImposeFieldOnNodes(FLUID_FRACTION, self.porosity_default_value,  self.porosity_formula, model_part, model_part.ProcessInfo, False)

    def ImposeForceField(self, model_part):
        self.force_field_utility.MarkNodesInside(model_part, model_part.ProcessInfo)
        self.force_field_utility.ImposeFieldOnNodes(BODY_FORCE, self.force_default_vector, self.force_formula, model_part, model_part.ProcessInfo, False)
