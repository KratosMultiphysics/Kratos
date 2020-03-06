# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================
import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
import KratosMultiphysics.scipy_conversion_tools as scipy_conversion_tools
from .make_consistent_aat import MakeConsistent

class ConsistentVertexMorphingMapper():
    """
    The ConsistentVertexMorphingMapper extends the standard Vertex Morphing approach.
    """

    def __init__(self, origin_model_part, destination_model_part, settings):
        self.settings = settings
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part

        extracted_vm_settings = settings.Clone()
        extracted_vm_settings["consistent_aat"].SetBool(False)
        self.vm_mapper = mapper_factory.CreateMapper(origin_model_part, destination_model_part, extracted_vm_settings)

        self.matrix = None

    def Initialize(self):
        self.vm_mapper.Initialize()
        self.matrix = self._GetScipyMatrix(self.vm_mapper)
        self.matrix = MakeConsistent(self.matrix)

    def Update(self):
        self.vm_mapper.Update()
        self.matrix = self._GetScipyMatrix(self.vm_mapper)
        self.matrix = MakeConsistent(self.matrix)

    def Map(self, origin_variable, destination_variable):
        for i in range(3):
            origin_input = self._AssembleVariableToVector(
                self.origin_model_part,
                origin_variable,
                i)

            destination_output = self.matrix @ origin_input

            self._AssignVectorToVariable(
                self.destination_model_part,
                destination_output,
                destination_variable,
                i)

    def InverseMap(self, destination_variable, origin_variable):
        for i in range(3):
            destination_input = self._AssembleVariableToVector(
                self.destination_model_part,
                destination_variable,
                i)

            origin_output = self.matrix.transpose() @ destination_input

            self._AssignVectorToVariable(
                self.origin_model_part,
                origin_output,
                origin_variable,
                i)

    @staticmethod
    def _AssembleVariableToVector(model_part, variable, dim):
        values = np.zeros(model_part.NumberOfNodes())
        for i, node in enumerate(model_part.Nodes):
            _v = node.GetSolutionStepValue(variable)
            values[i] = _v[dim]
        return values

    @staticmethod
    def _AssignVectorToVariable(model_part, values, variable, dim):
        for i, node in enumerate(model_part.Nodes):
            v = node.GetSolutionStepValue(variable)
            v[dim] = values[i]
            node.SetSolutionStepValue(variable, v)

    @staticmethod
    def _GetScipyMatrix(mapper):
        cpp_matrix = mapper.GetMatrix()
        return scipy_conversion_tools.to_csr(cpp_matrix)


