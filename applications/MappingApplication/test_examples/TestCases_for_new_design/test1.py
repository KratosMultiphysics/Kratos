from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

model_part_origin = KratosMultiphysics.ModelPart("NumberOne")
model_part_destination = KratosMultiphysics.ModelPart("NumberTwo")

model_part_origin.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
model_part_destination.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

model_part_origin.CreateNewNode(1, 1.1, 2.2, 3.3)
model_part_destination.CreateNewNode(1, 1.4, 2.5, 3.6)

mapper_settings = KratosMultiphysics.Parameters("""
        {
            "mapper_type": "NearestElement",
            "interface_submodel_part_origin": "interface_chimera_background",
            "interface_submodel_part_destination": "Inlet2D_inlet"
        }
""")

mapper = KratosMapping.NearestElementMapper(model_part_origin, model_part_destination, mapper_settings)

mapper.Map(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE, KratosMapping.Mapper.ADD_VALUES)

mapper.UpdateInterface()

KratosMapping.MapperFactoryNew.CreateMapper(model_part_origin, model_part_destination, mapper_settings)