from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.KratosUnittest as KratosUnittest

mdpa_file_name_beam    = "mdpa_files/blade_quad"
mdpa_file_name_surface = "mdpa_files/blade_tri"

class TestBeamMapper(KratosUnittest.TestCase):
    def setUp(self):
        self.current_model = KM.Model()
        self.model_part_beam = self.current_model.CreateModelPart("beam")
        self.model_part_surface = self.current_model.CreateModelPart("surface")

        # list of variables involved in the Mapper-Tests
        self.model_part_beam.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part_beam.AddNodalSolutionStepVariable(KM.FORCE)

        self.model_part_surface.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        self.model_part_surface.AddNodalSolutionStepVariable(KM.VELOCITY)

        KM.ModelPartIO(mdpa_file_name_beam).ReadModelPart(self.model_part_beam)
        KM.ModelPartIO(mdpa_file_name_surface).ReadModelPart(self.model_part_surface)

    def test_beam_mapper(self):
        mapper_settings = KM.Parameters("""{
            "mapper_type": "beam_mapper",
            "echo_level" : 0
        }""")

        self.mapper = KratosMapping.MapperFactory.CreateMapper(self.model_part_beam, self.model_part_surface, mapper_settings)

        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        # self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)

if __name__ == '__main__':
    KratosUnittest.main()