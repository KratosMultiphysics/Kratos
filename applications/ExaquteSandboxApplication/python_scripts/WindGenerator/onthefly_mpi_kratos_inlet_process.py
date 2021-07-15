from KratosMultiphysics.ExaquteSandboxApplication.WindGenerator.onthefly_kratos_inlet_process import ImposeWindInletProcess, Parameters, GenerateWind, get_extent
import KratosMultiphysics.MappingApplication

def Factory(settings, Model):
    return ImposeMPIWindInletProcess(Model, Parameters(settings['Parameters']))

class ImposeMPIWindInletProcess(ImposeWindInletProcess):

    def __init__(self, Model, settings):
        for name, value in settings.items():
            setattr(self, name, value)
        #mpi modelpart
        self.mpi_model_part = Model[self.inlet_model_part_name]
        # serial inlet
        self.model_part = Model.CreateModelPart("wind_inlet_serial")
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.wind_inlet_model_part_name = "wind_inlet_sub_model_part"
        KratosMultiphysics.ModelPartIO(self.wind_inlet_model_part_name).ReadModelPart(self.model_part)
        self.x0 = 0.0

        if len(self.inlet_nodes) > 0:
            # define dimensions of the inlet domain
            y_extent = get_extent(self.inlet_nodes, key=lambda node: node.Y)
            z_extent = get_extent(self.inlet_nodes, key=lambda node: node.Z)
            setattr(self, 'ly', y_extent.upper - y_extent.lower)
            setattr(self, 'lz', z_extent.upper - z_extent.lower)
            self.mean_profile = self.CreateMeanProfile()
            # x_extent = Extent(0, self.time_interval_length * self.mean_profile.bulk_wind_speed)
            # setattr(self, 'lx', x_extent.upper - x_extent.lower)
            self.inlet_position = 0.0

            # define wind field and map
            grid_dimensions = [self.lx, self.ly, self.lz]
            self.wind = GenerateWind(self.friction_velocity, self.reference_height, grid_dimensions, self.wind_grid_levels, self.seed)
            self.block_num=0
            self.blend_region, self.mappers = self.Create3DMappers(self.wind)

    def ExecuteInitialize(self):
        for node in self.mpi_model_part.Nodes:
            for var, _ in self.mappers:
                node.Fix(var)

    def ExecuteInitializeSolutionStep(self):
        super().ExecuteInitializeSolutionStep()
        print("PERFORMING MAPPING")
        # map from current model part of interest to reference model part
        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
            }""")
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMPIMapper(self.model_part, self.mpi_model_part,mapping_parameters)
        mapper.Map(KratosMultiphysics.VELOCITY, \
            KratosMultiphysics.VELOCITY)
