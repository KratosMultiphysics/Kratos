from KratosMultiphysics.ExaquteSandboxApplication.WindGenerator.onthefly_kratos_inlet_process import ImposeWindInletProcess, Parameters, GenerateWind, get_extent
import KratosMultiphysics.MappingApplication
import KratosMultiphysics.mpi as KratosMPI
import time
def Factory(settings, Model):
    return ImposeMPIWindInletProcess(Model, Parameters(settings['Parameters']))

class ImposeMPIWindInletProcess(ImposeWindInletProcess):

    def __init__(self, Model, settings):
        for name, value in settings.items():
            setattr(self, name, value)
        #mpi modelpart
        self.mpi_model_part = Model[self.inlet_model_part_name]
        self.data_comm = self.mpi_model_part.GetCommunicator().GetDataCommunicator()
        # serial inlet
        if not Model.HasModelPart("wind_inlet_serial"):
            self.model_part = Model.CreateModelPart("wind_inlet_serial")
        else:
            self.model_part = Model.GetModelPart("wind_inlet_serial")
        # if self.data_comm.Rank() == 0:
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        self.wind_inlet_model_part_name = "problem_settings/wind_inlet_sub_model_part"
        self.x0 = 0.0

        if self.data_comm.Rank() == 0:
            KratosMultiphysics.ModelPartIO(self.wind_inlet_model_part_name).ReadModelPart(self.model_part)
        ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.model_part, self.data_comm)
        ParallelFillCommunicator.Execute()
        if self.data_comm.Rank() == 0:
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


    def UpdateInletPosition(self):
        dt = self.mpi_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.inlet_position += dt * self.mean_profile.bulk_wind_speed

        #If the end of the time block is reached, generate a new wind block and reset map
        if self.inlet_position >= self.x0 + self.lx:
            self.x0 += self.lx
            self.blend_region, self.mappers = self.Create3DMappers(self.wind, self.blend_region)

    def ApplyRamp(self):
        time = self.mpi_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if time < self.ramp_time:
            scal = time / self.ramp_time
            for node in self.inlet_nodes:
                for var, _ in self.mappers:
                    vel = node.GetSolutionStepValue(var)
                    node.SetSolutionStepValue(var, scal * vel)

    def ExecuteInitialize(self):
        for node in self.mpi_model_part.Nodes:
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)

    def ExecuteInitializeSolutionStep(self):
        if self.data_comm.Rank() == 0:
            super().ExecuteInitializeSolutionStep()

        KratosMultiphysics.Logger.PrintInfo("PERFORMING MAPPING")
        # map from current model part of interest to reference model part
        ini_time = time.time()
        mapping_parameters = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
            }""")
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMPIMapper(self.model_part, self.mpi_model_part,mapping_parameters)
        mapper.Map(KratosMultiphysics.VELOCITY, \
            KratosMultiphysics.VELOCITY)
        mapping_time = time.time() - ini_time
        KratosMultiphysics.Logger.PrintInfo("TIME SPENT MAPPING:", self.data_comm.MaxAll(mapping_time))
