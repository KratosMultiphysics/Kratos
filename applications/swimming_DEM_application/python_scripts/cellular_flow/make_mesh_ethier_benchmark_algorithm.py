from KratosMultiphysics import *
import ethier_benchmark_algorithm
BaseAlgorithm = ethier_benchmark_algorithm.Algorithm


class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)

    def SetBetaParameters(self):
        BaseAlgorithm.SetBetaParameters(self)
        self.pp.CFD_DEM.AddEmptyValue("pressure_grad_recovery_type")
        self.pp.CFD_DEM.AddEmptyValue("size_parameter").SetInt(1)

    def SetCustomBetaParameters(self, custom_parameters): # TO DO: remove and make all calls to .size_parameter calls to Parameters object
        BaseAlgorithm.SetCustomBetaParameters(self, custom_parameters)
        self.pp.CFD_DEM.size_parameter = self.pp.CFD_DEM["size_parameter"].GetInt()

    def ReadFluidModelParts(self):
        from meshing import meshing_utilities
        self.mesh_generator = meshing_utilities.ParallelepipedRegularMesher(
                                                model_part_to_be_filled = self.fluid_solution.fluid_model_part,
                                                lower_corner_coordinates = [0.0, 0.0, 0.0],
                                                higher_corner_coordinates = [0.1, 0.1, 0.1],
                                                number_of_divisions_per_dimension = 10)

        self.mesh_generator.FillModelPartWithNewMesh()
