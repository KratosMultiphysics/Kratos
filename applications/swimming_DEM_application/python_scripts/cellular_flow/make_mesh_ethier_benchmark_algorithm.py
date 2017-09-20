from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import ethier_benchmark_algorithm
import swimming_DEM_procedures as SDP
import math
import numpy as np
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
        os.chdir(self.main_path)
        if ethier_benchmark_algorithm.num_type(self.pp.CFD_DEM.size_parameter) == 'int':
            model_part_io_fluid = ModelPartIO(self.pp.problem_name.replace('ethier', 'ethier_ndiv_' + str(self.pp.CFD_DEM.size_parameter) + 'GiD'))
        elif ethier_benchmark_algorithm.num_type(self.pp.CFD_DEM.size_parameter) == 'float':
            model_part_io_fluid = ModelPartIO(self.pp.problem_name.replace('ethier', 'ethier_h_' + str(self.pp.CFD_DEM.size_parameter)))
        # model_part_io_fluid.ReadModelPart(self.fluid_algorithm.fluid_model_part)
        print(self.fluid_algorithm.fluid_model_part)
        parameters = Parameters("{}")
        parameters.AddEmptyValue("element_name").SetString('VMS3D')
        parameters.AddEmptyValue("condition_name").SetString('WallCondition3D')
        parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        parameters.AddEmptyValue("number_of_divisions").SetInt(self.pp.CFD_DEM.size_parameter)
        # sample_geometry = Tetrahedra3D4(Node(1, 0.,0.,0.),
        #                                 Node(2, 1.,0.,0.),
        #                                 Node(3, 1.,1.,0.),
        #                                 Node(4, 0.,1.,0.))
        print(parameters)
        sample_geometry = Hexahedra3D8(Node(1, 0.1, 0.1, 0.0),
                                       Node(2, 0.0, 0.1, 0.0),
                                       Node(3, 0.0, 0.0, 0.0),
                                       Node(4, 0.1, 0.0, 0.0),
                                       Node(5, 0.1, 0.1, 0.1),
                                       Node(6, 0.0, 0.1, 0.1),
                                       Node(7, 0.0, 0.0, 0.1),
                                       Node(8, 0.1, 0.0, 0.1))
        # sample_geometry = Tetrahedra3D4(Point3D(0.,0.,0.),
        #                                 Point3D(1.,0.,0.),
        #                                 Point3D(1.,1.,0.),
        #                                 Point3D(0.,1.,0.))
        mesh_generator_process = StructuredMeshGeneratorProcess(sample_geometry, self.fluid_algorithm.fluid_model_part, parameters)
        print(parameters)
        print(mesh_generator_process)
        self.fluid_algorithm.fluid_model_part.CreateSubModelPart("Skin")
        # __del__(sample_geometry)
        mesh_generator_process.Execute()
        area_calculator = CalculateNodalAreaProcess(self.fluid_algorithm.fluid_model_part, 3)
        area_calculator.Execute()
        for node in self.fluid_algorithm.fluid_model_part.Nodes:
            print(node.GetSolutionStepValue(NODAL_AREA))
