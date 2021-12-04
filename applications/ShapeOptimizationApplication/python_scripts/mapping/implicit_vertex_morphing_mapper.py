# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl
#
# ==============================================================================

import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

class ImplicitVertexMorphingMapper():
    """

    """

    def __init__(self, origin_model_part, destination_model_part, settings):

        self.settings = settings
        self.origin_model_part = origin_model_part
        self.destination_model_part = destination_model_part

        implicit_settings = self.settings["implicit_vm_settings"]
        implicit_settings.RecursivelyValidateAndAssignDefaults(self.GetDefaultImplicitVMSettings())

        ## Construct the linear solver
        import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
        self.linear_solver = python_linear_solver_factory.ConstructSolver(implicit_settings["linear_solver_settings"])

        self.im_vm_mapper = KSO.MapperImplicitVertexMorphing(self.origin_model_part, self.linear_solver, self.settings)


    @classmethod
    def GetDefaultImplicitVMSettings(cls):
        return KM.Parameters("""{
            "linear_solver_settings" : {
                "solver_type" : "amgcl",
                "smoother_type":"ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation",
                "max_iteration": 200,
                "provide_coordinates": false,
                "gmres_krylov_space_dimension": 100,
                "verbosity" : 0,
                "tolerance": 1e-7,
                "scaling": false,
                "block_size": 1,
                "use_block_matrices_if_possible" : true,
                "coarse_enough" : 5000
            }
        }""")

    def Initialize(self):
        self.im_vm_mapper.Initialize()

    def Update(self):
        self.im_vm_mapper.Update()

    def Map(self, origin_variable, destination_variable):
        self.im_vm_mapper.Map(origin_variable, destination_variable)

    def InverseMap(self, destination_variable, origin_variable):
        self.im_vm_mapper.InverseMap(destination_variable, origin_variable)

