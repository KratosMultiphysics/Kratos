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

        self.implicit_settings = self.settings["implicit_vm_settings"]
        self.implicit_settings.RecursivelyValidateAndAssignDefaults(self.GetDefaultImplicitVMSettings())
        self.implicit_settings.AddString("design_surface_sub_model_part_name",str(self.destination_model_part.Name))
        self.implicit_settings.AddBool("export_linear_system",bool(self.settings["export_linear_system"]))


        ## Construct the linear solver
        import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
        self.linear_solver = python_linear_solver_factory.ConstructSolver(self.implicit_settings["linear_solver_settings"])

        self.im_vm_mapper = KSO.MapperImplicitVertexMorphing(self.origin_model_part, self.linear_solver, self.implicit_settings)


    @classmethod
    def GetDefaultImplicitVMSettings(cls):
        return KM.Parameters("""{
            "element_type" : "helmholtz_vec_element",
            "surface_element_type" : "helmholtz_surf_element",
            "only_design_surface_parameterization" : true,
            "formulate_on_the_undeformed_configuration" : true,
            "automatic_filter_size" : true,
            "adaptive_filter_size" : false,
            "surface_filter_radius" : 0.000000000001,            
            "bulk_filter_radius" : 0.000000000001,
            "poisson_ratio" : 0.3,            
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

