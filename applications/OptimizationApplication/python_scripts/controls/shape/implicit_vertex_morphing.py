# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.OptimizationApplication.controls.shape.shape_control import ShapeControl

class ImplicitVertexMorphing(ShapeControl):

    def __init__(self, name, model, settings):
        super().__init__(name,model,settings)
        self.technique_settings = self.settings["technique_settings"]

        self.default_technique_settings = KM.Parameters("""{
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

        self.technique_settings.RecursivelyValidateAndAssignDefaults(self.default_technique_settings)

        ## Construct the linear solver
        import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
        self.linear_solver = python_linear_solver_factory.ConstructSolver(self.technique_settings["linear_solver_settings"])
        self.implicit_vertex_morphing = KOA.ImplicitVertexMorphing(self.name,self.model,self.linear_solver,self.settings)

    def Initialize(self):

        super().Initialize()

        self.implicit_vertex_morphing.Initialize()

        hkhk

        # self.ex_vm_mapper = {}
        # for model_part_name in self.controlling_objects:
        #     if not self.model.HasModelPart(model_part_name):
        #         raise RuntimeError("ImplicitVertexMorphing: Model part {} from control {} does not exist in the input model parts".format(model_part_name,self.name))
        #     ex_mapper = mapper_factory.CreateMapper(self.model.GetModelPart(model_part_name), self.model.GetModelPart(model_part_name), self.technique_settings)
        #     ex_mapper.Initialize()
        #     self.ex_vm_mapper[model_part_name] = ex_mapper
    

    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):
        pass
        # for mapper in self.ex_vm_mapper.values():
        #     mapper.InverseMap(derivative_variable_name,mapped_derivative_variable_name)

    def Compute(self):
        pass
        # for mapper in self.ex_vm_mapper.values():
        #     mapper.Map(KOA.D_CX,KOA.D_X)   

    def Update(self):
        pass
        # for model_part_name in self.controlling_objects:
        #     model_part = self.model.GetModelPart(model_part_name)
        #     for node in model_part.Nodes:
        #         shape_update = node.GetSolutionStepValue(KOA.D_X)
        #         node.X0 += shape_update[0]
        #         node.Y0 += shape_update[1]
        #         node.Z0 += shape_update[2]
        #         node.X += shape_update[0]
        #         node.Y += shape_update[1]
        #         node.Z += shape_update[2]   

        # for mapper in self.ex_vm_mapper.values():
        #     mapper.Update()               
            


