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
import math

import numpy as np

class ImplicitVertexMorphing(ShapeControl):
    def __init__(self, name, model, settings):
        self.technique_settings = settings["technique_settings"]
        self.default_technique_settings = KM.Parameters("""{
                    "fixed_model_parts"  : [],
                    "fixed_model_parts_X"  : [],
                    "fixed_model_parts_Y"  : [],
                    "fixed_model_parts_Z"  : [],
                    "only_design_surface_parameterization" : true,
                    "formulate_on_the_undeformed_configuration" : true,
                    "automatic_filter_size" : true,
                    "adaptive_filter_size" : false,
                    "surface_filter_radius" : 0.000000000001,
                    "surface_bulk_ratio" : 2,
                    "project_to_normal" : false,
                    "poisson_ratio" : 0.3,
                    "utilities":[],
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

        self.technique_settings.ValidateAndAssignDefaults(self.default_technique_settings)
        self.project_to_normal = self.technique_settings["project_to_normal"].GetBool()

        super().__init__(name,model,settings)

        fixed_model_parts = self.technique_settings["fixed_model_parts"]
        fixed_model_parts_X = self.technique_settings["fixed_model_parts_X"]
        fixed_model_parts_Y = self.technique_settings["fixed_model_parts_Y"]
        fixed_model_parts_Z = self.technique_settings["fixed_model_parts_Z"]

        if not (fixed_model_parts.size() == fixed_model_parts_X.size() == fixed_model_parts_Y.size() == fixed_model_parts_Z.size()):
            raise RuntimeError("ImplicitVertexMorphing:__init__: fixed_model_parts & fixed_model_parts_X & fixed_model_parts_Y & fixed_model_parts_Z should have the same size")

        for i in range(fixed_model_parts.size()):
            if not fixed_model_parts[i].IsString():
                raise RuntimeError("ImplicitVertexMorphing:__init__: entry {} of 'fixed_model_parts' of control '{}' should be a string .".format(i+1,self.name))
            if not fixed_model_parts_X[i].IsBool():
                raise RuntimeError("ImplicitVertexMorphing:__init__: entry {} of 'fixed_model_parts_X' of control '{}' should be a bool .".format(i+1,self.name))
            if not fixed_model_parts_Y[i].IsBool():
                raise RuntimeError("ImplicitVertexMorphing:__init__: entry {} of 'fixed_model_parts_Y' of control '{}' should be a bool .".format(i+1,self.name))
            if not fixed_model_parts_Z[i].IsBool():
                raise RuntimeError("ImplicitVertexMorphing:__init__: entry {} of 'fixed_model_parts_Z' of control '{}' should be a bool .".format(i+1,self.name))
            fixed_model_part = fixed_model_parts[i].GetString()
            root_fixed_model_part = fixed_model_part.split(".")[0]
            if not self.model.HasModelPart(root_fixed_model_part):
                raise RuntimeError("ImplicitVertexMorphing:__init__: model_part {} in fixed_model_parts of control '{}' does not exist .".format(fixed_model_part,self.name))
        # add vars
        for model_part_name in self.controlling_objects:
            root_model = model_part_name.split(".")[0]
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.HELMHOLTZ_VARS_SHAPE)
            self.model.GetModelPart(root_model).AddNodalSolutionStepVariable(KOA.HELMHOLTZ_SOURCE_SHAPE)

        ## Construct the linear solver
        import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
        self.linear_solvers = []
        root_model_parts = []
        for model_part_name in self.controlling_objects:
            extracted_root_model_part_name = model_part_name.split(".")[0]
            if not extracted_root_model_part_name in root_model_parts:
                root_model_parts.append(extracted_root_model_part_name)
                self.linear_solvers.append(python_linear_solver_factory.ConstructSolver(self.technique_settings["linear_solver_settings"]))


        self.implicit_vertex_morphing = KOA.ImplicitVertexMorphing(self.name,self.model,self.linear_solvers,self.settings)

        # add utils
        self.utils = []
        if self.technique_settings["utilities"].size():
            for itr in range(self.technique_settings["utilities"].size()):
                util_settings = self.technique_settings["utilities"][itr]
                util_type = util_settings["type"].GetString()
                if  util_type== "plane_symmetry" or util_type== "rotational_symmetry":
                    for model_part_name in self.controlling_objects:
                        self.utils.append(KOA.SymmetryUtility(util_settings["name"].GetString(),self.model.GetModelPart(model_part_name),util_settings))

    def Initialize(self):
        super().Initialize()
        self.implicit_vertex_morphing.Initialize()
        for util in self.utils:
            util.Initialize()

        self.penetration_expected = False

    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):
        for util in self.utils:
            util.ApplyOnVectorField(derivative_variable_name)


        self.implicit_vertex_morphing.MapFirstDerivative(derivative_variable_name,mapped_derivative_variable_name)

        if self.penetration_expected:
            for node, distance in zip(self.model.GetModelPart(self.controlling_objects[0]).Nodes, self.distances):
                projection_direction = node.GetSolutionStepValue(KOA.AUXILIARY_FIELD)
                if self.ComputeNodalNorm(projection_direction)>1e-10:
                    shape_update = node.GetSolutionStepValue(mapped_derivative_variable_name)
                    projected_shape_update = projection_direction[0]*shape_update[0] + projection_direction[1]*shape_update[1] + projection_direction[2]*shape_update[2]
                    projected_shape_update_vec = projected_shape_update*projection_direction
                    if distance < self.ref_val and projected_shape_update>0.0:
                        projection_direction /= self.ComputeNodalNorm(projection_direction)
                        projected_shape_update = projection_direction[0]*shape_update[0] + projection_direction[1]*shape_update[1] + projection_direction[2]*shape_update[2]
                        projected_shape_update_vec = projected_shape_update*projection_direction
                        tangent_shape_update = shape_update-projected_shape_update_vec
                        node.SetSolutionStepValue(mapped_derivative_variable_name, tangent_shape_update)
                    elif projected_shape_update>0.0:
                        tangent_shape_update = shape_update-projected_shape_update_vec
                        node.SetSolutionStepValue(mapped_derivative_variable_name, tangent_shape_update)

    def ComputeNodalNorm(self,kratos_nodal_vector):
        return math.sqrt(kratos_nodal_vector[0]*kratos_nodal_vector[0]+kratos_nodal_vector[1]*kratos_nodal_vector[1]+kratos_nodal_vector[2]*kratos_nodal_vector[2])

    def Compute(self):
        self.implicit_vertex_morphing.MapControlUpdate(KOA.D_CX,KOA.D_X)

        try:
            import QuESo_PythonApplication as QuESoApp
            from QuESo_PythonApplication.PyQuESo import PyQuESo
        except ImportError:
            raise Exception("QuESo python library is not available")

        pyqueso = PyQuESo("QUESOParameters.json")

        nodes = QuESoApp.PointVector()
        directions = QuESoApp.PointVector()
        for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
            shape_update = node.GetSolutionStepValue(KOA.D_X)
            nodal_normal = node.GetSolutionStepValue(KM.NORMAL)
            nodal_normal *= -1
            projection = nodal_normal[0] * shape_update[0] + nodal_normal[1] * shape_update[1] + nodal_normal[2] * shape_update[2]
            direction = projection*nodal_normal
            nodes.append( QuESoApp.Point(node.X0, node.Y0, node.Z0) )
            directions.append( QuESoApp.Point(nodal_normal[0], nodal_normal[1], nodal_normal[2]) )

        pyqueso.Run()

        self.ref_val = 20
        self.distances = pyqueso.ClosestDistances(nodes, directions)
        self.penetration_expected = False
        for node, distance in zip(self.model.GetModelPart(self.controlling_objects[0]).Nodes, self.distances):
            if abs(distance) < self.ref_val:
                self.penetration_expected = True
                shape_update = node.GetSolutionStepValue(KOA.D_X)
                nodal_normal = node.GetSolutionStepValue(KM.NORMAL)
                nodal_normal *= -1
                projection = nodal_normal[0] * shape_update[0] + nodal_normal[1] * shape_update[1] + nodal_normal[2] * shape_update[2]
                direction = projection*nodal_normal
                direction /= self.ComputeNodalNorm(direction)
                node.SetSolutionStepValue(KM.KratosGlobals.GetVariable("ADJOINT_DISPLACEMENT"), nodal_normal)
            else:
                node.SetSolutionStepValue(KM.KratosGlobals.GetVariable("ADJOINT_DISPLACEMENT"), [0.0, 0.0, 0.0])

        if self.penetration_expected:

            fixed_model_parts = self.technique_settings["fixed_model_parts"]

            for model_part_name in fixed_model_parts:
                for node in self.model.GetModelPart(model_part_name.GetString()).Nodes:
                    node.Free(KOA.HELMHOLTZ_VARS_SHAPE_X)
                    node.Free(KOA.HELMHOLTZ_VARS_SHAPE_Y)
                    node.Free(KOA.HELMHOLTZ_VARS_SHAPE_Z)
            self.implicit_vertex_morphing.MapControlUpdate(KM.KratosGlobals.GetVariable("ADJOINT_DISPLACEMENT"),KOA.AUXILIARY_FIELD)
            for model_part_name in fixed_model_parts:
                for node in self.model.GetModelPart(model_part_name.GetString()).Nodes:
                    node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_X)
                    node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_Y)
                    node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_Z)


            self.max_norm = 0
            for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
                projection_direction = node.GetSolutionStepValue(KOA.AUXILIARY_FIELD)
                norm = self.ComputeNodalNorm(projection_direction)
                if norm>self.max_norm:
                    self.max_norm = norm

            for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
                projection_direction = node.GetSolutionStepValue(KOA.AUXILIARY_FIELD)
                node.SetSolutionStepValue(KOA.AUXILIARY_FIELD, projection_direction/self.max_norm)

            for node, distance in zip(self.model.GetModelPart(self.controlling_objects[0]).Nodes, self.distances):
                projection_direction = node.GetSolutionStepValue(KOA.AUXILIARY_FIELD)
                if self.ComputeNodalNorm(projection_direction)>1e-10:
                    shape_update = node.GetSolutionStepValue(KOA.D_CX)
                    projected_shape_update = projection_direction[0]*shape_update[0] + projection_direction[1]*shape_update[1] + projection_direction[2]*shape_update[2]
                    projected_shape_update_vec = projected_shape_update*projection_direction
                    if distance < self.ref_val and projected_shape_update>0.0:
                        projection_direction /= self.ComputeNodalNorm(projection_direction)
                        projected_shape_update = projection_direction[0]*shape_update[0] + projection_direction[1]*shape_update[1] + projection_direction[2]*shape_update[2]
                        projected_shape_update_vec = projected_shape_update*projection_direction
                        tangent_shape_update = shape_update-projected_shape_update_vec
                        node.SetSolutionStepValue(KOA.D_CX, tangent_shape_update)
                    elif projected_shape_update>0.0:
                        tangent_shape_update = shape_update-projected_shape_update_vec
                        node.SetSolutionStepValue(KOA.D_CX, tangent_shape_update)


        self.implicit_vertex_morphing.MapControlUpdate(KOA.D_CX,KOA.D_X)

    def Update(self):
        self.implicit_vertex_morphing.Update()

    def GetControllingObjects(self):
        if self.technique_settings["only_design_surface_parameterization"].GetBool():
            return self.controlling_objects
        else:
            root_controlling_names=[]
            for model_part_name in self.controlling_objects:
                root_model = model_part_name.split(".")[0]
                root_controlling_names.append(root_model)
            return root_controlling_names



