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

    def Initialize(self, opt_parameters):
        super().Initialize(opt_parameters)
        self.implicit_vertex_morphing.Initialize()
        for util in self.utils:
            util.Initialize()

        self.Update()

    def NormalizeField(self, variable_name):
        max_value = 0.0
        for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
            values = node.GetSolutionStepValue(variable_name)
            max_value = max(max_value, self.ComputeNodalNorm(values))

        if(max_value > 1e-10):
            for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
                normalized_values = 1.0/max_value * node.GetSolutionStepValue(variable_name)
                node.SetSolutionStepValue(variable_name, normalized_values)

    def GetActiveSensitivityDirection(self, sensitivity_variable_name):
        for objective in self.opt_parameters["objectives"]:
            if( objective["control_gradient_names"][0].GetString() == sensitivity_variable_name  ):
                return -1
        for constraint in self.opt_parameters["constraints"]:
            if( constraint["control_gradient_names"][0].GetString() == sensitivity_variable_name  ):
                if( "equality" in constraint["type"].GetString() ):
                    current_value = constraint["value"].GetDouble()
                    ref_value = constraint["ref_value"].GetDouble()
                    if( current_value > ref_value ):
                        return -1
                    else:
                        return 1

    def MapFirstDerivative(self,derivative_variable_name,mapped_derivative_variable_name):

        for util in self.utils:
            util.ApplyOnVectorField(derivative_variable_name)

        self.implicit_vertex_morphing.MapFirstDerivative(derivative_variable_name, mapped_derivative_variable_name)

        #direction = self.GetActiveSensitivityDirection(mapped_derivative_variable_name.Name())
        self.RemoveActiveFilteredIntersectionNormals(mapped_derivative_variable_name)
        self.RemoveActiveFilteredBCNormals(mapped_derivative_variable_name)

    def ComputeNodalNorm(self,kratos_nodal_vector):
        return math.sqrt(kratos_nodal_vector[0]*kratos_nodal_vector[0]+kratos_nodal_vector[1]*kratos_nodal_vector[1]+kratos_nodal_vector[2]*kratos_nodal_vector[2])

    def MapField(self, variable_raw, variable_filtered, dirichlet_is_active=True):
        if dirichlet_is_active:
            self.implicit_vertex_morphing.MapControlUpdate(variable_raw, variable_filtered)
        else:
            fixed_model_parts = self.technique_settings["fixed_model_parts"]
            fixed_model_parts_X = self.technique_settings["fixed_model_parts_X"]
            fixed_model_parts_Y = self.technique_settings["fixed_model_parts_Y"]
            fixed_model_parts_Z = self.technique_settings["fixed_model_parts_Z"]

            for model_part_name, fixed_x, fixed_y, fixed_z in zip(fixed_model_parts, fixed_model_parts_X, fixed_model_parts_Y, fixed_model_parts_Z):
                for node in self.model.GetModelPart(model_part_name.GetString()).Nodes:
                    if(fixed_x):
                        node.Free(KOA.HELMHOLTZ_VARS_SHAPE_X)
                    if(fixed_y):
                        node.Free(KOA.HELMHOLTZ_VARS_SHAPE_Y)
                    if(fixed_z):
                        node.Free(KOA.HELMHOLTZ_VARS_SHAPE_Z)

            self.implicit_vertex_morphing.MapControlUpdate(variable_raw, variable_filtered)

            for model_part_name, fixed_x, fixed_y, fixed_z in zip(fixed_model_parts, fixed_model_parts_X, fixed_model_parts_Y, fixed_model_parts_Z):
                for node in self.model.GetModelPart(model_part_name.GetString()).Nodes:
                    if(fixed_x):
                        node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_X)
                    if(fixed_y):
                        node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_Y)
                    if(fixed_z):
                        node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_Z)
                    if(fixed_x or fixed_y or fixed_z):
                        node.SetSolutionStepValue(variable_filtered, [0.0, 0.0, 0.0])


    def ComputeActiveIntersectionNormal(self, min_distance, positive_dir=True):
        try:
            import QuESo_PythonApplication as QuESoApp
            from QuESo_PythonApplication.PyQuESo import PyQuESo
        except ImportError:
            raise Exception("QuESo python library is not available")

        pyqueso = PyQuESo("QUESOParameters_vertical.json")

        nodes = QuESoApp.PointVector()
        directions = QuESoApp.PointVector()
        for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
            nodal_normal = node.GetSolutionStepValue(KM.NORMAL)
            if positive_dir:
                nodal_normal *= -1
            nodes.append( QuESoApp.Point(node.X0, node.Y0, node.Z0) )
            directions.append( QuESoApp.Point(nodal_normal[0], nodal_normal[1], nodal_normal[2]) )

        pyqueso.Run(self.model.GetModelPart(self.controlling_objects[0]))
        self.distances = pyqueso.ClosestDistances(nodes, directions)

        if(positive_dir):
            variable = KOA.ACTIVE_POSITIVE_INTERSECTION_NORMAL
        else:
            variable = KOA.ACTIVE_NEGATIVE_INTERSECTION_NORMAL

        intersection_active = False
        for node, distance in zip(self.model.GetModelPart(self.controlling_objects[0]).Nodes, self.distances):
            nodal_normal = node.GetSolutionStepValue(KM.NORMAL)
            if positive_dir:
                nodal_normal *= -1

            if (abs(distance) > min_distance and self.nodal_max_control_update>abs(distance-min_distance)) or (abs(distance) < min_distance ):
                intersection_active = True
                node.SetSolutionStepValue(variable, nodal_normal)
            else:
                node.SetSolutionStepValue(variable, [0.0, 0.0, 0.0])

        return intersection_active

    def ComputeConstructionSpaceConstraint(self):
        try:
            import QuESo_PythonApplication as QuESoApp
            from QuESo_PythonApplication.PyQuESo import PyQuESo
        except ImportError:
            raise Exception("QuESo python library is not available")

        pyqueso = PyQuESo("QUESOParameters_con_space.json")

        nodes = QuESoApp.PointVector()
        for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
            nodal_normal = self.nodal_max_control_update*node.GetSolutionStepValue(KM.NORMAL)
            nodes.append( QuESoApp.Point(node.X0 + nodal_normal[0], node.Y0 + nodal_normal[1], node.Z0 + nodal_normal[2]) )

        pyqueso.Run()
        are_inside = pyqueso.IsInside(nodes)

        variable = KOA.ACTIVE_NEGATIVE_INTERSECTION_NORMAL

        intersection_active = False
        for node, is_inside in zip(self.model.GetModelPart(self.controlling_objects[0]).Nodes, are_inside):
            nodal_normal = node.GetSolutionStepValue(KM.NORMAL)
            if not is_inside:
                intersection_active = True
                node.SetSolutionStepValue(variable, nodal_normal)

        return intersection_active

    def RemoveActiveFilteredIntersectionNormals(self, field_variable_name):

        pos_active = self.positive_intersection_normals_active
        neg_active = self.negative_intersection_normals_active or self.construction_space_constraint_active

        for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
            nodal_field_variable = node.GetSolutionStepValue(field_variable_name)

            C = np.empty((0, 3))
            variable_vector = np.zeros(3)

            variable_vector[0] = nodal_field_variable[0]
            variable_vector[1] = nodal_field_variable[1]
            variable_vector[2] = nodal_field_variable[2]

            if pos_active or neg_active:
                raw_value_pos = node.GetSolutionStepValue(KOA.ACTIVE_POSITIVE_INTERSECTION_NORMAL)
                raw_value_neg = node.GetSolutionStepValue(KOA.ACTIVE_NEGATIVE_INTERSECTION_NORMAL)
                if self.ComputeNodalNorm(raw_value_pos)>1e-10 or self.ComputeNodalNorm(raw_value_neg)>1e-10:
                    remove_direction = node.GetSolutionStepValue(KM.NORMAL)
                    C = np.append(C, np.array([remove_direction]), axis=0)

            if pos_active:
                avg_value_pos = node.GetSolutionStepValue(KOA.FILTERED_ACTIVE_POSITIVE_INTERSECTION_NORMAL)
                norm = self.ComputeNodalNorm(avg_value_pos)
                if(norm > 0.01):
                    avg_value_pos /= norm
                    C = np.append(C, np.array([avg_value_pos]), axis=0)

            if neg_active:
                avg_value_neg = node.GetSolutionStepValue(KOA.FILTERED_ACTIVE_NEGATIVE_INTERSECTION_NORMAL)
                norm = self.ComputeNodalNorm(avg_value_neg)
                if(norm > 0.01):
                    avg_value_neg /= norm
                    C = np.append(C, np.array([avg_value_neg]), axis=0)

            if C.shape[0] > 0:
                inverse = np.linalg.inv( np.matmul(C, (np.transpose(C))) )
                A = np.matmul( np.transpose(C), inverse)
                B = np.matmul( C, variable_vector )
                C = A @ B
                new_field = variable_vector - C
                node.SetSolutionStepValue(field_variable_name, [new_field[0], new_field[1], new_field[2]])



    def ComputeFilteredActiveBCNormals(self):
        fixed_model_parts = self.technique_settings["fixed_model_parts"]
        fixed_model_parts_X = self.technique_settings["fixed_model_parts_X"]
        fixed_model_parts_Y = self.technique_settings["fixed_model_parts_Y"]
        fixed_model_parts_Z = self.technique_settings["fixed_model_parts_Z"]

        for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
            node.SetSolutionStepValue(KOA.ACTIVE_BC_NORMAL, [0.0, 0.0, 0.0])
            node.SetSolutionStepValue(KOA.ADJACENT_TRIANGLE_COUNT, 0)
            node.SetSolutionStepValue(KOA.FILTERED_ACTIVE_BC_NORMAL, [0.0, 0.0, 0.0])

        for el in self.model.GetModelPart(self.controlling_objects[0]).Elements:
            count_fixed = 0
            count_not_fixed = 0
            for node in el.GetGeometry():
                if( node.IsFixed(KOA.HELMHOLTZ_VARS_SHAPE_X) ):
                    count_fixed += 1
                else:
                    count_not_fixed += 1
            #if( count_fixed != 3 and count_fixed > 0):
            for node in el.GetGeometry():
                if( node.IsFixed(KOA.HELMHOLTZ_VARS_SHAPE_X) ):
                    node.SetSolutionStepValue(KOA.ADJACENT_TRIANGLE_COUNT, 1)
                    node.SetSolutionStepValue(KOA.ACTIVE_BC_NORMAL, [1.0, 0.0, 0.0])

        for model_part_name, fixed_x, fixed_y, fixed_z in zip(fixed_model_parts, fixed_model_parts_X, fixed_model_parts_Y, fixed_model_parts_Z):
            for node in self.model.GetModelPart(model_part_name.GetString()).Nodes:
                if(fixed_x):
                    node.Free(KOA.HELMHOLTZ_VARS_SHAPE_X)
                if(fixed_y):
                    node.Free(KOA.HELMHOLTZ_VARS_SHAPE_Y)
                if(fixed_z):
                    node.Free(KOA.HELMHOLTZ_VARS_SHAPE_Z)

        self.implicit_vertex_morphing.SetFilterRadius(1)
        self.implicit_vertex_morphing.MapControlUpdate(KOA.ACTIVE_BC_NORMAL, KOA.FILTERED_ACTIVE_BC_NORMAL)
        self.implicit_vertex_morphing.SetFilterRadius(5)

        for model_part_name, fixed_x, fixed_y, fixed_z in zip(fixed_model_parts, fixed_model_parts_X, fixed_model_parts_Y, fixed_model_parts_Z):
            for node in self.model.GetModelPart(model_part_name.GetString()).Nodes:
                if(fixed_x):
                    node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_X)
                if(fixed_y):
                    node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_Y)
                if(fixed_z):
                    node.Fix(KOA.HELMHOLTZ_VARS_SHAPE_Z)
                if(fixed_x or fixed_y or fixed_z):
                    node.SetSolutionStepValue(KOA.FILTERED_ACTIVE_BC_NORMAL, [0.0, 0.0, 0.0])

        self.NormalizeField(KOA.FILTERED_ACTIVE_BC_NORMAL)


    def RemoveActiveFilteredBCNormals(self, field_variable):
        for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
            filtered_active_bc_normal = node.GetSolutionStepValue(KOA.FILTERED_ACTIVE_BC_NORMAL)
            norm_filtered_active_bc_normal = self.ComputeNodalNorm(filtered_active_bc_normal)
            if norm_filtered_active_bc_normal>0.05:
                node.SetSolutionStepValue(field_variable, [0.0, 0.0, 0.0])

    def ComputeConstrainingFields(self):
        min_distance = 1

        self.nodal_max_control_update = 0.1
        # for node in self.model.GetModelPart(self.controlling_objects[0]).Nodes:
        #     control_update = node.GetSolutionStepValue(KOA.D_X)
        #     self.nodal_max_control_update = max(self.nodal_max_control_update, self.ComputeNodalNorm(control_update))

        self.positive_intersection_normals_active = self.ComputeActiveIntersectionNormal(min_distance,  True)
        self.negative_intersection_normals_active = self.ComputeActiveIntersectionNormal(min_distance, False)
        self.construction_space_constraint_active = self.ComputeConstructionSpaceConstraint()

        active_bcs = True
        if self.positive_intersection_normals_active:
            self.MapField(KOA.ACTIVE_POSITIVE_INTERSECTION_NORMAL, KOA.FILTERED_ACTIVE_POSITIVE_INTERSECTION_NORMAL, active_bcs)
        if self.negative_intersection_normals_active or self.construction_space_constraint_active:
            self.MapField(KOA.ACTIVE_NEGATIVE_INTERSECTION_NORMAL, KOA.FILTERED_ACTIVE_NEGATIVE_INTERSECTION_NORMAL, active_bcs)

        self.ComputeFilteredActiveBCNormals()

    def Compute(self):

        self.RemoveActiveFilteredIntersectionNormals(KOA.D_CX)
        self.RemoveActiveFilteredBCNormals(KOA.D_CX)
        active_bcs = True
        self.MapField(KOA.D_CX, KOA.D_X, active_bcs)

    def Update(self):
        self.implicit_vertex_morphing.Update()

        self.ComputeConstrainingFields()

    def GetControllingObjects(self):
        if self.technique_settings["only_design_surface_parameterization"].GetBool():
            return self.controlling_objects
        else:
            root_controlling_names=[]
            for model_part_name in self.controlling_objects:
                root_model = model_part_name.split(".")[0]
                root_controlling_names.append(root_model)
            return root_controlling_names



