# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Schmölz David, https://github.com/dschmoelz
#
# ==============================================================================

# Kratos Core and Apps
import KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.ShapeOptimizationApplication.utilities import custom_math as cm
# Import additional libraries

# some response are continuously formulated
continuous_response_gradients = ("plane_based_packaging", "mesh_based_packaging")

# ==============================================================================
# --------------------------------------------------------------------------
class SensitivityHeatmapCalculator(object):

    # --------------------------------------------------------------------------
    def __init__(self, design_surface, objectives, constraints, settings):
        self.design_surface = design_surface
        self.objectives = objectives
        self.constraints = constraints

        if objectives.size() > 1:
            raise RuntimeError("Sensitivity Heatmap only supports one objective function!")

        self.settings = settings
        self.norm_type = settings["norm_type"].GetString()
        if self.norm_type not in ("max", "l2"):
            raise NameError("""The following norm type is not supported by the
                            sensitivity heatmap logger (name may be misspelled): """ + self.norm_type)

        self.sensitivity_weighting = settings["sensitivity_weighting"].GetBool()
        self.map_sensitivities = settings["mapping"].GetBool()
        self.relaxation_method = settings["relaxation_method"].GetString()
        if self.relaxation_method not in ("reciprocal", "constant"):
            raise NameError("""The following relaxation method is not supported by the
                            sensitivity heatmap logger (name may be misspelled): """ + self.relaxation_method)

        self.design_variable_name = settings["design_variable_name"].GetString()
        self.design_variable_dimension = settings["design_variable_dimension"].GetInt()

        self.optimization_utilities = KSO.OptimizationUtilities

    # --------------------------------------------------------------------------
    def ComputeHeatmaps(self, optimization_iteration, mapper):

        # calculate nodal areas for weighting
        if self.sensitivity_weighting:
            KSO.GeometryUtilities(self.design_surface).CalculateNodalAreasFromConditions()

        self.__ComputeIndividualHeatmaps(optimization_iteration, mapper)
        if self.constraints.size() != 0:
            self.__ComputeAggregatedHeatmap(optimization_iteration, mapper)

    # --------------------------------------------------------------------------
    def __ComputeIndividualHeatmaps(self, optimization_iteration, mapper):

        objective_gradient_name = f"DF1D{self.design_variable_name}"
        objective_type = self.objectives[0]["response_settings"]["response_type"]
        self.__ComputeResponseHeatmap(objective_gradient_name, objective_type, optimization_iteration, mapper)

        if self.constraints.size() != 0:
            for itr in range(self.constraints.size()):
                constraint_gradient_name = f"DC{(itr+1)}D{self.design_variable_name}"
                constraint_type = self.constraints[itr]["response_settings"]["response_type"]
                self.__ComputeResponseHeatmap(constraint_gradient_name, constraint_type, optimization_iteration, mapper)

    # --------------------------------------------------------------------------
    def __ComputeAggregatedHeatmap(self, optimization_iteration, mapper):

        # normalize objective gradient
        objective_gradient_name = f"DF1D{self.design_variable_name}"
        objective_type = self.objectives[0]["response_settings"]["response_type"]
        df_dx_normalized = self.__NormalizeResponseGradient(objective_gradient_name, objective_type)

        # normalize constraints
        dc_dx_normalized = []
        for itr in range(self.constraints.size()):
            constraint_gradient_name = f"DC{(itr+1)}D{self.design_variable_name}"
            constraint_type = self.constraints[itr]["response_settings"]["response_type"]
            dc_dx_normalized.append(self.__NormalizeResponseGradient(constraint_gradient_name, constraint_type))

        # fill heat map for each node
        heat = Kratos.Vector(len(self.design_surface.Nodes))
        for i in range(len(self.design_surface.Nodes)):
            index = self.design_variable_dimension*i
            df_dx_i = df_dx_normalized[index:index+self.design_variable_dimension]
            df_dx_i_norm = cm.Norm2(df_dx_i)

            heat_i = df_dx_i_norm
            for dc_dx in dc_dx_normalized:
                dc_dx_i = dc_dx[index:index+self.design_variable_dimension]
                dc_dx_i_norm = cm.Norm2(dc_dx_i)
                heat_i = max(heat_i, dc_dx_i_norm)

            heat[i] = heat_i

        if self.norm_type == "max":
            heat_map_name = "HEATMAP_MAX"
        elif self.norm_type == "l2":
            heat_map_name = "HEATMAP_L2"

        # relax heatmap
        relaxation_coefficient = self.__GetRelaxationCoefficient(optimization_iteration)
        if optimization_iteration == 1:
            heat_relaxed = heat
        else:
            prev_heat = Kratos.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, prev_heat, Kratos.KratosGlobals.GetVariable(heat_map_name))
            heat_relaxed = Kratos.Vector(len(self.design_surface.Nodes))
            for i in range(len(self.design_surface.Nodes)):
                heat_relaxed[i] = relaxation_coefficient * heat[i] + (1 - relaxation_coefficient) * prev_heat[i]

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, heat_relaxed, Kratos.KratosGlobals.GetVariable(heat_map_name))

    # --------------------------------------------------------------------------
    def __ComputeResponseHeatmap( self, response_gradient_name, response_type, optimization_iteration, mapper):

        response_is_weighted = False
        if self.sensitivity_weighting and \
            response_type not in continuous_response_gradients:
            response_is_weighted = True

        if response_is_weighted:
            self.__WeightResponseGradient(response_gradient_name)
            response_gradient_name += "_WEIGHTED"
        if self.map_sensitivities:
            # reponse has been already mapped by the optimization algorithm if it is not weighted
            if response_is_weighted:
                self.__MapResponseGradient(response_gradient_name, mapper)
            response_gradient_name += "_MAPPED"

        self.__RelaxResponseHeatmap(response_gradient_name, optimization_iteration)

    # --------------------------------------------------------------------------
    def __WeightResponseGradient( self, response_gradient_name ):

        nodal_area = Kratos.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nodal_area, Kratos.KratosGlobals.GetVariable("NODAL_AREA"))

        dg_dx = Kratos.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, dg_dx, Kratos.KratosGlobals.GetVariable(response_gradient_name))

        # weight response gradients with nodal area
        dg_dx_weighted = Kratos.Vector(self.design_variable_dimension*len(self.design_surface.Nodes))
        for i in range(len(self.design_surface.Nodes)):
            for dim in range(self.design_variable_dimension):
                index = self.design_variable_dimension * i + dim
                if nodal_area[i] > 1e-8:
                    dg_dx_weighted[index] = dg_dx[index] / nodal_area[i]
                else:
                    dg_dx_weighted[index] = dg_dx[index]

        response_gradient_weighted_name = f"{response_gradient_name}_WEIGHTED"
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, dg_dx_weighted,
                                                           Kratos.KratosGlobals.GetVariable(response_gradient_weighted_name))

    # --------------------------------------------------------------------------
    def __MapResponseGradient( self, response_gradient_name, mapper ):

        response_gradient_mapped_name = f"{response_gradient_name}_MAPPED"
        mapper.InverseMap(Kratos.KratosGlobals.GetVariable(response_gradient_name),
                               Kratos.KratosGlobals.GetVariable(response_gradient_mapped_name))

    # --------------------------------------------------------------------------
    def __RelaxResponseHeatmap( self, response_gradient_name, optimization_iteration ):

        relaxation_coefficient = self.__GetRelaxationCoefficient(optimization_iteration)

        dg_dx = Kratos.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, dg_dx, Kratos.KratosGlobals.GetVariable(response_gradient_name))

        g_name = f"{response_gradient_name[0:4]}"
        heatmap_dgdx_name = f"HEATMAP_{g_name}{self.design_variable_name}"
        if optimization_iteration == 1:
            heat_dfdx_relaxed = dg_dx
        else:
            prev_heat_dfdx = Kratos.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, prev_heat_dfdx, Kratos.KratosGlobals.GetVariable(heatmap_dgdx_name))
            heat_dfdx_relaxed = Kratos.Vector(self.design_variable_dimension*len(self.design_surface.Nodes))
            for i in range(len(self.design_surface.Nodes)):
                for dim in range(self.design_variable_dimension):
                    index = self.design_variable_dimension * i + dim
                    heat_dfdx_relaxed[index] = relaxation_coefficient * dg_dx[index] + (1 - relaxation_coefficient) * prev_heat_dfdx[index]

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, heat_dfdx_relaxed, Kratos.KratosGlobals.GetVariable(heatmap_dgdx_name))

    # --------------------------------------------------------------------------
    def __NormalizeResponseGradient( self, response_gradient_name, response_type ):

        response_is_weighted = False
        if self.sensitivity_weighting and \
            response_type not in continuous_response_gradients:
            response_is_weighted = True

        if response_is_weighted:
            response_gradient_name += "_WEIGHTED"
        if self.map_sensitivities:
            response_gradient_name += "_MAPPED"

        optimization_utilities = KSO.OptimizationUtilities

        dg_dx = Kratos.Vector()
        optimization_utilities.AssembleVector(self.design_surface, dg_dx, Kratos.KratosGlobals.GetVariable(response_gradient_name))

        if self.norm_type == "max":
            dg_dx_norm = cm.NormInf3D(dg_dx)
        elif self.norm_type == "l2":
            dg_dx_norm = cm.Norm2(dg_dx)
        else:
            raise RuntimeError("Sensitivity Heatmap only supports 'max' or 'l2' norm!")

        if dg_dx_norm != 0.0:
            dg_dx_normalized = (1/dg_dx_norm) * dg_dx
        else:
            dg_dx_normalized = 0 * dg_dx

        return dg_dx_normalized

    # --------------------------------------------------------------------------
    def __GetRelaxationCoefficient( self, optimization_iteration ):

        if self.relaxation_method == "reciprocal":
            return 1 / optimization_iteration
        elif self.relaxation_method == "constant":
            return self.settings["relaxation_coefficient"].GetDouble()

# ==============================================================================
