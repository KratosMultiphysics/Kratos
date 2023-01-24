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

# ==============================================================================
# --------------------------------------------------------------------------
def ComputeSensitivityHeatmap( design_surface, objectives, constraints,
                               optimization_iteration, mapper, design_variable_name, design_variable_dimension,
                               settings ):

    continuous_response_gradients = ("plane_based_packaging", "mesh_based_packaging")

    if objectives.size() > 1:
        raise RuntimeError("Sensitivity Heatmap only supports one objective function!")

    optimization_utilities = KSO.OptimizationUtilities

    norm_type = settings["norm_type"].GetString()
    sensitivity_weighting = settings["sensitivity_weighting"].GetBool()
    map_sensitivities = settings["mapping"].GetBool()

    # relaxation
    if settings["relaxation_coefficient"].IsString():
        if settings["relaxation_coefficient"].GetString() == "reciprocal":
            relax_coeff = 1 / optimization_iteration
        else:
            raise NameError("The following relaxation coefficient for the sensitivity heatmap is not supported: " + settings["relaxation_coefficient"].GetString())
    else:
        relax_coeff = settings["relaxation_coefficient"].GetDouble()

    # calculate nodal areas for weighting
    if sensitivity_weighting:
        KSO.GeometryUtilities(design_surface).CalculateNodalAreasFromConditions()

    objective_gradient_name = f"DF1D{design_variable_name}"
    objective_is_weigthed = False
    if sensitivity_weighting and \
        not objectives[0]["response_settings"]["response_type"].GetString() in continuous_response_gradients:
            objective_is_weigthed = True

    if objective_is_weigthed:
        WeightResponseGradient(design_surface, objective_gradient_name, design_variable_dimension)
        objective_gradient_name += "_WEIGHTED"
    if map_sensitivities:
        if objective_is_weigthed:
            MapResponseGradient(objective_gradient_name, mapper)
        objective_gradient_name += "_MAPPED"
    ComputeResponseHeatmap(design_surface, objective_gradient_name, optimization_iteration,
                           relax_coeff, design_variable_name, design_variable_dimension)

    if constraints.size() != 0:

        def ___ComputeHeatmapWithNorm(norm_type):
            # normalize objective gradient
            df_dx_normalized = NormalizeResponseGradient(design_surface, objective_gradient_name, norm_type)

            dc_dx_normalized = []
            for itr in range(constraints.size()):
                # constraint gradients
                constraint_gradient_name = f"DC{(itr+1)}D{design_variable_name}"
                constraint_is_weigthed = False
                if sensitivity_weighting and \
                    not constraints[itr]["response_settings"]["response_type"] in continuous_response_gradients:
                        constraint_is_weigthed = True
                if constraint_is_weigthed:
                    WeightResponseGradient(design_surface, constraint_gradient_name, design_variable_dimension)
                    constraint_gradient_name += "_WEIGHTED"
                if map_sensitivities:
                    if constraint_is_weigthed:
                        MapResponseGradient(constraint_gradient_name, mapper)
                    constraint_gradient_name += "_MAPPED"
                # DCiDX individual heatmap
                ComputeResponseHeatmap(design_surface, constraint_gradient_name, optimization_iteration,
                                       relax_coeff, design_variable_name, design_variable_dimension)

                # normalize constraints
                dc_dx_normalized.append(NormalizeResponseGradient(design_surface, constraint_gradient_name, norm_type))

            heat = Kratos.Vector(len(design_surface.Nodes))
            # fill heat map for each node
            for i in range(len(design_surface.Nodes)):
                index = design_variable_dimension*i
                df_dx_i = df_dx_normalized[index:index+design_variable_dimension]
                df_dx_i_norm = cm.Norm2(df_dx_i)

                heat_i = df_dx_i_norm
                for dc_dx in dc_dx_normalized:
                    dc_dx_i = dc_dx[index:index+design_variable_dimension]
                    dc_dx_i_norm = cm.Norm2(dc_dx_i)
                    heat_i = max(heat_i, dc_dx_i_norm)

                heat[i] = heat_i

            heat_map_name = f"HEATMAP_{norm_type}"

            # Heatmap Relaxed
            if optimization_iteration == 1:
                heat_relaxed = heat
            else:
                prev_heat = Kratos.Vector()
                optimization_utilities.AssembleVector(design_surface, prev_heat, Kratos.KratosGlobals.GetVariable(heat_map_name))
                heat_relaxed = Kratos.Vector(len(design_surface.Nodes))
                for i in range(len(design_surface.Nodes)):
                    heat_relaxed[i] = relax_coeff * heat[i] + (1 - relax_coeff) * prev_heat[i]

            optimization_utilities.AssignVectorToVariable(design_surface, heat_relaxed, Kratos.KratosGlobals.GetVariable(heat_map_name))

        if norm_type == "max":
            ___ComputeHeatmapWithNorm(norm_type="MAX")
        elif norm_type == "l2":
            ___ComputeHeatmapWithNorm(norm_type="L2")
        else:
            raise NameError("The following norm type is not supported by the sensitivity heatmap logger (name may be misspelled): " + norm_type)

def ComputeResponseHeatmap( design_surface, response_gradient_name, optimization_iteration,
                            relaxation_coefficient, design_variable_name, design_variable_dimension ):

    optimization_utilities = KSO.OptimizationUtilities

    dg_dx = Kratos.Vector()
    optimization_utilities.AssembleVector(design_surface, dg_dx, Kratos.KratosGlobals.GetVariable(response_gradient_name))

    # DF1DX individual heatmap
    g_name = f"{response_gradient_name[0:4]}"
    heatmap_dgdx_name = f"HEATMAP_{g_name}{design_variable_name}"
    if optimization_iteration == 1:
        heat_dfdx_relaxed = dg_dx
    else:
        prev_heat_dfdx = Kratos.Vector()
        optimization_utilities.AssembleVector(design_surface, prev_heat_dfdx, Kratos.KratosGlobals.GetVariable(heatmap_dgdx_name))
        heat_dfdx_relaxed = Kratos.Vector(design_variable_dimension*len(design_surface.Nodes))
        for i in range(len(design_surface.Nodes)):
            for dim in range(design_variable_dimension):
                index = design_variable_dimension * i + dim
                heat_dfdx_relaxed[index] = relaxation_coefficient * dg_dx[index] + (1 - relaxation_coefficient) * prev_heat_dfdx[index]

    optimization_utilities.AssignVectorToVariable(design_surface, heat_dfdx_relaxed, Kratos.KratosGlobals.GetVariable(heatmap_dgdx_name))

def NormalizeResponseGradient( design_surface, response_gradient_name, norm_type ):

    optimization_utilities = KSO.OptimizationUtilities

    dg_dx = Kratos.Vector()
    optimization_utilities.AssembleVector(design_surface, dg_dx, Kratos.KratosGlobals.GetVariable(response_gradient_name))

    if norm_type == "MAX":
        dg_dx_norm = cm.NormInf3D(dg_dx)
    elif norm_type == "L2":
        dg_dx_norm = cm.Norm2(dg_dx)
    else:
        raise RuntimeError("Sensitivity Heatmap only supports 'max' or 'l2' norm!")

    if dg_dx_norm != 0.0:
        dg_dx_normalized = (1/dg_dx_norm) * dg_dx
    else:
        dg_dx_normalized = 0 * dg_dx

    return dg_dx_normalized

def WeightResponseGradient( design_surface, response_gradient_name,
                            design_variable_dimension ):

    optimization_utilities = KSO.OptimizationUtilities

    nodal_area = Kratos.Vector()
    optimization_utilities.AssembleVector(design_surface, nodal_area, Kratos.KratosGlobals.GetVariable("NODAL_AREA"))

    dg_dx = Kratos.Vector()
    optimization_utilities.AssembleVector(design_surface, dg_dx, Kratos.KratosGlobals.GetVariable(response_gradient_name))

    # weight response gradients with nodal area
    dg_dx_weighted = Kratos.Vector(design_variable_dimension*len(design_surface.Nodes))
    for i in range(len(design_surface.Nodes)):
        for dim in range(design_variable_dimension):
            index = design_variable_dimension * i + dim
            if nodal_area[i] > 1e-8:
                dg_dx_weighted[index] = dg_dx[index] / nodal_area[i]
            else:
                dg_dx_weighted[index] = dg_dx[index]

    response_gradient_weighted_name = f"{response_gradient_name}_WEIGHTED"
    optimization_utilities.AssignVectorToVariable(design_surface, dg_dx_weighted,
                                                  Kratos.KratosGlobals.GetVariable(response_gradient_weighted_name))

def MapResponseGradient( response_gradient_name , mapper):

    response_gradient_mapped_name = f"{response_gradient_name}_MAPPED"
    mapper.InverseMap(Kratos.KratosGlobals.GetVariable(response_gradient_name),
                      Kratos.KratosGlobals.GetVariable(response_gradient_mapped_name))


# ==============================================================================
