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
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.ShapeOptimizationApplication.utilities import custom_math as cm
# Import additional libraries

# ==============================================================================
# --------------------------------------------------------------------------
def ComputeSensitivityHeatmap( design_surface, objectives, constraints,
                               optimization_iteration, design_variable_name="X", design_variable_dimension=3 ):

        if objectives.size() > 1:
            raise RuntimeError("Sensitivity Heatmap only supports one objective function!")

        optimization_utilities = KSO.OptimizationUtilities

        # reciprocal relaxation
        relax_coeff = 1 / optimization_iteration

        objective_gradient_name = f"DF1D{design_variable_name}_MAPPED"
        ComputeResponseHeatmap(design_surface, objective_gradient_name, optimization_iteration, design_variable_name, design_variable_dimension)

        if constraints.size() != 0:

            def ___ComputeHeatmapWithNorm(norm_type):
                # normalize objective gradient
                df_dx_normalized = NormalizeResponseGradient(design_surface, objective_gradient_name, norm_type)

                dc_dx_normalized = []
                for itr, constraint in enumerate(constraints):
                    # read constraint gradients
                    constraint_gradient_name = f"DC{(itr+1)}D{design_variable_name}_MAPPED"

                    # DCiDX individual heatmap
                    ComputeResponseHeatmap(design_surface, constraint_gradient_name, optimization_iteration, design_variable_name, design_variable_dimension)

                    # normalize constraints
                    dc_dx_normalized.append(NormalizeResponseGradient(design_surface, constraint_gradient_name, norm_type))

                heat = KM.Vector(len(design_surface.Nodes))
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
                    prev_heat = KM.Vector()
                    optimization_utilities.AssembleVector(design_surface, prev_heat, KM.KratosGlobals.GetVariable(heat_map_name))
                    heat_relaxed = KM.Vector(len(design_surface.Nodes))
                    for i in range(len(design_surface.Nodes)):
                        heat_relaxed[i] = relax_coeff * heat[i] + (1 - relax_coeff) * prev_heat[i]

                optimization_utilities.AssignVectorToVariable(design_surface, heat_relaxed, KM.KratosGlobals.GetVariable(heat_map_name))

            ___ComputeHeatmapWithNorm(norm_type="MAX")
            ___ComputeHeatmapWithNorm(norm_type="L2")

def ComputeResponseHeatmap( design_surface, response_gradient_name, optimization_iteration,
                            design_variable_name="X", design_variable_dimension=3 ):

    optimization_utilities = KSO.OptimizationUtilities

    # reciprocal relaxation
    relax_coeff = 1 / optimization_iteration

    dg_dx = KM.Vector()
    optimization_utilities.AssembleVector(design_surface, dg_dx, KM.KratosGlobals.GetVariable(response_gradient_name))

    # DF1DX individual heatmap
    g_name = f"{response_gradient_name[0:4]}"
    heatmap_dgdx_name = f"HEATMAP_{g_name}{design_variable_name}"
    if optimization_iteration == 1:
        heat_dfdx_relaxed = dg_dx
    else:
        prev_heat_dfdx = KM.Vector()
        optimization_utilities.AssembleVector(design_surface, prev_heat_dfdx, KM.KratosGlobals.GetVariable(heatmap_dgdx_name))
        heat_dfdx_relaxed = KM.Vector(design_variable_dimension*len(design_surface.Nodes))
        for i in range(len(design_surface.Nodes)):
            for dim in range(design_variable_dimension):
                index = design_variable_dimension * i + dim
                heat_dfdx_relaxed[index] = relax_coeff * dg_dx[index] + (1 - relax_coeff) * prev_heat_dfdx[index]

    optimization_utilities.AssignVectorToVariable(design_surface, heat_dfdx_relaxed, KM.KratosGlobals.GetVariable(heatmap_dgdx_name))

def NormalizeResponseGradient( design_surface, response_gradient_name, norm_type ):

    optimization_utilities = KSO.OptimizationUtilities

    dg_dx = KM.Vector()
    optimization_utilities.AssembleVector(design_surface, dg_dx, KM.KratosGlobals.GetVariable(response_gradient_name))

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

# ==============================================================================
