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
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteListToNodalVariable, ReadNodalVariableToList
from KratosMultiphysics.ShapeOptimizationApplication.utilities import custom_math as cm
# Import additional libraries

# ==============================================================================
# --------------------------------------------------------------------------
def ComputeSensitivityHeatmap( design_surface, objectives, constraints, constraint_gradient_variables,
                               optimization_iteration, design_variable_name="X", design_variable_dimension=3 ):

        if objectives.size() > 1:
            raise RuntimeError("Sensitivity Heatmap only supports one objective function!")

        optimization_utilities = KSO.OptimizationUtilities

        # reciprocal relaxation
        relax_coeff = 1 / optimization_iteration

        objective_gradient_name = f"DF1D{design_variable_name}_MAPPED"
        df_dx = ReadNodalVariableToList(design_surface, KM.KratosGlobals.GetVariable(objective_gradient_name), dimension=design_variable_dimension)

        # DF1DX individual heatmap
        heatmap_dfdx_name = f"HEATMAP_DF1D{design_variable_name}"
        if optimization_iteration == 1:
            heat_dfdx_relaxed = df_dx
        else:
            prev_heat_dfdx = KM.Vector()
            optimization_utilities.AssembleVector(design_surface, prev_heat_dfdx, KM.KratosGlobals.GetVariable(heatmap_dfdx_name))
            heat_dfdx_relaxed = []
            for i in range(len(design_surface.Nodes)):
                for dim in range(design_variable_dimension):
                    heat_dfdx_relaxed.append(relax_coeff * df_dx[design_variable_dimension*i+dim] + (1 - relax_coeff) * prev_heat_dfdx[design_variable_dimension*i+dim])

        WriteListToNodalVariable(heat_dfdx_relaxed, design_surface, KM.KratosGlobals.GetVariable(heatmap_dfdx_name), dimension=design_variable_dimension)

        if constraints.size() != 0:

            def ___ComputeHeatmapWithNorm(norm_type):
                heat = []

                # normalize objective gradient
                if norm_type == "MAX":
                    df_dx_norm = cm.NormInf3D(df_dx)
                elif norm_type == "L2":
                    df_dx_norm = cm.Norm2(df_dx)
                else:
                    raise RuntimeError("Sensitivity Heatmap only supports 'max' or 'l2' norm!")

                if df_dx_norm != 0.0:
                    df_dx_normalized = cm.ScalarVectorProduct(1/df_dx_norm, df_dx)
                else:
                    df_dx_normalized = [0] * len(df_dx)

                dc_dx_normalized = {}
                for itr, constraint in enumerate(constraints):
                    # read constraint gradients
                    con_id = constraint["identifier"].GetString()
                    constraint_gradient_name = f"DC{(itr+1)}D{design_variable_name}_MAPPED"

                    dci_dx = ReadNodalVariableToList(design_surface, KM.KratosGlobals.GetVariable(constraint_gradient_name), dimension=design_variable_dimension)

                    # normalize constraints
                    if norm_type == "MAX":
                        dci_dx_norm = cm.NormInf3D(dci_dx)
                    elif norm_type == "L2":
                        dci_dx_norm = cm.Norm2(dci_dx)

                    if dci_dx_norm != 0.0:
                        dc_dx_normalized.update({con_id : cm.ScalarVectorProduct(1/dci_dx_norm, dci_dx)})
                    else:
                        dc_dx_normalized.update({con_id : [0] * len(dci_dx)})

                    # DCiDX individual heatmap
                    heatmap_dcidx_name = f"HEATMAP_DC{itr+1}D{design_variable_name}"
                    if optimization_iteration == 1:
                        heat_dcidx_relaxed = dci_dx
                    else:
                        prev_heat_dcidx = KM.Vector()
                        optimization_utilities.AssembleVector(design_surface, prev_heat_dcidx, KM.KratosGlobals.GetVariable(heatmap_dcidx_name))
                        heat_dcidx_relaxed = []
                        for i in range(len(design_surface.Nodes)):
                            for dim in range(design_variable_dimension):
                                heat_dcidx_relaxed.append(relax_coeff * dci_dx[design_variable_dimension*i+dim] + (1 - relax_coeff) * prev_heat_dcidx[design_variable_dimension*i+dim])

                    WriteListToNodalVariable(heat_dcidx_relaxed, design_surface, KM.KratosGlobals.GetVariable(heatmap_dcidx_name), dimension=design_variable_dimension)

                # fill heat map for each node
                for i in range(len(design_surface.Nodes)):
                    df_dx_i = df_dx_normalized[design_variable_dimension*i:design_variable_dimension*i+design_variable_dimension]
                    df_dx_i_norm = cm.Norm2(df_dx_i)

                    heat_i = df_dx_i_norm
                    for dc_dx in dc_dx_normalized.values():
                        dc_dx_i = dc_dx[design_variable_dimension*i:design_variable_dimension*i+design_variable_dimension]
                        dc_dx_i_norm = cm.Norm2(dc_dx_i)
                        heat_i = max(heat_i, dc_dx_i_norm)

                    heat.append(heat_i)

                heat_map_name = f"HEATMAP_{norm_type}"


                # Heatmap Relaxed
                if optimization_iteration == 1:
                    heat_relaxed = heat
                else:
                    prev_heat = KM.Vector()
                    optimization_utilities.AssembleVector(design_surface, prev_heat, KM.KratosGlobals.GetVariable(heat_map_name))
                    heat_relaxed = []
                    for i in range(len(design_surface.Nodes)):
                        heat_relaxed.append(relax_coeff * heat[i] + (1 - relax_coeff) * prev_heat[i])

                WriteListToNodalVariable(heat_relaxed, design_surface, KM.KratosGlobals.GetVariable(heat_map_name), 1)

            ___ComputeHeatmapWithNorm(norm_type="MAX")
            ___ComputeHeatmapWithNorm(norm_type="L2")

# ==============================================================================
