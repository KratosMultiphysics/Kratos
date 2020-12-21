# Import Python libraries
import numpy as np

# Import Kratos
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication
try:
    import KratosMultiphysics.FluidDynamicsApplication
    have_fluid_dynamics_application = True
except:
    have_fluid_dynamics_application = False

from KratosMultiphysics.MultilevelMonteCarloApplication.tools import ParametersWrapper


class AdaptiveRefinement(object):
    """
    Class managing the adaptive refinement process, called when executing Multilevel Monte Carlo algorithms.
    This class handles the call to the MeshingApplication, and to other applications, if needed along the process.

    Input:
    - model_coarse: Kratos model class before refinement
    - parameters_coarse: Kratos parameters class before refinement
    - minimal_size_value: minimal size after remeshing
    - maximal_size_value: maximal size after remeshing
    - metric_param: Kratos parameters class containing metric custom settings
    - remesh_param: Kratos parameters class containing remeshing custom settings
    - metric_name: string defining the metring which the class will use to build the metric
    """
    def __init__(self,current_level,model_coarse,parameters_coarse,metric_param,remesh_param,metric_name="hessian"):
        self.model_coarse = model_coarse
        self.parameters_coarse = parameters_coarse
        self.metric_param = metric_param
        self.remesh_param = remesh_param
        self.problem_type = self.parameters_coarse["solver_settings"]["solver_type"].GetString()
        self.metric = metric_name
        self.current_level = current_level
        self.wrapper = ParametersWrapper(self.parameters_coarse)

    def ComputeAdaptiveRefinement(self):
        """
        Method computing the refinement of the model based on the solution on the coarse mesh,
        exploiting the hessian metric of the solution.

        Input:
        - self: an instance of the class

        Output:
        - current_model_refined : Kratos model class after refinement
        - current_parameters_refined : Kratos parameters class after refinement
        """
        parameters_coarse = self.parameters_coarse
        model_coarse = self.model_coarse
        metric_param = self.metric_param
        remesh_param = self.remesh_param
        problem_type = self.problem_type
        current_level = self.current_level

        if (self.metric is "hessian"):
            # initialize interpolation error
            original_interp_error = metric_param["hessian_strategy_parameters"]["interpolation_error"].GetDouble()

            # problem dependent section
            # we build the metric depending on the problem
            if (problem_type == "potential_flow"):
                model_part_name = parameters_coarse["solver_settings"]["model_part_name"].GetString()
                # set NODAL_AREA and NODAL_H as non historical variables
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
                # Setting Metric Tensor to 0
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,model_coarse.GetModelPart(model_part_name).Nodes)
                # calculate NODAL_H
                find_nodal_h = KratosMultiphysics.FindNodalHProcess(model_coarse.GetModelPart(model_part_name))
                find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
                find_nodal_h.Execute()
                custom_gradient = KratosMultiphysics.CompressiblePotentialFlowApplication.ComputeCustomNodalGradientProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.VELOCITY, KratosMultiphysics.NODAL_AREA)
                custom_gradient.Execute()

                if metric_param.Has("local_gradient_variable"):
                    metric_param.RemoveValue("local_gradient_variable")
                if current_level > 0:
                    coefficient_interp_error =  metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
                    metric_param["hessian_strategy_parameters"].RemoveValue("coefficient_interpolation_error")
                    # interp_error = original_interp_error*(coefficient_interp_error)**(-current_level)
                    interp_error = original_interp_error/(coefficient_interp_error*current_level)
                    metric_param["hessian_strategy_parameters"]["interpolation_error"].SetDouble(interp_error)

                local_gradient = KratosMultiphysics.MeshingApplication.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.VELOCITY_X,metric_param)
                local_gradient.Execute()
                local_gradient = KratosMultiphysics.MeshingApplication.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.VELOCITY_Y,metric_param)
                local_gradient.Execute()

            elif (problem_type == "monolithic"):
                if have_fluid_dynamics_application:
                    if metric_param.Has("local_gradient_variable"):
                        metric_param.RemoveValue("local_gradient_variable")
                    if current_level > 0:
                        coefficient_interp_error =  metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
                        metric_param["hessian_strategy_parameters"].RemoveValue("coefficient_interpolation_error")
                        # interp_error = original_interp_error*(coefficient_interp_error)**(-current_level)
                        interp_error = original_interp_error/(coefficient_interp_error*current_level)
                        metric_param["hessian_strategy_parameters"]["interpolation_error"].SetDouble(interp_error)
                    model_part_name = parameters_coarse["solver_settings"]["model_part_name"].GetString()

                    # Setting Metric Tensor to 0
                    KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,model_coarse.GetModelPart(model_part_name).Nodes)

                    # calculate NODAL_H
                    find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
                    find_nodal_h.Execute()
                    local_gradient = KratosMultiphysics.MeshingApplication.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.FluidDynamicsApplication.AVERAGE_VELOCITY_X,metric_param)
                    local_gradient.Execute()
                    local_gradient = KratosMultiphysics.MeshingApplication.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.FluidDynamicsApplication.AVERAGE_VELOCITY_Y,metric_param)
                    local_gradient.Execute()
                else:
                    raise Exception("[MultilevelMonteCarloApplication]: Tried to perform adaptive refinement with problem type \"monolithic\", but FluidDynamicsApplication import failed. This application is required to for problem type \"monolithic\".")

            elif (problem_type == "stationary"):
                if current_level > 0:
                    coefficient_interp_error =  metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
                    metric_param["hessian_strategy_parameters"].RemoveValue("coefficient_interpolation_error")
                    interp_error = original_interp_error*(coefficient_interp_error)**(-current_level)
                    metric_param["hessian_strategy_parameters"]["interpolation_error"].SetDouble(interp_error)
                model_part_name = parameters_coarse["solver_settings"]["model_part_name"].GetString()
                # Setting Metric Tensor to 0
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,model_coarse.GetModelPart(model_part_name).Nodes)
                # calculate NODAL_H
                find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
                find_nodal_h.Execute()
                local_gradient = KratosMultiphysics.MeshingApplication.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),KratosMultiphysics.TEMPERATURE,metric_param)
                local_gradient.Execute()

            # create the remeshing process and execute it
            domain_size = self.wrapper.GetDomainSize()
            if domain_size == 2:
                MmgProcess = KratosMultiphysics.MeshingApplication.MmgProcess2D(model_coarse.GetModelPart(model_part_name),remesh_param)
            elif domain_size == 3:
                MmgProcess = KratosMultiphysics.MeshingApplication.MmgProcess3D(model_coarse.GetModelPart(model_part_name),remesh_param)
            else:
                err_msg = "Domain size is {}. Supported values are 2 and 3.\n".format(domain_size)
            raise Exception(err_msg)
            MmgProcess.Execute()

            # reset variables if needed
            model_coarse.GetModelPart(model_part_name).ProcessInfo.SetValue(KratosMultiphysics.TIME , 0.0)
            model_coarse.GetModelPart(model_part_name).ProcessInfo.SetValue(KratosMultiphysics.STEP , 0)

            # the refinement process empties the coarse model part object and fill it with the refined model part
            # the solution on the refined grid is obtained from the interpolation of the coarse solution
            # there are not other operations, therefore to build the new model we just need to take the updated coarse model
            current_model_refined = model_coarse
            current_parameters_refined = parameters_coarse
            return current_model_refined,current_parameters_refined

        else:
            err_msg = "Metric passed to the AdaptiveRefinement class is {}, but the only supported metric is \"hessian\".\n".format(self.metric)
            raise Exception(err_msg)

    def ComputeMeshSizeCoarsestLevel(self):
        """
        Method computing the mesh size of coarsest level, estimated as minimum nodal_h.

        Input:
        - self: an instance of the class
        """
        model_coarse = self.model_coarse
        model_part_name = self.wrapper.GetModelPartName()
        # set NODAL_AREA and NODAL_H as non historical variables
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, model_coarse.GetModelPart(model_part_name).Nodes)
        # calculate NODAL_H
        find_nodal_h = KratosMultiphysics.FindNodalHProcess(model_coarse.GetModelPart(model_part_name))
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
        find_nodal_h.Execute()
        # compute average mesh size
        mesh_size = 10.0
        for node in model_coarse.GetModelPart(model_part_name).Nodes:
            if (node.GetValue(KratosMultiphysics.NODAL_H) < mesh_size):
                mesh_size = node.GetValue(KratosMultiphysics.NODAL_H)
        self.mesh_size_coarsest_level = mesh_size

    def EstimateMeshSizeCurrentLevel(self):
        """
        Method estimating the mesh size of current level.

        Input:
        - self: an instance of the class
        """
        self.ComputeMeshSizeCoarsestLevel()
        current_level = self.current_level
        if (self.metric is "hessian"):
            original_interp_error = self.metric_param["hessian_strategy_parameters"]["interpolation_error"].GetDouble()
            domain_size = self.wrapper.GetDomainSize()
            if (domain_size == 2):
                coefficient = 2/9 # 2d
            elif (domain_size == 3):
                coefficient = 9/32 # 3d
            # TODO: compute below interp error level more automatically
            coefficient_interp_error =  self.metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
            # interp_error_level = original_interp_error*(coefficient_interp_error)**(-current_level)
            if (current_level > 0):
                interp_error_level = original_interp_error/(coefficient_interp_error*current_level)
            else: # current_level == 0
                interp_error_level = original_interp_error
            mesh_size_level = self.mesh_size_coarsest_level*np.sqrt(interp_error_level/original_interp_error) # relation from [Alauzet] eqs. pag 34 and 35
            self.mesh_size = mesh_size_level
