# Import Python libraries
from math import sqrt

# Import Kratos
import KratosMultiphysics
from KratosMultiphysics.MultilevelMonteCarloApplication.tools import ParametersWrapper
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if KratosMultiphysics.IsDistributedRun():
    import KratosMultiphysics.mpi
if CheckIfApplicationsAvailable("MeshingApplication"):
    import KratosMultiphysics.MeshingApplication

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

        # check MeshingApplication is imported,
        # otherwise raise an error
        if not CheckIfApplicationsAvailable("MeshingApplication"):
            raise Exception("[MultilevelMonteCarloApplication]: MeshingApplication cannot be imported, but it is necessary to perform adaptive refinement.")

        if (self.metric is "hessian"):
            # initialize interpolation error
            original_interp_error = metric_param["hessian_strategy_parameters"]["interpolation_error"].GetDouble()
            # set interpolation error for current level
            if current_level > 0:
                coefficient_interp_error =  metric_param["hessian_strategy_parameters"]["coefficient_interpolation_error"].GetDouble()
                metric_param["hessian_strategy_parameters"].RemoveValue("coefficient_interpolation_error")
                interp_error = original_interp_error*(coefficient_interp_error)**(-current_level)
                # interp_error = original_interp_error/(coefficient_interp_error*current_level)
                metric_param["hessian_strategy_parameters"]["interpolation_error"].SetDouble(interp_error)
            # Setting metric tensor to 0
            domain_size = self.wrapper.GetDomainSize()
            model_part_name = parameters_coarse["solver_settings"]["model_part_name"].GetString()
            if domain_size == 2:
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,model_coarse.GetModelPart(model_part_name).Nodes)
            elif domain_size == 3:
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D,model_coarse.GetModelPart(model_part_name).Nodes)
            else:
                err_msg = "Domain size is {}. Supported values are 2 and 3.\n".format(domain_size)
            # calculate NODAL_H
            find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(model_coarse.GetModelPart(model_part_name))
            find_nodal_h.Execute()

            # build the metric
            if metric_param["hessian_strategy_parameters"].Has("metric_variable"):
                metric_variables = self.__generate_variable_list_from_input(metric_param["hessian_strategy_parameters"]["metric_variable"])
                # remove metric value from settings, since we pass it directly in the constructor
                metric_param["hessian_strategy_parameters"].RemoveValue("metric_variable")
            else:
                raise Exception("A list of variable is expected under the key [\"hessian_strategy_parameters\"][\"metric_variable\"] of the \"metric\" dictionary.")

            for mv in metric_variables:
                local_gradient = KratosMultiphysics.MeshingApplication.ComputeHessianSolMetricProcess(model_coarse.GetModelPart(model_part_name),mv,metric_param)
                local_gradient.Execute()

            # create the remeshing process and execute it
            if KratosMultiphysics.IsDistributedRun(): # MPI
                KratosMultiphysics.mpi.ParallelFillCommunicator(model_coarse.GetModelPart(model_part_name)).Execute()
                if domain_size == 3:
                    if(hasattr(KratosMultiphysics.MeshingApplication,"ParMmgProcess3D")):
                        pmmg_process = KratosMultiphysics.MeshingApplication.ParMmgProcess3D(model_coarse.GetModelPart(model_part_name),remesh_param)
                    else:
                        raise Exception("[MultilevelMonteCarloApplication]: ParMmgProcess3D atrribute not found within MeshingApplcation. It is required to perform remeshing.")
                else:
                    err_msg = "Domain size is {}. Supported value is 3.\n".format(domain_size)
                    raise Exception(err_msg)
                pmmg_process.Execute()
            else: # serial
                if domain_size == 2:
                    if(hasattr(KratosMultiphysics.MeshingApplication,"MmgProcess2D")):
                        mmg_process = KratosMultiphysics.MeshingApplication.MmgProcess2D(model_coarse.GetModelPart(model_part_name),remesh_param)
                    else:
                        raise Exception("[MultilevelMonteCarloApplication]: MmgProcess2D atrribute not found within MeshingApplcation. It is required to perform remeshing.")
                elif domain_size == 3:
                    if(hasattr(KratosMultiphysics.MeshingApplication,"MmgProcess3D")):
                        mmg_process = KratosMultiphysics.MeshingApplication.MmgProcess3D(model_coarse.GetModelPart(model_part_name),remesh_param)
                    else:
                        raise Exception("[MultilevelMonteCarloApplication]: MmgProcess3D atrribute not found within MeshingApplcation. It is required to perform remeshing.")
                else:
                    err_msg = "Domain size is {}. Supported values are 2 and 3.\n".format(domain_size)
                    raise Exception(err_msg)
                mmg_process.Execute()

            # reset variables if needed
            model_coarse.GetModelPart(model_part_name).ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
            model_coarse.GetModelPart(model_part_name).ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)
            model_coarse.GetModelPart(model_part_name).ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED, False)

            if (problem_type in ["monolithic", "FractionalStep", "potential_flow"]):
                model_coarse.GetModelPart(model_part_name).RemoveSubModelPart("fluid_computational_model_part")

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
            mesh_size_level = self.mesh_size_coarsest_level*sqrt(interp_error_level/original_interp_error) # relation from [Alauzet] eqs. pag 34 and 35
            self.mesh_size = mesh_size_level

    def __generate_variable_list_from_input(self, param):
        """
        Method taken and adapted from mmg_progess.py of MeshingApplication.
        Parse a list of variables from input.
        """
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        variable_list = []
        param_names = param.GetStringArray()
        for variable_name in param_names:
            varriable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
            if varriable_type == "Double" or varriable_type == "Component":
                variable_list.append(KratosMultiphysics.KratosGlobals.GetVariable(variable_name))
            else:
                variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( variable_name + "_X" ))
                variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( variable_name + "_Y" ))
                if self.wrapper.GetDomainSize() == 3:
                    variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( variable_name + "_Z" ))

        return variable_list
