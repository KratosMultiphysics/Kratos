from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

try:
    import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
    structural_dependencies = True
    missing_application = ''
except ImportError as e:
    structural_dependencies = False
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''',
                                    '{0}'.format(e)).group(1)

from json_utilities import *
import json
import os

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MmgProcess(Model, settings["Parameters"])

class MmgProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "help"                             : "This process remeshes using MMG library. This process uses different utilities and processes",
            "mesh_id"                          : 0,
            "filename"                         : "out",
            "model_part_name"                  : "MainModelPart",
            "strategy"                         : "LevelSet",
            "level_set_strategy_parameters"              :{
                "scalar_variable"                  : "DISTANCE",
                "gradient_variable"                : "DISTANCE_GRADIENT"
            },
            "error_strategy_parameters"              :{
                "compute_error_extra_parameters":
                {
                    "stress_vector_variable"              : "CAUCHY_STRESS_VECTOR"
                },
                "error_metric_parameters"                 :
                {
                    "error_threshold"                       : 0.05,
                    "interpolation_error"                   : 0.04
                },
                "set_target_number_of_elements"       : false,
                "target_number_of_elements"           : 1000,
                "perform_nodal_h_averaging"           : false
                "max_iterations"                      : 3
            },
            "framework"                            : "Eulerian",
            "internal_variables_parameters"        :
            {
                "allocation_size"                      : 1000,
                "bucket_size"                          : 4,
                "search_factor"                        : 2,
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" :[]
            },
            "hessian_strategy_parameters"              :{
                "metric_variable"                  : ["DISTANCE"],
                "estimate_interpolation_error"     : false,
                "interpolation_error"              : 0.04,
                "mesh_dependent_constant"          : 0.28125
            },
            "enforce_current"                  : true,
            "initial_step"                     : 1,
            "step_frequency"                   : 0,
            "automatic_remesh"                 : true,
            "automatic_remesh_parameters"      :{
                "automatic_remesh_type"            : "Ratio",
                "min_size_ratio"                   : 1.0,
                "max_size_ratio"                   : 3.0,
                "refer_type"                       : "Mean",
                "min_size_current_percentage"      : 50.0,
                "max_size_current_percentage"      : 98.0
            },
            "initial_remeshing"                : false,
            "fix_contour_model_parts"          : [],
            "force_min"                        : false,
            "minimal_size"                     : 0.1,
            "force_max"                        : false,
            "maximal_size"                     : 10.0,
            "advanced_parameters"                  :
            {
                "force_hausdorff_value"               : false,
                "hausdorff_value"                     : 0.0001,
                "no_move_mesh"                        : false,
                "no_surf_mesh"                        : false,
                "no_insert_mesh"                      : false,
                "no_swap_mesh"                        : false,
                "deactivate_detect_angle"             : false,
                "force_gradation_value"               : false,
                "gradation_value"                     : 1.3
            },
            "anisotropy_remeshing"             : true,
            "anisotropy_parameters":{
                "reference_variable_name"          : "DISTANCE",
                "hmin_over_hmax_anisotropic_ratio" : 0.01,
                "boundary_layer_max_distance"      : 1.0,
                "boundary_layer_min_size_ratio"    : 2.0,
                "interpolation"                    : "Linear"
            },
            "save_external_files"              : false,
            "max_number_of_searchs"            : 1000,
            "debug_mode"                       : false,
            "echo_level"                       : 3
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.model_part= Model[self.settings["model_part_name"].GetString()]
        self.dim = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        self.enforce_current = self.settings["enforce_current"].GetBool()

        self.initial_remeshing = self.settings["initial_remeshing"].GetBool()

        self.initial_step = self.settings["initial_step"].GetInt()
        self.step_frequency = self.settings["step_frequency"].GetInt()

    def ExecuteInitialize(self):
        # Calculate NODAL_H
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.model_part)
        self.find_nodal_h.Execute()

        # Calculate the parameters of automatic remeshing
        if (self.settings["automatic_remesh"].GetBool() == True):
            import statistics as stat
            nodal_h_values = []
            for node in self.model_part.Nodes:
                nodal_h_values.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_H))

            # Calculate the minimum size
            if (self.settings["automatic_remesh_parameters"]["automatic_remesh_type"].GetString() == "Ratio"):
                # NOTE: For mode: https://docs.python.org/3/library/statistics.html
                if (self.settings["automatic_remesh_parameters"]["refer_type"].GetString() == "Mean"):
                    ref = stat.mean(nodal_h_values)
                elif (self.settings["automatic_remesh_parameters"]["refer_type"].GetString() == "Median"):
                    ref = stat.median(nodal_h_values)

                self.settings["minimal_size"].SetDouble(ref * (self.settings["automatic_remesh_parameters"]["min_size_ratio"].GetDouble()))
                self.settings["maximal_size"].SetDouble(ref * (self.settings["automatic_remesh_parameters"]["max_size_ratio"].GetDouble()))
            elif (self.settings["automatic_remesh_parameters"]["automatic_remesh_type"].GetString() == "Percentage"):
                mean = stat.mean(nodal_h_values)
                stdev = stat.stdev(nodal_h_values)
                prob = (self.settings["automatic_remesh_parameters"]["min_size_current_percentage"].GetDouble())/100
                self.settings["minimal_size"].SetDouble(_normvalf(prob, mean, stdev)) # Using normal normal distribution to get the minimal size as a stadistical meaninful value

                prob = (self.settings["automatic_remesh_parameters"]["max_size_current_percentage"].GetDouble())/100
                self.settings["maximal_size"].SetDouble(_normvalf(prob, mean, stdev)) # Using normal normal distribution to get the maximal size as a stadistical meaninful value

        # Anisotropic remeshing parameters
        self.anisotropy_remeshing = self.settings["anisotropy_remeshing"].GetBool()
        if (self.anisotropy_remeshing == True):
            if (self.settings["automatic_remesh"].GetBool() == True):
                self.settings["anisotropy_parameters"]["boundary_layer_max_distance"].SetDouble(self.settings["minimal_size"].GetDouble() * self.settings["anisotropy_parameters"]["boundary_layer_min_size_ratio"].GetDouble())

        # Select the remeshing strategy
        self.strategy = self.settings["strategy"].GetString()
        if (self.strategy == "LevelSet"):
            self.scalar_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.settings["level_set_strategy_parameters"]["scalar_variable"].GetString() )
            self.gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.settings["level_set_strategy_parameters"]["gradient_variable"].GetString() )
        elif (self.strategy == "Hessian"):
            self.metric_variable = self.__generate_variable_list_from_input(self.settings["hessian_strategy_parameters"]["metric_variable"])
            mesh_dependent_constant = self.settings["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble()
            if (mesh_dependent_constant == 0.0):
                self.settings["hessian_strategy_parameters"]["mesh_dependent_constant"].SetDouble(0.5 * (self.dim/(self.dim + 1))**2.0)
        elif (self.strategy == "superconvergent_patch_recovery"):
            self.error_threshold = self.settings["error_strategy_parameters"]["error_metric_parameters"]["error_threshold"].GetDouble()
            self.estimated_error = 0
            self.remeshing_cycle = 0
            self.model_part.ProcessInfo[MeshingApplication.EXECUTE_REMESHING] = True

        self.internal_variable_interpolation_list = self.__generate_internal_variable_list_from_input(self.settings["internal_variables_parameters"]["internal_variable_interpolation_list"])

        # NOTE: Add more model part if interested
        submodelpartslist = self.__generate_submodelparts_list_from_input(self.settings["fix_contour_model_parts"])

        for submodelpart in submodelpartslist:
            for node in submodelpart.Nodes:
                node.Set(KratosMultiphysics.BLOCKED, True)

        if (self.strategy == "LevelSet"):
            self._CreateGradientProcess()

        if (self.dim == 2):
            self.initialize_metric = MeshingApplication.MetricFastInit2D(self.model_part)
        else:
            self.initialize_metric = MeshingApplication.MetricFastInit3D(self.model_part)

        self.initialize_metric.Execute()

        self._CreateMetricsProcess()

        mmg_parameters = KratosMultiphysics.Parameters("""{"force_sizes":{}}""")
        mmg_parameters.AddValue("filename",self.settings["filename"])
        mmg_parameters.AddValue("framework",self.settings["framework"])
        mmg_parameters.AddValue("internal_variables_parameters",self.settings["internal_variables_parameters"])
        mmg_parameters.AddValue("save_external_files",self.settings["save_external_files"])
        mmg_parameters.AddValue("max_number_of_searchs",self.settings["max_number_of_searchs"])
        mmg_parameters["force_sizes"].AddValue("force_min",self.settings["force_min"])
        mmg_parameters["force_sizes"].AddValue("minimal_size",self.settings["maximal_size"])
        mmg_parameters["force_sizes"].AddValue("force_max",self.settings["force_max"])
        mmg_parameters["force_sizes"].AddValue("maximal_size",self.settings["maximal_size"])
        mmg_parameters.AddValue("advanced_parameters",self.settings["advanced_parameters"])
        mmg_parameters.AddValue("echo_level",self.settings["echo_level"])
        if (self.dim == 2):
            self.mmg_process = MeshingApplication.MmgProcess2D(self.model_part, mmg_parameters)
        else:
            self.mmg_process = MeshingApplication.MmgProcess3D(self.model_part, mmg_parameters)

        # We reset the step
        self.step = 0

        # We compute initial remeshing is desired
        if (self.initial_remeshing == True):
            if (self.model_part.Is(KratosMultiphysics.MODIFIED) == False):
                self._ExecuteRefinement()
            else:
                self.model_part.Set(KratosMultiphysics.MODIFIED, False)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        if (self.initial_remeshing == False):
            # We need to check if the model part has been modified recently
            if (self.model_part.Is(KratosMultiphysics.MODIFIED) == True):
                self.model_part.Set(KratosMultiphysics.MODIFIED, False)
                self.step = 0  # Reset (just to be sure)
            else:
                self.step += 1
                if self.step_frequency > 0:
                    if self.step >= self.step_frequency:
                        if self.model_part.ProcessInfo[KratosMultiphysics.STEP] >= self.initial_step:
                            self._ExecuteRefinement()
                            self.step = 0  # Reset

    def ExecuteFinalizeSolutionStep(self):
        if (self.strategy == "superconvergent_patch_recovery"):
            self._ErrorCalculation()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        if (self.strategy == "superconvergent_patch_recovery"):
            if (self.model_part.ProcessInfo[MeshingApplication.ERROR_ESTIMATE] > self.error_threshold):
                self.__execute_refinement()
            self.remeshing_cycle += 1
            if (self.model_part.ProcessInfo[MeshingApplication.ERROR_ESTIMATE] <= self.error_threshold or self.remeshing_cycle > self.params["max_iterations"].GetInt()):
                self.model_part.ProcessInfo[MeshingApplication.EXECUTE_REMESHING] = False

    def ExecuteFinalize(self):
        pass

    def _CreateMetricsProcess(self):
        self.metric_processes = []
        if (self.strategy == "LevelSet"):
            level_set_parameters = KratosMultiphysics.Parameters("""{}""")
            level_set_parameters.AddValue("minimal_size",self.settings["minimal_size"])
            level_set_parameters.AddValue("enforce_current",self.settings["enforce_current"])
            level_set_parameters.AddValue("anisotropy_remeshing",self.settings["anisotropy_remeshing"])
            level_set_parameters.AddValue("anisotropy_parameters",self.settings["anisotropy_parameters"])
            if (self.dim == 2):
                self.metric_processes.append(MeshingApplication.ComputeLevelSetSolMetricProcess2D(
                    self.model_part,
                    self.gradient_variable,
                    level_set_parameters))

            else:
                self.metric_processes.append(MeshingApplication.ComputeLevelSetSolMetricProcess3D(
                    self.model_part,
                    self.gradient_variable,
                    level_set_parameters))

        elif (self.strategy == "Hessian"):
            hessian_parameters = KratosMultiphysics.Parameters("""{}""")
            hessian_parameters.AddValue("minimal_size",self.settings["minimal_size"])
            hessian_parameters.AddValue("maximal_size",self.settings["maximal_size"])
            hessian_parameters.AddValue("enforce_current",self.settings["enforce_current"])
            hessian_parameters.AddValue("hessian_strategy_parameters",self.settings["hessian_strategy_parameters"])
            hessian_parameters.AddValue("anisotropy_remeshing",self.settings["anisotropy_remeshing"])
            hessian_parameters.AddValue("anisotropy_parameters",self.settings["anisotropy_parameters"])
            for current_metric_variable in self.metric_variable:
                if (type(current_metric_variable) is KratosMultiphysics.Array1DComponentVariable):
                    if (self.dim == 2):
                        self.metric_processes.append(MeshingApplication.ComputeHessianSolMetricProcessComp2D(
                            self.model_part,
                            current_metric_variable,
                            hessian_parameters))
                    else:
                        self.metric_processes.append(MeshingApplication.ComputeHessianSolMetricProcessComp3D(
                            self.model_part,
                            current_metric_variable,
                            hessian_parameters))
                else:
                    if (self.dim == 2):
                        self.metric_processes.append(MeshingApplication.ComputeHessianSolMetricProcess2D(
                            self.model_part,
                            current_metric_variable,
                            hessian_parameters))
                    else:
                        self.metric_processes.append(MeshingApplication.ComputeHessianSolMetricProcess3D(
                            self.model_part,
                            current_metric_variable,
                            hessian_parameters))
        elif (self.strategy == "superconvergent_patch_recovery"):
            if not structural_dependencies:
                raise Exception("You need to compile the StructuralMechanicsApplication in order to use this criteria")

            # We compute the error
            error_compute_parameters = KratosMultiphysics.Parameters("""{}""")
            error_compute_parameters.AddValue("stress_vector_variable", self.settings["compute_error_extra_parameters"]["stress_vector_variable"])
            error_compute_parameters.AddValue("echo_level", self.settings["echo_level"])
            if (self.dim == 2):
                self.error_compute = StructuralMechanicsApplication.SPRErrorProcess2D(
                    self.model_part,
                    error_compute_parameters
                    )
            else:
                self.error_compute = StructuralMechanicsApplication.SPRErrorProcess3D(
                    self.model_part,
                    error_compute_parameters
                    )

            # Now we compute the metric
            error_metric_parameters = KratosMultiphysics.Parameters("""{}""")
            error_metric_parameters.AddValue("minimal_size",self.settings["minimal_size"])
            error_metric_parameters.AddValue("maximal_size",self.settings["maximal_size"])
            error_metric_parameters.AddValue("target_error",self.settings["error_strategy_parameters"]["error_metric_parameters"]["interpolation_error"])
            error_metric_parameters.AddValue("set_target_number_of_elements", self.settings["error_strategy_parameters"]["set_target_number_of_elements"])
            error_metric_parameters.AddValue("target_number_of_elements", self.settings["error_strategy_parameters"]["target_number_of_elements"])
            error_metric_parameters.AddValue("perform_nodal_h_averaging", self.settings["error_strategy_parameters"]["perform_nodal_h_averaging"])
            error_metric_parameters.AddValue("echo_level", self.settings["echo_level"])

            if (self.dim == 2):
                self.metric_process = MeshingApplication.MetricErrorProcess2D(
                    self.model_part,
                    error_metric_parameters
                    )
            else:
                self.metric_process = MeshingApplication.MetricErrorProcess3D(
                    self.model_part,
                    error_metric_parameters
                    )

    def _CreateGradientProcess(self):
        # We compute the scalar value gradient
        if (self.dim == 2):
            self.local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.model_part, self.scalar_variable, self.gradient_variable, KratosMultiphysics.NODAL_AREA)
        else:
            self.local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.model_part, self.scalar_variable, self.gradient_variable, KratosMultiphysics.NODAL_AREA)

    def _ExecuteRefinement(self):
        if (self.strategy == "LevelSet"):
            # Calculate the gradient
            self.local_gradient.Execute()

        # Recalculate NODAL_H
        self.find_nodal_h.Execute()

        # Initialize metric
        self.initialize_metric.Execute()

        KratosMultiphysics.Logger.PrintInfo("MMG Remeshing Process", "Calculating the metrics")
        # Execute metric computation
        for metric_process in self.metric_processes:
            metric_process.Execute()

        KratosMultiphysics.Logger.PrintInfo("MMG Remeshing Process", "Remeshing")
        self.mmg_process.Execute()

        if (self.settings["debug_mode"].GetBool() == True):
            self.gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
            self.singlefile = KratosMultiphysics.MultiFileFlag.SingleFile
            self.deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
            self.write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions
            self._debug_output(self.step, "")

        if (self.strategy == "LevelSet"):
            self.local_gradient.Execute() # Recalculate gradient after remeshing

        # Recalculate NODAL_H
        self.find_nodal_h.Execute()

        # We need to set that the model part has been modified (later on we will act in consequence)
        self.model_part.Set(KratosMultiphysics.MODIFIED, True)

        KratosMultiphysics.Logger.PrintInfo("MMG Remeshing Process", "Remesh finished")

    def _ErrorCalculation(self):

        # Initialize metric
        self.initialize_metric.Execute()

        KratosMultiphysics.Logger.PrintInfo("MMG Remeshing Process", "Calculating the metrics")
        # Execute error computation
        self.error_compute.Execute()
        # Execute metric computation
        self.metric_process.Execute()
        self.estimated_error = self.model_part.ProcessInfo[MeshingApplication.ERROR_ESTIMATE]

    def __generate_submodelparts_list_from_input(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve submodelparts name from input (a string) and request the corresponding C++ object to the kernel
        return [self.model_part.GetSubModelPart(param[i].GetString()) for i in range(0, param.size())]

    def __generate_variable_list_from_input(self,param):
      '''Parse a list of variables from input.'''
      # At least verify that the input is a string
      if not param.IsArray():
          raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

      # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel

      variable_list = []
      if (len(self.model_part.Nodes) > 0):
          node = (self.model_part.Nodes)[1]
          for i in range( 0,param.size()):
              aux_var = KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() )
              val = node.GetSolutionStepValue(aux_var, 0)
              if isinstance(val,float):
                  variable_list.append(aux_var)
              else:
                  variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString()+"_X" ))
                  variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString()+"_Y" ))
                  if (self.dim == 3):
                      variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString()+"_Z" ))

      return variable_list

    def __generate_internal_variable_list_from_input(self,param):
      '''Parse a list of variables from input.'''
      # At least verify that the input is a string
      if not param.IsArray():
          raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

      # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel

      variable_list = []

      for i in range( 0,param.size()):
          aux_var = KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() )
          variable_list.append(aux_var)

      return variable_list

    def _debug_output(self, label, name):

        gid_io = KratosMultiphysics.GidIO("REMESHING_"+name+"_STEP_"+str(label), self.gid_mode, self.singlefile, self.deformed_mesh_flag, self.write_conditions)

        gid_io.InitializeMesh(label)
        gid_io.WriteMesh(self.model_part.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(label, self.model_part.GetMesh())
        if (self.settings["framework"].GetString() ==  "Lagrangian"):
            gid_io.WriteNodalResults(KratosMultiphysics.DISPLACEMENT, self.model_part.Nodes, label, 0)
            for var in self.internal_variable_interpolation_list:
                gid_io.PrintOnGaussPoints(var, self.model_part, label)
        else:
            gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY, self.model_part.Nodes, label, 0)
        gid_io.FinalizeResults()

        #raise NameError("DEBUG")

def _linear_interpolation(x, x_list, y_list):
    tb = KratosMultiphysics.PiecewiseLinearTable()
    for i in range(len(x_list)):
        tb.AddRow(x_list[i], y_list[i])

    return tb.GetNearestValue(x)

def _normpdf(x, mean, sd):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data = read_external_json(dir_path+"/normal_distribution.json")
    z = (x-mean)/sd
    z_list = data["Z"]
    prob_list = data["Prob"]
    if (z > 0):
        prob = _linear_interpolation(z, z_list, prob_list)
    else:
        prob = 1.0 - _linear_interpolation(-z, z_list, prob_list)
    return prob


def _normvalf(prob, mean, sd):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data = read_external_json(dir_path+"/normal_distribution.json")
    z_list = data["Z"]
    prob_list = data["Prob"]
    if (prob >= 0.5):
        z = _linear_interpolation(prob, prob_list, z_list)
    else:
        z = - _linear_interpolation(1.0 - prob, prob_list, z_list)
    x = z * sd + mean
    return x
