from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication

from json_utilities import *
import json
import os

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MmgProcess(Model, settings["Parameters"])

class MmgProcess(KratosMultiphysics.Process):

    def __init__(self,Model,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                          : 0,
            "filename"                         : "out",
            "model_part_name"                  : "MainModelPart",
            "strategy"                         : "LevelSet",
            "level_set_strategy_parameters"              :{
                "scalar_variable"                  : "DISTANCE",
                "gradient_variable"                : "DISTANCE_GRADIENT"
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
                "interpolation_error"              : 0.04,
                "mesh_dependent_constant"          : 0.0
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
            "initial_remeshing"                : true,
            "fix_contour_model_parts"          : [],
            "minimal_size"                     : 0.1,
            "maximal_size"                     : 10.0,
            "anisotropy_remeshing"             : true,
            "anisotropy_parameters":{
                "hmin_over_hmax_anisotropic_ratio" : 0.01,
                "boundary_layer_max_distance"      : 1.0,
                "boundary_layer_min_size_ratio"    : 2.0,
                "interpolation"                    : "Linear"
            },
            "save_external_files"              : false,
            "max_number_of_searchs"            : 1000,
            "debug_mode"                  : false,
            "echo_level"                       : 3
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.Model= Model
        self.model_part_name = self.params["model_part_name"].GetString()
        self.dim = self.Model[self.model_part_name].ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.params = params

        self.enforce_current = self.params["enforce_current"].GetBool()

        self.initial_remeshing = self.params["initial_remeshing"].GetBool()

        # Select the remeshing strategy
        self.strategy = self.params["strategy"].GetString()
        if (self.strategy == "LevelSet"):
            self.scalar_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["level_set_strategy_parameters"]["scalar_variable"].GetString() )
            self.gradient_variable = KratosMultiphysics.KratosGlobals.GetVariable( self.params["level_set_strategy_parameters"]["gradient_variable"].GetString() )
        elif (self.strategy == "Hessian"):
            self.metric_variable = self.__generate_variable_list_from_input(self.params["hessian_strategy_parameters"]["metric_variable"])
            mesh_dependent_constant = self.params["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble()
            if (mesh_dependent_constant == 0.0):
                self.params["hessian_strategy_parameters"]["mesh_dependent_constant"].SetDouble(0.5 * (self.dim/(self.dim + 1))**2.0)
        
        # Calculate NODAL_H
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.Model[self.model_part_name])
        self.find_nodal_h.Execute()

        # Calculate the parameters of automatic remeshing
        if (self.params["automatic_remesh"].GetBool() == True):
            import statistics as stat
            nodal_h_values = []
            for node in self.Model[self.model_part_name].Nodes:
                nodal_h_values.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_H))

            # Calculate the minimum size
            if (self.params["automatic_remesh_parameters"]["automatic_remesh_type"].GetString() == "Ratio"):
                # NOTE: For mode: https://docs.python.org/3/library/statistics.html
                if (self.params["automatic_remesh_parameters"]["refer_type"].GetString() == "Mean"):
                    ref = stat.mean(nodal_h_values)
                elif (self.params["automatic_remesh_parameters"]["refer_type"].GetString() == "Median"):
                    ref = stat.median(nodal_h_values)

                self.params["minimal_size"].SetDouble(ref * (self.params["automatic_remesh_parameters"]["min_size_ratio"].GetDouble()))
                self.params["maximal_size"].SetDouble(ref * (self.params["automatic_remesh_parameters"]["max_size_ratio"].GetDouble()))
            elif (self.params["automatic_remesh_parameters"]["automatic_remesh_type"].GetString() == "Percentage"):
                mean = stat.mean(nodal_h_values)
                stdev = stat.stdev(nodal_h_values)
                prob = (self.params["automatic_remesh_parameters"]["min_size_current_percentage"].GetDouble())/100
                self.params["minimal_size"].SetDouble(_normvalf(prob, mean, stdev)) # Using normal normal distribution to get the minimal size as a stadistical meaninful value

                prob = (self.params["automatic_remesh_parameters"]["max_size_current_percentage"].GetDouble())/100
                self.params["maximal_size"].SetDouble(_normvalf(prob, mean, stdev)) # Using normal normal distribution to get the maximal size as a stadistical meaninful value

        # Anisotropic remeshing parameters
        self.anisotropy_remeshing = self.params["anisotropy_remeshing"].GetBool()
        if (self.anisotropy_remeshing == True):
            if (self.params["automatic_remesh"].GetBool() == True):
                self.params["anisotropy_parameters"]["boundary_layer_max_distance"].SetDouble(self.params["minimal_size"].GetDouble() * self.params["anisotropy_parameters"]["boundary_layer_min_size_ratio"].GetDouble())

        self.initial_step = self.params["initial_step"].GetInt()
        self.step_frequency = self.params["step_frequency"].GetInt()

    def ExecuteInitialize(self):

        # NOTE: Add more model part if interested
        submodelpartslist = self.__generate_submodelparts_list_from_input(self.params["fix_contour_model_parts"])

        for submodelpart in submodelpartslist:
            for node in submodelpart.Nodes:
                node.Set(KratosMultiphysics.BLOCKED, True)

        if (self.strategy == "LevelSet"):
            self._CreateGradientProcess()

        if (self.dim == 2):
            self.initialize_metric = MeshingApplication.MetricFastInit2D(self.Model[self.model_part_name])
        else:
            self.initialize_metric = MeshingApplication.MetricFastInit3D(self.Model[self.model_part_name])
            
        self.initialize_metric.Execute()

        self._CreateMetricsProcess()

        mmg_parameters = KratosMultiphysics.Parameters("""{}""")
        mmg_parameters.AddValue("filename",self.params["filename"])
        mmg_parameters.AddValue("framework",self.params["framework"])
        mmg_parameters.AddValue("internal_variables_parameters",self.params["internal_variables_parameters"])
        mmg_parameters.AddValue("save_external_files",self.params["save_external_files"])
        mmg_parameters.AddValue("max_number_of_searchs",self.params["max_number_of_searchs"])
        mmg_parameters.AddValue("echo_level",self.params["echo_level"])
        if (self.dim == 2):
            self.MmgProcess = MeshingApplication.MmgProcess2D(self.Model[self.model_part_name], mmg_parameters)
        else:
            self.MmgProcess = MeshingApplication.MmgProcess3D(self.Model[self.model_part_name], mmg_parameters)

        if (self.initial_remeshing == True):
            self._ExecuteRefinement()

    def ExecuteBeforeSolutionLoop(self):
        self.step = 0

    def ExecuteInitializeSolutionStep(self):
        # We need to check if the model part has been modified recently
        if (self.Model[self.model_part_name].Is(KratosMultiphysics.MODIFIED) == True):
            self.Model[self.model_part_name].Set(KratosMultiphysics.MODIFIED, False)
            self.step = 0  # Reset (just to be sure)
        else:
            self.step += 1
            if self.step_frequency > 0:
                if self.step >= self.step_frequency:
                    if self.Model[self.model_part_name].ProcessInfo[KratosMultiphysics.TIME_STEPS] >= self.initial_step:
                        self._ExecuteRefinement()
                        self.step = 0  # Reset

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def _CreateMetricsProcess(self):
        self.MetricsProcess = []
        if (self.strategy == "LevelSet"):
            level_set_parameters = KratosMultiphysics.Parameters("""{}""")
            level_set_parameters.AddValue("minimal_size",self.params["minimal_size"])
            level_set_parameters.AddValue("enforce_current",self.params["enforce_current"])
            level_set_parameters.AddValue("anisotropy_remeshing",self.params["anisotropy_remeshing"])
            level_set_parameters.AddValue("anisotropy_parameters",self.params["anisotropy_parameters"])
            if (self.dim == 2):
                self.MetricsProcess.append(MeshingApplication.ComputeLevelSetSolMetricProcess2D(
                    self.Model[self.model_part_name],
                    self.gradient_variable,
                    level_set_parameters))

            else:
                self.MetricsProcess.append(MeshingApplication.ComputeLevelSetSolMetricProcess3D(
                    self.Model[self.model_part_name],
                    self.gradient_variable,
                    level_set_parameters))

        elif (self.strategy == "Hessian"):
            hessian_parameters = KratosMultiphysics.Parameters("""{}""")
            hessian_parameters.AddValue("minimal_size",self.params["minimal_size"])
            hessian_parameters.AddValue("maximal_size",self.params["maximal_size"])
            hessian_parameters.AddValue("enforce_current",self.params["enforce_current"])
            hessian_parameters.AddValue("hessian_strategy_parameters",self.params["hessian_strategy_parameters"])
            hessian_parameters.AddValue("anisotropy_remeshing",self.params["anisotropy_remeshing"])
            hessian_parameters.AddValue("anisotropy_parameters",self.params["anisotropy_parameters"])
            for current_metric_variable in self.metric_variable:
                if (type(current_metric_variable) is KratosMultiphysics.Array1DComponentVariable):
                    if (self.dim == 2):
                        self.MetricsProcess.append(MeshingApplication.ComputeHessianSolMetricProcessComp2D(
                            self.Model[self.model_part_name],
                            current_metric_variable,
                            hessian_parameters))
                    else:
                        self.MetricsProcess.append(MeshingApplication.ComputeHessianSolMetricProcessComp3D(
                            self.Model[self.model_part_name],
                            current_metric_variable,
                            hessian_parameters))
                else:
                    if (self.dim == 2):
                        self.MetricsProcess.append(MeshingApplication.ComputeHessianSolMetricProcess2D(
                            self.Model[self.model_part_name],
                            current_metric_variable,
                            hessian_parameters))
                    else:
                        self.MetricsProcess.append(MeshingApplication.ComputeHessianSolMetricProcess3D(
                            self.Model[self.model_part_name],
                            current_metric_variable,
                            hessian_parameters))

    def _CreateGradientProcess(self):
        # We compute the scalar value gradient
        if (self.dim == 2):
            self.local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.Model[self.model_part_name], self.scalar_variable, self.gradient_variable, KratosMultiphysics.NODAL_AREA)
        else:
            self.local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.Model[self.model_part_name], self.scalar_variable, self.gradient_variable, KratosMultiphysics.NODAL_AREA)

    def _ExecuteRefinement(self):
        if (self.strategy == "LevelSet"):
            # Calculate the gradient
            self.local_gradient.Execute()

        # Recalculate NODAL_H
        self.find_nodal_h.Execute()

        # Initialize metric
        self.initialize_metric.Execute()

        print("Calculating the metrics")
        # Execute metric computation
        for metric_process in self.MetricsProcess:
            metric_process.Execute()

        print("Remeshing")
        self.MmgProcess.Execute()

        if (self.params["debug_mode"].GetBool() == True):
            self.gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
            self.singlefile = KratosMultiphysics.MultiFileFlag.SingleFile
            self.deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
            self.write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions
            self._debug_output(self.step, "")

        if (self.strategy == "LevelSet"):
            self.local_gradient.Execute() # Recalculate gradient after remeshing

        # We need to set that the model part has been modified (later on we will act in consequence)
        self.Model[self.model_part_name].Set(KratosMultiphysics.MODIFIED, True)

        print("Remesh finished")

    def __generate_submodelparts_list_from_input(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve submodelparts name from input (a string) and request the corresponding C++ object to the kernel
        return [self.Model[self.model_part_name].GetSubModelPart(param[i].GetString()) for i in range(0, param.size())]

    def __generate_variable_list_from_input(self,param):
      '''Parse a list of variables from input.'''
      # At least verify that the input is a string
      if not param.IsArray():
          raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

      # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel

      variable_list = []

      for i in range( 0,param.size()):
          aux_var = KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() )
          for node in self.Model[self.model_part_name].Nodes:
            val = node.GetSolutionStepValue(aux_var, 0)
            break
          if isinstance(val,float):
              variable_list.append(aux_var)
          else:
              variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString()+"_X" ))
              variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString()+"_Y" ))
              if (self.dim == 3):
                variable_list.append( KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString()+"_Z" ))

      return variable_list

    def _debug_output(self, label, name):

        gid_io = KratosMultiphysics.GidIO("REMESHING_"+name+"_STEP_"+str(label), self.gid_mode, self.singlefile, self.deformed_mesh_flag, self.write_conditions)
        
        gid_io.InitializeMesh(label)
        gid_io.WriteMesh(self.Model[self.model_part_name].GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(label, self.Model[self.model_part_name].GetMesh())
        if (self.params["framework"].GetString() ==  "Lagrangian"):
            gid_io.WriteNodalResults(KratosMultiphysics.DISPLACEMENT, self.Model[self.model_part_name].Nodes, label, 0)
        else:
            gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY, self.Model[self.model_part_name].Nodes, label, 0)
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
