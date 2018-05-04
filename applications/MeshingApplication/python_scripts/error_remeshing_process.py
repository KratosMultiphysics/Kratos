from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.MeshingApplication as MeshingApplication
KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ErrorRemeshingProcess(Model, settings["Parameters"])

class ErrorRemeshingProcess(KratosMultiphysics.Process):

    def __init__(self,Model,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                          : 0,
            "filename"                         : "out",
            "model_part_name"                  : "MainModelPart",
            "framework"                        : "Lagrangian",
            "internal_variables_parameters"    :
            {
                "allocation_size"                      : 1000, 
                "bucket_size"                          : 4, 
                "search_factor"                        : 2, 
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" :[]
            },
            "error_parameters"                 :
            {
                "error_threshold"                       : 0.05,
                "interpolation_error"                   : 0.04
            },
            "average_nodal_h"                  : false,
            "fix_contour_model_parts"          : [],
            "minimal_size"                     : 0.01,
            "maximal_size"                     : 10.0,
            "save_external_files"              : false,
            "save_mdpa_file"                   : false,
            "max_number_of_searchs"            : 1000,
            "debug_mode"                       : false,
            "echo_level"                       : 3,
            "set_number_of_elements"              : false,
            "number_of_elements"                  : 1000,
            "max_iterations"                      : 3
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.Model= Model
        self.model_part_name = self.params["model_part_name"].GetString()
        self.main_model_part = self.Model[self.model_part_name]
        self.dim = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.params = params
        
        self.error_threshold = self.params["error_parameters"]["error_threshold"].GetDouble()
        self.estimated_error = 0
        self.remeshing_cycle = 0
        self.main_model_part.ProcessInfo[MeshingApplication.EXECUTE_REMESHING] = True

    def ExecuteInitialize(self):
        # NOTE: Add more model part if interested
        submodelpartslist = self.__generate_submodelparts_list_from_input(self.params["fix_contour_model_parts"])

        for submodelpart in submodelpartslist:
            for node in submodelpart.Nodes:
                node.Set(KratosMultiphysics.BLOCKED, True)
            del(node)
        
        if (self.dim == 2):
            self.initialize_metric = MeshingApplication.MetricFastInit2D(self.main_model_part)
        else:
            self.initialize_metric = MeshingApplication.MetricFastInit3D(self.main_model_part)
            
        self.initialize_metric.Execute()

        self.__create_metric_process()

        mmg_parameters = KratosMultiphysics.Parameters("""{}""")
        mmg_parameters.AddValue("filename",self.params["filename"])
        mmg_parameters.AddValue("framework",self.params["framework"])
        mmg_parameters.AddValue("internal_variables_parameters",self.params["internal_variables_parameters"])
        mmg_parameters.AddValue("save_external_files",self.params["save_external_files"])
        mmg_parameters.AddValue("max_number_of_searchs",self.params["max_number_of_searchs"])
        mmg_parameters.AddValue("echo_level",self.params["echo_level"])
        mmg_parameters.AddValue("save_mdpa_file",self.params["save_mdpa_file"])
        if (self.dim == 2):
            self.MmgProcess = MeshingApplication.MmgProcess2D(self.main_model_part, mmg_parameters)
        else:
            self.MmgProcess = MeshingApplication.MmgProcess3D(self.main_model_part, mmg_parameters)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass
             
    def ExecuteFinalizeSolutionStep(self):
        self.__error_calculation()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        if (self.main_model_part.ProcessInfo[MeshingApplication.ERROR_ESTIMATE] > self.error_threshold):
            self.__execute_refinement()
        self.remeshing_cycle += 1
        if (self.main_model_part.ProcessInfo[MeshingApplication.ERROR_ESTIMATE] <= self.error_threshold or self.remeshing_cycle > self.params["max_iterations"].GetInt()):
            self.main_model_part.ProcessInfo[MeshingApplication.EXECUTE_REMESHING] = False

    def ExecuteFinalize(self):
        pass

    def __create_metric_process(self):
        spr_parameters = KratosMultiphysics.Parameters("""{}""")
        spr_parameters.AddValue("minimal_size",self.params["minimal_size"])
        spr_parameters.AddValue("maximal_size",self.params["maximal_size"])
        spr_parameters.AddValue("error",self.params["error_parameters"]["interpolation_error"])
        spr_parameters.AddValue("echo_level", self.params["echo_level"])
        spr_parameters.AddValue("set_number_of_elements", self.params["set_number_of_elements"])
        spr_parameters.AddValue("number_of_elements", self.params["number_of_elements"])
        spr_parameters.AddValue("average_nodal_h", self.params["average_nodal_h"])
            
        if (self.dim == 2):
            self.metric_process = MeshingApplication.SPRMetricProcess2D(
                self.main_model_part, 
                spr_parameters)
        else:
            self.metric_process = MeshingApplication.SPRMetricProcess3D(
                self.main_model_part, 
                spr_parameters)                     

    def __execute_refinement(self):

        print("Remeshing")
        self.MmgProcess.Execute()

        if (self.params["debug_mode"].GetBool() == True):
            step = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS]
            self._debug_output(step, "")

        # We need to set that the model part has been modified (later on we will act in consequence)
        self.main_model_part.Set(KratosMultiphysics.MODIFIED, True)
        print("Remesh finished")

    def __error_calculation(self):

        # Initialize metric
        self.initialize_metric.Execute()

        print("Calculating the metrics")
        # Execute metric computation
        self.metric_process.Execute()
        self.estimated_error = self.main_model_part.ProcessInfo[MeshingApplication.ERROR_ESTIMATE]
        
        
    def __generate_submodelparts_list_from_input(self, param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve submodelparts name from input (a string) and request the corresponding C++ object to the kernel
        return [self.main_model_part.GetSubModelPart(param[i].GetString()) for i in range(0, param.size())]

    def __generate_variable_list_from_input(self, param):
      '''Parse a list of variables from input.'''
      # At least verify that the input is a string
      if not param.IsArray():
          raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

      # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel

      variable_list = []

      for i in range( 0,param.size()):
          aux_var = KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() )
          for node in self.main_model_part.Nodes:
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
        gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
        singlefile = KratosMultiphysics.MultiFileFlag.SingleFile
        deformed = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions
        gid_io = KratosMultiphysics.GidIO("REMESHING_"+name+"_STEP_"+str(label), gid_mode, singlefile, deformed, write_conditions)
        
        gid_io.InitializeMesh(label)
        gid_io.WriteMesh(self.main_model_part.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(label, self.main_model_part.GetMesh())
        if (self.params["framework"].GetString() ==  "Lagrangian"):
            gid_io.WriteNodalResults(KratosMultiphysics.DISPLACEMENT, self.main_model_part.Nodes, label, 0)
        else:
            gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY, self.main_model_part.Nodes, label, 0)
        gid_io.FinalizeResults()
        
        #raise NameError("DEBUG")
