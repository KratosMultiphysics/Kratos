"""This module contains the available structural response functions and their base class"""
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import structural_mechanics_analysis

import time as timer

def _GetModelPart(model, solver_settings):
    #TODO can be removed once model is fully available
    model_part_name = solver_settings["model_part_name"].GetString()
    if not model.HasModelPart(model_part_name):
        model_part = model.CreateModelPart(model_part_name, 2)
        domain_size = solver_settings["domain_size"].GetInt()
        if domain_size < 0:
            raise Exception('Please specify a "domain_size" >= 0!')
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
    else:
        model_part = model.GetModelPart(model_part_name)

    return model_part


# ==============================================================================
class DirectSensitivityAnalysis(object):
    """Linear static direct sensitivity response function.
    - runs the primal analysis (writes the primal results to an .h5 file)
    - reads the primal results from the .h5 file into the direct sensitivity model part
    - uses primal results to calculate value
    - uses primal results to calculate gradient by running the direct sensitivity analysis

    Attributes
    ----------
    primal_analysis : Primal analysis object of the response function
    direct_sensitivity_analysis : Direct sensitivity analysis object of the response function
    """

    def RunCalculation(self, calculate_gradient):
        self.Initialize()
        self.InitializeSolutionStep()
        if calculate_gradient:
            self.CalculateGradient()
        self.FinalizeSolutionStep()
        self.Finalize()        
  

    def __init__(self, direct_settings, model):        
        self.direct_settings = direct_settings

        # Create the primal solverstatic
        with open(self.direct_settings["primal_settings"].GetString(),'r') as parameter_file:
            primal_parameters = Parameters( parameter_file.read() )
        
        self.primal_model_part = _GetModelPart(model, primal_parameters["solver_settings"])        

        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, primal_parameters)

        # Create the direct sensitivity solver
        direct_sensitivity_parameters = self._GetDirectSensitivityParameters() 
        print(str(direct_sensitivity_parameters))              
        direct_sensitivity_model = KratosMultiphysics.Model()
        self.direct_sensitivity_model_part = _GetModelPart(direct_sensitivity_model, direct_sensitivity_parameters["solver_settings"])

        # TODO find out why it is not possible to use the same model_part
        self.direct_sensitivity_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(direct_sensitivity_model, direct_sensitivity_parameters)

        self.primal_state_variables = [KratosMultiphysics.DISPLACEMENT]
        if primal_parameters["solver_settings"].Has("rotation_dofs"):
            if primal_parameters["solver_settings"]["rotation_dofs"].GetBool():
                self.primal_state_variables.append(KratosMultiphysics.ROTATION)

        

    def Initialize(self):
        
        self.primal_analysis.Initialize()
        self.direct_sensitivity_analysis.Initialize()
        self._InitializePostProcessAndResponseFunctions()

    def InitializeSolutionStep(self):

        # Run the primal analysis.
        # TODO if primal_analysis.status==solved: return
        Logger.PrintInfo("\n> Starting primal analysis")
        startTime = timer.time()
        if not self.primal_analysis.time < self.primal_analysis.end_time:
            self.primal_analysis.end_time += 1
        self.primal_analysis.RunSolutionLoop()
        Logger.PrintInfo("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")

        # TODO the response value calculation for stresses currently only works on the adjoint modelpart
        # this needs to be improved, also the response value should be calculated on the PRIMAL modelpart!!
        self.direct_sensitivity_analysis.time = self.direct_sensitivity_analysis._GetSolver().AdvanceInTime(self.direct_sensitivity_analysis.time)
        
        # Put primal solution on direct sensitivity model - for "auto" setting, else it has to be done by the user e.g. using hdf5 process
        if self.direct_settings["direct_sensitivity_settings"].GetString() == "auto":
            Logger.PrintInfo("> Transfer primal state to direct sensitivity model part.")
            variable_utils = KratosMultiphysics.VariableUtils()
            for variable in self.primal_state_variables:
                variable_utils.CopyModelPartNodalVar(variable, self.primal_model_part, self.direct_sensitivity_model_part, 0)
        
        self.direct_sensitivity_analysis.InitializeSolutionStep()  

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting direct sensitivity analysis")
        startTime = timer.time()
        self.direct_sensitivity_analysis._GetSolver().Predict()
        self.direct_sensitivity_analysis._GetSolver().SolveSolutionStep()
        Logger.PrintInfo("> Time needed for solving the direct sensitivity analysis = ",round(timer.time() - startTime,2),"s")
        for i in self.response_list:        
            self.direct_sensitivity_postprocess.UpdateSensitivities(i)

    def FinalizeSolutionStep(self):
        self.direct_sensitivity_analysis.FinalizeSolutionStep()
        self.direct_sensitivity_analysis.OutputSolutionStep()


    def Finalize(self):
        self.primal_analysis.Finalize()
        self.direct_sensitivity_analysis.Finalize()    
    

    def _GetDirectSensitivityParameters(self):

        direct_sensitivity_settings = self.direct_settings["direct_sensitivity_settings"].GetString()

        if direct_sensitivity_settings == "auto":
            Logger.PrintInfo("\n> Automatic set up direct sensitivity parameters:")

            with open(self.direct_settings["primal_settings"].GetString(),'r') as parameter_file:
                primal_parameters = Parameters( parameter_file.read() )

            #check that HDF5 process is not there
            if primal_parameters["processes"].Has("list_other_processes"):
                for i in range(0,primal_parameters["processes"]["list_other_processes"].size()):
                    process = primal_parameters["processes"]["list_other_processes"][i]
                    raise Exception("Auto setup of direct sensitivity parameters does not support {} in list_other_processes".format(process["python_module"].GetString()))

            # clone primal settings as base for direct sensitivity
            direct_sensitivity_parameters = primal_parameters.Clone()

            # direct sensitivity settings
            solver_settings = direct_sensitivity_parameters["solver_settings"]
            primal_solver_type = solver_settings["solver_type"].GetString()
            if primal_solver_type != "static":
                raise Exception("Auto setup of direct sensitivity parameters does not support {} solver_type. Only available for 'static'".format(primal_solver_type))
            solver_settings["solver_type"].SetString("direct_"+primal_solver_type)

            if not solver_settings.Has("compute_reactions"):
                solver_settings.AddEmptyValue("compute_reactions")
            solver_settings["compute_reactions"].SetBool(False)

            if not solver_settings.Has("move_mesh_flag"):
                solver_settings.AddEmptyValue("move_mesh_flag")
            solver_settings["move_mesh_flag"].SetBool(False)

            if not solver_settings.Has("scheme_settings"):
                tmp = solver_settings.AddEmptyValue("scheme_settings")
                if not tmp.Has("scheme_type"):
                    tmp.AddEmptyValue("scheme_type")
            solver_settings["scheme_settings"]["scheme_type"].SetString("direct_structural")

            # Dirichlet conditions: change variables
            for i in range(0,primal_parameters["processes"]["constraints_process_list"].size()):
                process = direct_sensitivity_parameters["processes"]["constraints_process_list"][i]
                variable_name = process["Parameters"]["variable_name"].GetString()
            process["Parameters"]["variable_name"].SetString("ADJOINT_"+variable_name)

            # Neumann conditions - do not modify to read the same load values as in primal:

            # Output process:
            # TODO how to add the output process? How find out about the variables?
            if direct_sensitivity_parameters.Has("output_configuration"):
                Logger.PrintInfo("> Output process is removed for direct sensitivity analysis. To enable it define direct sensitivity parameters yourself.")
                direct_sensitivity_parameters.RemoveValue("output_configuration")

            # variable settings
            direct_sensitivity_parameters["solver_settings"].AddValue("variable_settings", self.direct_settings["variable_settings"])            

        else: # direct_sensitivity parameters file is explicitely given - do not change it.
            with open(self.direct_settings["direct_sensitivity_settings"].GetString(),'r') as parameter_file:
                direct_sensitivity_parameters = Parameters( parameter_file.read() )

        return direct_sensitivity_parameters

    def _GetVariable(self):
        return self.direct_sensitivity_analysis._GetSolver().variable

    def _InitializePostProcessAndResponseFunctions(self):   

        # Initialize list of response functions
        self.response_list = list()

        direct_response_settings = self.direct_settings["response_function_settings"].Clone()
        direct_sensitivity_settings = self.direct_settings["sensitivity_settings"].Clone() 
        
        for i in range(0, direct_response_settings["response_functions"]["local_stress"].size()):            
            variable_name = direct_response_settings["response_functions"]["local_stress"][i].GetString()            
            response_function = StructuralMechanicsApplication.DirectSensitivityLocalStressResponseFunction(self.direct_sensitivity_model_part, direct_response_settings, variable_name)
            self.response_list.append(response_function)

        for i in range(0, direct_response_settings["response_functions"]["nodal_displacement"].size()):            
            variable_name = direct_response_settings["response_functions"]["nodal_displacement"][i].GetString()            
            response_function = StructuralMechanicsApplication.DirectSensitivityNodalDisplacementResponseFunction(self.direct_sensitivity_model_part, direct_response_settings, variable_name)
            self.response_list.append(response_function)
        
        # Initialize the postprocess
        self.direct_sensitivity_postprocess = StructuralMechanicsApplication.DirectSensitivityPostprocess(self.direct_sensitivity_model_part, self._GetVariable(), direct_sensitivity_settings)
        self.direct_sensitivity_postprocess.Initialize()            

        
    
    
        

