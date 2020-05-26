# Import Kratos
import KratosMultiphysics
try:  #TODO remove after migration
    from KratosMultiphysics.MultilevelMonteCarloApplication.compressible_adaptive_refinement_utilities import AdaptiveRefinement
except:
    from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement
from KratosMultiphysics.MultilevelMonteCarloApplication.tools import ParametersWrapper

# Import Kratos problem
from simulation_definition import SimulationScenario
import xmc.classDefs_solverWrapper.methodDefs_KratosSolverWrapper.solve as mds
import xmc.classDefs_solverWrapper.methodDefs_KratosSolverWrapper.utilities as mdu

# Import xmc classes
import xmc.solverWrapper as sw

# Import cpickle to pickle the serializer
import pickle

class KratosSolverWrapper(sw.SolverWrapper):

    # TODO: solverWrapperIndex will be removed from here and will have an indicator about the level we are at and not which algorithm we are using
    def __init__(self,**keywordArgs):
        super().__init__(**keywordArgs)
        self.analysis = SimulationScenario
        self.adaptive_refinement_jump_to_finest_level = keywordArgs.get("adaptiveRefinementJumpToFinestLevel",False)
        self.asynchronous = keywordArgs.get("asynchronous",False)
        self.fake_sample_to_serialize =  keywordArgs.get('fakeRandomVariable')
        self.mapping_output_quantities = keywordArgs.get("mappingOutputQuantities",False)
        self.number_contributions_per_instance = keywordArgs.get("numberContributionsPerInstance",2)
        self.number_qoi = keywordArgs.get("numberQoI",1)
        self.number_combined_qoi = keywordArgs.get("numberCombinedQoi",0)
        self.print_to_file = keywordArgs.get("printToFile",False)
        self.project_parameters_path = keywordArgs.get('projectParametersPath')
        self.refinement_parameters_path = keywordArgs.get('refinementParametersPath')
        self.refinement_strategy = keywordArgs.get('refinementStrategy')
        self.different_tasks = not keywordArgs.get('taskAllAtOnce',False)

        # TODO: remove this hard code to run MC
        if (self.solverWrapperIndex == []):
            self.solverWrapperIndex.append(0)

        if (self.solverWrapperIndex[0] >= 0): # for index < 0 nothing is needed
            if (self.asynchronous is not True): # synchronous framework
                self.serialize()
            else: # asynchronous framework
                pass


    def serialize(self):
        if (self.refinement_strategy == "stochastic_adaptive_refinement"):
            # serialization
            self.is_project_parameters_pickled = False
            self.is_model_pickled = False
            self.is_custom_settings_metric_refinement_pickled = False
            self.is_custom_settings_remesh_refinement_pickled = False
            self.SetRefinementParameters()
            self.SerializeRefinementParameters()
            self.SerializeModelParameters()
            # estimate mesh size of current index
            self.ComputeMeshParameters()

        elif (self.refinement_strategy == "deterministic_adaptive_refinement"):
            # serialization
            self.is_project_parameters_pickled = False
            self.is_model_pickled = False
            self.is_custom_settings_metric_refinement_pickled = False
            self.is_custom_settings_remesh_refinement_pickled = False
            self.SetRefinementParameters()
            self.SerializeRefinementParameters()
            self.SerializeModelParameters()
            # estimate mesh size of current index
            self.ComputeMeshParameters()

        elif (self.refinement_strategy == "reading_from_file"):
            # TODO - Change this implementation for later
            # If solverWrapperIndex == [], then projectParametersPath should indicate the absolute file
            # Else, projectParametersPath should be a prefix such that (number).json will be appeneded
            if (self.solverWrapperIndex[0] >= 0):
                # reading of json file
                self.project_parameters_path = self.project_parameters_path+'_'+str(self.solverWrapperIndex[0])+'.json'
                # serialization
                self.is_project_parameters_pickled = False
                self.is_model_pickled = False
                self.SerializeModelParameters()
            else:
                pass

        else:
            raise Exception ("Select KratosMultiphysics refinement stategy.\nOptions:\
                \n   i) stochastic_adaptive_refinement\
                \n  ii) deterministic_adaptive_refinement\
                \n iii) reading_from_file")


    def solve(self,random_variable):
        if all([component>=0 for component in self.solverWrapperIndex]):
            aux_qoi_array = [[] for _ in range (0,self.number_qoi+self.number_combined_qoi)] # to store each qoi
            for _ in range (0,self.number_contributions_per_instance):
                if (self.refinement_strategy == "stochastic_adaptive_refinement"):
                    qoi,time_for_qoi = self.executeInstanceStochasticAdaptiveRefinement(random_variable)
                elif (self.refinement_strategy == "deterministic_adaptive_refinement"):
                    qoi,time_for_qoi = self.executeInstanceDeterministicAdaptiveRefinement(random_variable)
                elif (self.refinement_strategy == "reading_from_file"):
                    qoi,time_for_qoi = self.executeInstanceReadingFromFile(random_variable)
                # unfold qoi into its components
                if ((self.number_qoi + self.number_combined_qoi) == 1):
                    qoi_list = mdu.unfoldVales_Wrapper(self.number_qoi,qoi)
                else:
                    qoi_list = mdu.UnfolderManager(self.number_qoi + self.number_combined_qoi).UnfoldNValues_Task(qoi)
                # append components to aux array
                for qoi_counter in range (0,self.number_qoi+self.number_combined_qoi):
                    aux_qoi_array[qoi_counter].append(qoi_list[qoi_counter])
            # postprocess components
            qoi_postprocessed = mdu.PostprocessContributionsPerInstance(aux_qoi_array,self.number_qoi,self.number_combined_qoi)
            # unfold qoi into its components
            if ((self.number_qoi + self.number_combined_qoi) == 1):
                qoi_list = mdu.unfoldVales_Wrapper(self.number_qoi,qoi_postprocessed)
            else:
                qoi_list = mdu.UnfolderManager(self.number_qoi + self.number_combined_qoi).UnfoldNValues_Task(qoi_postprocessed)
        else:
            qoi,time_for_qoi = mds.returnZeroQoIandTime_Task()
            qoi_list = [qoi]*(self.number_qoi + self.number_combined_qoi)
        return qoi_list,time_for_qoi


    ####################################################################################################
    ######################################### EXECUTION TOOLS ##########################################
    ####################################################################################################


    """
    function executing an instance of the UQ algorithm, i.e. a single MC simulation and eventually the refinement (that occurs before the simulation run)
    input:  self: an instance of the class
    output: qoi: quantity of interest value
    """
    def executeInstanceStochasticAdaptiveRefinement(self,random_variable):
        current_index = self.solverWrapperIndex[0]
        # local variables
        pickled_coarse_model = self.pickled_model[0]
        pickled_reference_model_mapping = pickled_coarse_model
        pickled_coarse_project_parameters = self.pickled_project_parameters[0]
        pickled_custom_metric_refinement_parameters = self.pickled_custom_metric_refinement_parameters
        pickled_custom_remesh_refinement_parameters = self.pickled_custom_remesh_refinement_parameters
        current_analysis = self.analysis
        different_tasks = self.different_tasks
        mapping_flag = self.mapping_output_quantities
        adaptive_refinement_jump_to_finest_level = self.adaptive_refinement_jump_to_finest_level
        print_to_file = self.print_to_file
        time_for_qoi = 0.0
        if (different_tasks is False): # tasks all at once
            qoi,time_for_qoi = \
                mds.executeInstanceStochasticAdaptiveRefinementAllAtOnce_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file)
        elif (different_tasks is True): # multiple tasks
            if (current_index == 0): # index = 0
                current_local_index = 0
                qoi,pickled_current_model,time_for_qoi = \
                    mds.executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file)
            else: # index > 0
                for current_local_index in range(current_index+1):
                    if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
                        if (mapping_flag is False):
                            qoi,pickled_current_model,time_for_qoi = \
                                mds.executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file)
                        elif (mapping_flag is True):
                            qoi,pickled_current_model,time_for_qoi = \
                                mds.executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file,pickled_reference_model_mapping=pickled_mapping_reference_model)
                            del(pickled_coarse_model)
                        pickled_coarse_model = pickled_current_model
                        del(pickled_current_model)
                    else: # not running since we jump from coarsest to finest level
                        pass
        else:
            raise Exception ("Boolean variable different task is not a boolean, instead is equal to",different_tasks)
        return qoi,time_for_qoi

    """
    function executing an instance of the UQ algorithm, i.e. a single MC simulation and eventually the refinement (that occurs before the simulation run)
    input:  self: an instance of the class
    output: qoi: quantity of interest value
    """
    def executeInstanceDeterministicAdaptiveRefinement(self,random_variable):
        # local variables
        current_index = self.solverWrapperIndex[0]
        pickled_model = self.pickled_model[current_index]
        pickled_project_parameters = self.pickled_project_parameters[0]
        current_analysis = self.analysis
        time_for_qoi = 0.0
        qoi,time_for_qoi = mds.executeInstanceDeterministicAdaptiveRefinement_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi)
        return qoi,time_for_qoi

    """
    function executing an instance of the UQ algorithm, i.e. a single MC simulation and eventually the refinement (that occurs before the simulation run)
    input:  self: an instance of the class
    output: qoi: quantity of interest value
    """
    def executeInstanceReadingFromFile(self,random_variable):
        # local variables
        current_index = self.solverWrapperIndex[0]
        pickled_model = self.pickled_model[0]
        pickled_project_parameters = self.pickled_project_parameters[0]
        current_analysis = self.analysis
        time_for_qoi = 0.0
        # TODO - Change this to be more general
        qoi,time_for_qoi = mds.executeInstanceReadingFromFile_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi)
        return qoi,time_for_qoi


    ####################################################################################################
    ####################################### SERIALIZATION TOOLS ########################################
    ####################################################################################################


    """
    function serializing and pickling the model and the project parameters of the problem
    the serialization-pickling process is the following:
    i)   from Model/Parameters Kratos object to StreamSerializer Kratos object
    ii)  from StreamSerializer Kratos object to pickle string
    iii) from pickle string to StreamSerializer Kratos object
    iv)  from StreamSerializer Kratos object to Model/Parameters Kratos object
    requires: self.project_parameters_path: path of the Project Parameters file
    builds: self.pickled_model:              pickled model
            self.pickled_project_parameters: pickled project parameters
    input:  self: an instance of the class
    """
    # def SerializeModelParameters(self):
    #     with open(self.project_parameters_path,'r') as parameter_file:
    #         parameters = KratosMultiphysics.Parameters(parameter_file.read())
    #     model = KratosMultiphysics.Model()
    #     fake_sample = [1.0] # TODO: make application-indipendent
    #     simulation = self.analysis(model,parameters,fake_sample)
    #     simulation.Initialize()
    #     serialized_model = KratosMultiphysics.StreamSerializer()
    #     serialized_model.Save("ModelSerialization",simulation.model)
    #     serialized_project_parameters = KratosMultiphysics.StreamSerializer()
    #     serialized_project_parameters.Save("ParametersSerialization",simulation.project_parameters)
    #     self.serialized_model = serialized_model
    #     self.serialized_project_parameters = serialized_project_parameters
    #     # pickle dataserialized_data
    #     pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    #     pickled_project_parameters = pickle.dumps(serialized_project_parameters, 2)
    #     self.pickled_model = pickled_model
    #     self.pickled_project_parameters = pickled_project_parameters
    #     self.is_project_parameters_pickled = True
    #     self.is_model_pickled = True
    #     print("\n","#"*50," SERIALIZATION MODEL PART AND PARAMETERS COMPLETED ","#"*50,"\n")

    def SerializeModelParameters(self):
        self.serialized_model = []
        self.serialized_project_parameters = []
        self.pickled_model = []
        self.pickled_project_parameters = []

        if (self.refinement_strategy == "stochastic_adaptive_refinement" or self.refinement_strategy == "reading_from_file"):
            self.SerializeModelParametersStochasticAdaptiveRefinement()
        elif (self.refinement_strategy == "deterministic_adaptive_refinement"):
            self.SerializeModelParametersDeterministicAdaptiveRefinement()
        else:
            raise Exception ("Specify refinement_strategy: stochastic_adaptive_refinement or deterministic_adaptive_refinement")
        self.is_project_parameters_pickled = True
        self.is_model_pickled = True
        print("\n","#"*50," SERIALIZATION MODEL AND PROJECT PARAMETERS COMPLETED ","#"*50,"\n")

    def SerializeModelParametersStochasticAdaptiveRefinement(self):
        with open(self.project_parameters_path,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        # create wrapper instance to modify current project parameters
        self.wrapper = ParametersWrapper(parameters)
        # serialize parmeters (to avoid adding new data dependent on the application)
        parameters = self.wrapper.SetModelImportSettingsInputType("use_input_model_part")
        serialized_project_parameters = KratosMultiphysics.StreamSerializer()
        serialized_project_parameters.Save("ParametersSerialization",parameters)
        self.serialized_project_parameters.append(serialized_project_parameters)
        # reset to read the model part
        parameters = self.wrapper.SetModelImportSettingsInputType("mdpa")
        # prepare the model to serialize
        model = KratosMultiphysics.Model()
        fake_sample = self.fake_sample_to_serialize
        simulation = self.analysis(model,parameters,fake_sample)
        simulation.Initialize()
        # reset general flags
        main_model_part_name = self.wrapper.GetModelPartName()
        simulation.model.GetModelPart(main_model_part_name).ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED,True)
        # serialize model
        serialized_model = KratosMultiphysics.StreamSerializer()
        serialized_model.Save("ModelSerialization",simulation.model)
        self.serialized_model.append(serialized_model)
        # pickle dataserialized_data
        pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
        pickled_project_parameters = pickle.dumps(serialized_project_parameters, 2)
        self.pickled_model.append(pickled_model)
        self.pickled_project_parameters.append(pickled_project_parameters)

    def SerializeModelParametersDeterministicAdaptiveRefinement(self):
        self.SerializeModelParametersStochasticAdaptiveRefinement() # to prepare parameters and model part of coarsest level
        number_levels_to_serialize = self.solverWrapperIndex[0]
        # same routine of ExecuteInstanceConcurrentAdaptiveRefinemnt() to build models and parameters, but here we save models and parameters
        pickled_coarse_model = self.pickled_model[0]
        pickled_coarse_project_parameters = self.pickled_project_parameters[0]
        pickled_custom_metric_refinement_parameters = self.pickled_custom_metric_refinement_parameters
        pickled_custom_remesh_refinement_parameters = self.pickled_custom_remesh_refinement_parameters
        current_analysis = self.analysis
        # generate the sample and prepare auxiliary variables we need
        fake_sample = self.fake_sample_to_serialize
        fake_computational_time = 0.0
        if (number_levels_to_serialize > 0):
            for current_level in range(number_levels_to_serialize+1):
                if (current_level < number_levels_to_serialize):
                    fake_qoi,pickled_current_model,fake_computational_time = \
                        mds.executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(number_levels_to_serialize,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,fake_sample,current_level,current_analysis,fake_computational_time)
                        # mds.executeInstanceStochasticAdaptiveRefinementAux_Task(number_levels_to_serialize,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,fake_sample,current_level,current_analysis,fake_computational_time)
                    del(pickled_coarse_model)
                    pickled_coarse_model = pickled_current_model
                elif (current_level == number_levels_to_serialize):
                    pickled_current_model, fake_computational_time = mds.executeInstanceOnlyAdaptiveRefinement_Wrapper(pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,fake_sample,current_level,current_analysis,fake_computational_time)
                # save if current level > 0
                if (current_level>0):
                    # save pickled and serialized model and parameters
                    self.pickled_model.append(pickled_current_model)
                    self.serialized_model.append(pickle.loads(pickled_current_model))
                del(pickled_current_model)

    """
    function serializing and pickling the custom setting metric and remeshing for the refinement
    requires: self.custom_metric_refinement_parameters
              self.custom_remesh_refinement_parameters
    builds: self.pickled_custom_metric_refinement_parameters
            self.pickled_custom_remesh_refinement_parameters
    input:  self: an instance of the class
    """
    def SerializeRefinementParameters(self):
        metric_refinement_parameters = self.custom_metric_refinement_parameters
        remeshing_refinement_parameters = self.custom_remesh_refinement_parameters
        # save parameters as StreamSerializer Kratos objects
        serialized_metric_refinement_parameters = KratosMultiphysics.StreamSerializer()
        serialized_metric_refinement_parameters.Save("MetricRefinementParametersSerialization",metric_refinement_parameters)
        serialized_remesh_refinement_parameters = KratosMultiphysics.StreamSerializer()
        serialized_remesh_refinement_parameters.Save("RemeshRefinementParametersSerialization",remeshing_refinement_parameters)
        # pickle parameters
        pickled_metric_refinement_parameters = pickle.dumps(serialized_metric_refinement_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
        pickled_remesh_refinement_parameters = pickle.dumps(serialized_remesh_refinement_parameters, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
        self.pickled_custom_metric_refinement_parameters = pickled_metric_refinement_parameters
        self.pickled_custom_remesh_refinement_parameters = pickled_remesh_refinement_parameters
        self.is_custom_settings_metric_refinement_pickled = True
        self.is_custom_settings_remesh_refinement_pickled = True
        print("\n","#"*50," SERIALIZATION REFINEMENT PARAMETERS COMPLETED ","#"*50,"\n")


    ####################################################################################################
    ######################################### AUXILIARY TOOLS ##########################################
    ####################################################################################################


    """
    function reading the refinement parameters passed from json file
    input:  self: an instance of the class
    """
    def SetRefinementParameters(self):
        with open(self.refinement_parameters_path,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        self.custom_metric_refinement_parameters = parameters["hessian_metric"]
        self.custom_remesh_refinement_parameters = parameters["refinement_mmg"]

    """
    function giving as output the mesh discretization parameter
    the mesh parameter is the reciprocal of the minimum mesh size of the grid
    input:  self: an instance of the class
    """
    def ComputeMeshParameters(self):
        # unpickle and unserialize model and build Kratos Model object
        serialized_model = pickle.loads(self.pickled_model[0])
        current_model = KratosMultiphysics.Model()
        serialized_model.Load("ModelSerialization",current_model)
        # unpickle and unserialize parameters and build Kratos Parameters object
        serialized_project_parameters = pickle.loads(self.pickled_project_parameters[0])
        current_project_parameters = KratosMultiphysics.Parameters()
        serialized_project_parameters.Load("ParametersSerialization",current_project_parameters)
        # unpickle and unserialize metric refinement parameters and build Kratos Parameters objects
        serialized_custom_metric_refinement_parameters = pickle.loads(self.pickled_custom_metric_refinement_parameters)
        current_custom_metric_refinement_parameters = KratosMultiphysics.Parameters()
        serialized_custom_metric_refinement_parameters.Load("MetricRefinementParametersSerialization",current_custom_metric_refinement_parameters)
        self.mesh_sizes = []
        self.mesh_parameters = []
        level = self.solverWrapperIndex[0]
        adaptive_refinement_manager = AdaptiveRefinement(level,current_model,current_project_parameters,current_custom_metric_refinement_parameters,None)
        adaptive_refinement_manager.EstimateMeshSizeCurrentLevel()
        h_current_level = adaptive_refinement_manager.mesh_size
        mesh_parameter_current_level = h_current_level**(-1)
        self.mesh_sizes.append(h_current_level)
        self.mesh_parameters.append(mesh_parameter_current_level)
