# Import Kratos
import KratosMultiphysics
from KratosMultiphysics.MultilevelMonteCarloApplication.adaptive_refinement_utilities import AdaptiveRefinement
from KratosMultiphysics.MultilevelMonteCarloApplication.tools import ParametersWrapper

# Import Kratos problem
from simulation_definition import SimulationScenario
import xmc.classDefs_solverWrapper.methodDefs_KratosSolverWrapper.solve as mds
import xmc.classDefs_solverWrapper.methodDefs_KratosSolverWrapper.utilities as mdu

# Import xmc classes
import xmc.solverWrapper as sw

# Math utilities
import math

# Import cpickle to pickle the serializer
import pickle

# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

class KratosSolverWrapper(sw.SolverWrapper):

    # TODO: solverWrapperIndex will be removed from here and will have an indicator about the level we are at and not which algorithm we are using
    def __init__(self,**keywordArgs):
        super().__init__(**keywordArgs)
        self.analysis = SimulationScenario
        self.adaptive_refinement_jump_to_finest_level = keywordArgs.get("adaptiveRefinementJumpToFinestLevel",False)
        self.asynchronous = keywordArgs.get("asynchronous",False)
        self.different_tasks = not keywordArgs.get('taskAllAtOnce',False)
        self.fake_sample_to_serialize =  keywordArgs.get('fakeRandomVariable')
        self.mapping_output_quantities = keywordArgs.get("mappingOutputQuantities",False)
        self.number_contributions_per_instance = keywordArgs.get("numberContributionsPerInstance",1)
        self.number_qoi = keywordArgs.get("numberQoI",1)
        self.number_combined_qoi = keywordArgs.get("numberCombinedQoi",0)
        self.outputBatchSize = keywordArgs.get('outputBatchSize',1)
        self.print_to_file = keywordArgs.get("printToFile",False)
        self.project_parameters_path = keywordArgs.get('projectParametersPath')
        self.refinement_parameters_path = keywordArgs.get('refinementParametersPath')
        self.refinement_strategy = keywordArgs.get('refinementStrategy')

        # Set outputDimension, as per SolverWrapper specifications (see doc)
        self.outputDimension = keywordArgs.get('outputDimension',None)
        # If not given, compute from self.outputBatchSize for backward compatibility
        if self.outputDimension is None:
            outputNb = self.number_qoi + self.number_combined_qoi
            # Total number of output splits, including (possibly) a last one of smaller size
            batchNb = int(math.ceil(outputNb/self.outputBatchSize))
            # Assemble the list of sizes of each split
            # They are all equal to outputBatchSize, except perhaps the last one
            # E.g. outputBatchSize=2 and outputNb=5 gives [2,2,1]
            self.outputDimension = [min(self.outputBatchSize, outputNb-i*self.outputBatchSize)
                                    for i in range(batchNb)]

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

        else:
            raise Exception ("Select KratosMultiphysics refinement stategy.\nOptions:\
                \n   i) stochastic_adaptive_refinement\
                \n  ii) deterministic_adaptive_refinement\
                \n iii) reading_from_file")


    def solve(self,random_variable):
        if all([component>=0 for component in self.solverWrapperIndex]):
            aux_qoi_array = []
            for contribution_counter in range (0,self.number_contributions_per_instance):
                self.current_local_contribution = contribution_counter
                if (self.refinement_strategy == "stochastic_adaptive_refinement"):
                    qoi,time_for_qoi = self.executeInstanceStochasticAdaptiveRefinement(random_variable)
                elif (self.refinement_strategy == "deterministic_adaptive_refinement"):
                    qoi,time_for_qoi = self.executeInstanceDeterministicAdaptiveRefinement(random_variable)
                elif (self.refinement_strategy == "reading_from_file"):
                    qoi,time_for_qoi = self.executeInstanceReadingFromFile(random_variable)
                # append components to aux array
                aux_qoi_array.append(qoi)
            # delete COMPSs future objects no longer needed
            delete_object(random_variable)

            # postprocess components
            if self.number_contributions_per_instance > 1:
                ppm = mdu.PostProcessManager(self.numberOfOutputs(),self.outputBatchSize)
                if (self.numberOfOutputs() == self.outputBatchSize):
                    qoi_list = [ppm.PostprocessContributionsPerInstance(aux_qoi_array,self.number_qoi,self.number_combined_qoi)]
                else:
                    qoi_list = ppm.PostprocessContributionsPerInstance(aux_qoi_array,self.number_qoi,self.number_combined_qoi)
                delete_object(ppm)
            else:
                # unfold qoi into its components of fixed size
                unm = mdu.UnfolderManager(self.numberOfOutputs(), self.outputBatchSize)
                if (self.numberOfOutputs() == self.outputBatchSize):
                    qoi_list = [unm.UnfoldNValues_Task(aux_qoi_array[0])]
                else:
                    qoi_list = unm.UnfoldNValues_Task(aux_qoi_array[0])
                # delete COMPSs future objects no longer needed
                delete_object(unm)

            # delete COMPSs future objects no longer needed
            for contribution_counter in range (0,self.number_contributions_per_instance):
                delete_object(aux_qoi_array[contribution_counter])
            delete_object(qoi)
            del(aux_qoi_array)

        else:
            qoi,time_for_qoi = mds.returnZeroQoiAndTime_Task(self.numberOfOutputs())
            # unfold qoi into its components of fixed size
            unm = mdu.UnfolderManager(self.numberOfOutputs(), self.outputBatchSize)
            if (self.numberOfOutputs() == self.outputBatchSize):
                qoi_list = [unm.UnfoldNValues_Task(qoi)]
            else:
                qoi_list = unm.UnfoldNValues_Task(qoi)
            # delete COMPSs future objects no longer needed
            delete_object(unm)

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
        current_local_contribution = self.current_local_contribution
        time_for_qoi = 0.0
        if (different_tasks is False): # tasks all at once
            qoi,time_for_qoi = \
                mds.executeInstanceStochasticAdaptiveRefinementAllAtOnce_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_analysis,time_for_qoi,mapping_flag,adaptive_refinement_jump_to_finest_level,print_to_file,current_local_contribution)
        elif (different_tasks is True): # multiple tasks
            if (current_index == 0): # index = 0
                current_local_index = 0
                qoi,pickled_current_model,time_for_qoi = \
                    mds.executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file,current_local_contribution)
                delete_object(pickled_current_model)
            else: # index > 0
                for current_local_index in range(current_index+1):
                    if ((adaptive_refinement_jump_to_finest_level is False) or (adaptive_refinement_jump_to_finest_level is True and (current_local_index == 0 or current_local_index == current_index))):
                        if (mapping_flag is False):
                            qoi,pickled_current_model,time_for_qoi = \
                                mds.executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file,current_local_contribution)
                        elif (mapping_flag is True):
                            qoi,pickled_current_model,time_for_qoi = \
                                mds.executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(current_index,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,random_variable,current_local_index,current_analysis,time_for_qoi,mapping_flag,print_to_file,current_local_contribution,pickled_mapping_reference_model=pickled_reference_model_mapping)
                            delete_object(pickled_coarse_model)
                            del(pickled_coarse_model)
                        pickled_coarse_model = pickled_current_model
                        del(pickled_current_model)
                    else: # not running since we jump from coarsest to finest level
                        pass
                delete_object(pickled_coarse_model)
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
        pickled_model = self.pickled_model[current_index]
        pickled_mapping_reference_model = self.pickled_model[0]
        pickled_project_parameters = self.pickled_project_parameters[current_index]
        mapping_flag = self.mapping_output_quantities
        print_to_file = self.print_to_file
        current_local_contribution = self.current_local_contribution
        current_analysis = self.analysis
        time_for_qoi = 0.0
        # TODO - Change this to be more general
        qoi,time_for_qoi = mds.executeInstanceReadingFromFile_Wrapper(current_index,pickled_model,pickled_project_parameters,current_analysis,random_variable,time_for_qoi,mapping_flag,pickled_mapping_reference_model,print_to_file,current_local_contribution)
        return qoi,time_for_qoi


    ####################################################################################################
    ####################################### SERIALIZATION TOOLS ########################################
    ####################################################################################################


    """
    function serializing and pickling the model and the project parameters of the problem
    the serialization-pickling process is the following:
    i)   from Model/Parameters Kratos object to MpiSerializer Kratos object
    ii)  from MpiSerializer Kratos object to pickle string
    iii) from pickle string to MpiSerializer Kratos object
    iv)  from MpiSerializer Kratos object to Model/Parameters Kratos object
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
    #     serialized_model = KratosMultiphysics.MpiSerializer()
    #     serialized_model.Save("ModelSerialization",simulation.model)
    #     serialized_project_parameters = KratosMultiphysics.MpiSerializer()
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

        if (self.refinement_strategy == "stochastic_adaptive_refinement"):
            self.SerializeModelParametersStochasticAdaptiveRefinement()
        elif (self.refinement_strategy == "deterministic_adaptive_refinement"):
            self.SerializeModelParametersDeterministicAdaptiveRefinement()
        elif (self.refinement_strategy == "reading_from_file"):
            self.SerializeModelParametersReadingFromFile()
        else:
            raise Exception ("Specify refinement_strategy: stochastic_adaptive_refinement or deterministic_adaptive_refinement or reading_from_file")
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
        serialized_project_parameters = KratosMultiphysics.MpiSerializer()
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
        serialized_model = KratosMultiphysics.MpiSerializer()
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
                fake_qoi,pickled_current_model,fake_computational_time = \
                    mds.executeInstanceStochasticAdaptiveRefinementMultipleTasks_Wrapper(number_levels_to_serialize,pickled_coarse_model,pickled_coarse_project_parameters,pickled_custom_metric_refinement_parameters,pickled_custom_remesh_refinement_parameters,fake_sample,current_level,current_analysis,fake_computational_time,mapping_flag=False,print_to_file=False,current_contribution=0)
                del(pickled_coarse_model)
                pickled_coarse_model = pickled_current_model
                # save if current level > 0 (level = 0 has already been saved)
                if (current_level>0):
                    # save pickled and serialized model and parameters
                    self.pickled_model.append(pickled_current_model)
                    self.serialized_model.append(pickle.loads(pickled_current_model))
                del(pickled_current_model)

    def SerializeModelParametersReadingFromFile(self):
        for parameters_path in self.project_parameters_path:
            with open(parameters_path,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            # create wrapper instance to modify current project parameters
            self.wrapper = ParametersWrapper(parameters)
            # serialize parmeters (to avoid adding new data dependent on the application)
            parameters = self.wrapper.SetModelImportSettingsInputType("use_input_model_part")
            serialized_project_parameters = KratosMultiphysics.MpiSerializer()
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
            serialized_model = KratosMultiphysics.MpiSerializer()
            serialized_model.Save("ModelSerialization",simulation.model)
            self.serialized_model.append(serialized_model)
            # pickle dataserialized_data
            pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
            pickled_project_parameters = pickle.dumps(serialized_project_parameters, 2)
            self.pickled_model.append(pickled_model)
            self.pickled_project_parameters.append(pickled_project_parameters)

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
        # save parameters as MpiSerializer Kratos objects
        serialized_metric_refinement_parameters = KratosMultiphysics.MpiSerializer()
        serialized_metric_refinement_parameters.Save("MetricRefinementParametersSerialization",metric_refinement_parameters)
        serialized_remesh_refinement_parameters = KratosMultiphysics.MpiSerializer()
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


    def numberOfOutputs(self):
        """
        Returns the total number of _scalar_ outputs.
        """
        if isinstance(self.outputDimension,int):
            return self.outputDimension
        else: #is a list of integers
            return sum(self.outputDimension)
