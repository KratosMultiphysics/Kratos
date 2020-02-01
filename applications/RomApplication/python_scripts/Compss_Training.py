from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#importing base classes
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

#import processes to vary parameters ( entries of the vector mu)
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess
from KratosMultiphysics.ConvectionDiffusionApplication.apply_thermal_face_process import ApplyThermalFaceProcess

#import the SVD packages 
from scipy import linalg
from KratosMultiphysics.RomApplication.RSVDT_Library import rsvdt

#import the concurrent package
import concurrent.futures

# Import other packages
import numpy as np
import h5py
import json
import time
import sys
from matplotlib import pyplot as plt

# Import cpickle to pickle the serializer
try:
    import cpickle as pickle  # Use cPickle on Python 2.7
except ImportError:
    import pickle

from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from pycompss.api.task import task
# from pycompss.api.api import compss_wait_on
# from pycompss.api.parameter import *


class ConvectionDiffusionAnalysisWithFlush(ConvectionDiffusionAnalysis):

    def __init__(self,model,project_parameters, sample ):
        super(ConvectionDiffusionAnalysisWithFlush,self).__init__(model,project_parameters)
        self.flush_frequency = 10
        self.last_flush = time.time()
        self.sample = sample

    def FinalizeSolutionStep(self):
        super(ConvectionDiffusionAnalysisWithFlush,self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now


    def ModifyInitialProperties(self):  
        self.processes = []
        
        ################################################################
        #                     TEMPERATURE PROPERTY                     #
        ################################################################
        
        ImposedTemperature = self.sample[0]
        TemperatureSettings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"             : "ThermalModelPart.ImposedTemperature3D_Imposed_temperature_Auto4",
                "variable_name"               : "TEMPERATURE",
                "constrained"                 : true,
                "interval"                    : [0.0,"End"]
            }
            """
            )
        TemperatureSettings.AddEmptyValue("value").SetDouble(ImposedTemperature)
        self.processes.append(AssignScalarVariableProcess(self.model, TemperatureSettings))
        ################################################################



        # ################################################################
        # #                      RADIATION PROPERTY                      #
        # ################################################################
        

        RadiationSettings = KratosMultiphysics.Parameters("""
            {
                "model_part_name": "ThermalModelPart.ThermalFace3D_Thermal_face_conditions_Auto2",
                "add_ambient_radiation": true,
                "emissivity": 0.8,
                "add_ambient_convection": true,
                "convection_coefficient": 100.0
            }
            """
            )

        AmbientTemperature = self.sample[1]
        RadiationSettings.AddEmptyValue("ambient_temperature").SetDouble(AmbientTemperature)
        self.processes.append(ApplyThermalFaceProcess(self.model, RadiationSettings))

        

        RadiationSettings = KratosMultiphysics.Parameters("""
            {
                "model_part_name": "ThermalModelPart.ThermalFace3D_Thermal_face_conditions_Auto3",
                "add_ambient_radiation": true,
                "emissivity": 0.8,
                "add_ambient_convection": true,
                "convection_coefficient": 100.0
            }
            """
            )

        AmbientTemperature = self.sample[2]
        RadiationSettings.AddEmptyValue("ambient_temperature").SetDouble(AmbientTemperature)
        self.processes.append(ApplyThermalFaceProcess(self.model, RadiationSettings))



        ################################################################
        for process in self.processes:
            process.ExecuteBeforeSolutionLoop()
        ################################################################


    def EvaluateQuantityOfInterest(self):
        ##############################################################################################
        # Functions evaluating the QoI of the problem: Array of temperature at every node on the mesh #
        #                    and nodal area                                                           #
        ##############################################################################################
        ArrayOfTemperatures = []
        #for node in self._GetSolver().GetComputingModelPart().Nodes:
        for node in self._GetSolver().model["ThermalModelPart"].Nodes:       
            ArrayOfTemperatures.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        return ArrayOfTemperatures        


    def SaveHyperReductionDataInJSON(self):
        OriginalNumberOfElements = self.model["ThermalModelPart"].NumberOfElements()
        OriginalNumberOfConditions = self.model["ThermalModelPart"].NumberOfConditions()
        
        HRData = {}
        HRData["OriginalNumberOfElements"] = OriginalNumberOfElements
        HRData["OriginalNumberOfConditions"] = OriginalNumberOfConditions
        with open('HRData.json', 'w' ) as f:
            json.dump(HRData,f, indent=2)    


        ##############################################################################################


###############################################################################################################################################################################


"""
function executing a set of instances of the problem
input:
        pickled_model:      serialization of the model
        pickled_parameters: serialization of the Project Parameters
        Cases:              list of values for the variables of interest
output:
        QoI: U_hat = U @ Sigma
"""

@ExaquteTask(returns=1)
def Get_Basis_From_Simulations(pickled_model, pickled_parameters, Cases):
    #Run a batch of simulations, and obtain its basis
    qoi_loop = []
    #### When using COMPSs, only send one entry in Case per time-dependent scenario
    for sample in Cases:
        qoi_loop.append(Single_Simulation(pickled_model,pickled_parameters,sample))

    # Build snapshot matrix        
    SnapshotMatrix = np.zeros((len(qoi_loop[0]),len(qoi_loop)))
    for i in range (len(qoi_loop)):
        Snapshot_i=np.array(qoi_loop[i])
        SnapshotMatrix[:,i]=Snapshot_i.transpose()
    DATA = {}
    DATA['TypeOfSVD'] = 0   
    u,s,_,_=rsvdt(SnapshotMatrix,0,0,0, DATA)    
    #u, s, _ = linalg.svd(SnapshotMatrix, full_matrices=False)
    diagSigma = np.diag(s)
    u_hat = u @ diagSigma
    del (u, s)
    return u_hat

@ExaquteTask(returns=1)
def Single_Simulation(pickled_model,pickled_parameters,sample):
    model_serializer = pickle.loads(pickled_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    # overwrite the old parameters serializer with the unpickled one
    serialized_parameters = pickle.loads(pickled_parameters)    
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    simulation = ConvectionDiffusionAnalysisWithFlush(current_model, current_parameters, sample)
    simulation.Run()
    QoI = simulation.EvaluateQuantityOfInterest()
    return QoI

@ExaquteTask(returns=1,MatricesList=COLLECTION_IN)
def build_snapshots_matrix(MatricesList):
    for i in range(len(MatricesList)):
        if i ==0:
            concatenated_matrix = MatricesList[i]
        else:            
            concatenated_matrix = np.c_[ concatenated_matrix , MatricesList[i]]
    #concatenated_matrix = np.c_[args]        
    return concatenated_matrix


@ExaquteTask(returns=1,MatricesList=COLLECTION_IN)
def Get_Basis_From_Basis(MatricesList):
    concatenated_matrix = build_snapshots_matrix(MatricesList)
    DATA = {}
    DATA['TypeOfSVD'] = 0   
    u,s,_,_=rsvdt(concatenated_matrix,0,0,0, DATA)
    diagSigma = np.diag(s)
    u_hat = u @ diagSigma
    return u_hat

@ExaquteTask(returns=1,MatricesList=COLLECTION_IN)
def Get_Final_Data_From_Basis(MatricesList):  
    concatenated_matrix = build_snapshots_matrix(MatricesList)
    DATA = {}
    DATA['TypeOfSVD'] = 0   
    u,s,_,_=rsvdt(concatenated_matrix,0,0,0, DATA)
    return [u,s]


"""
function serializing and pickling the model and the parameters of the problem
input:
        parameter_file_name: path of the Project Parameters file
output:
        pickled_model:      model serializaton
        pickled_parameters: project parameters serialization
"""
def SerializeModelParameters_Task(parameter_file_name):
    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()
    fake_sample = [25.0,25.0,25.0]
    simulation = ConvectionDiffusionAnalysisWithFlush(model,parameters,fake_sample)
    simulation.Initialize()
    serialized_model = KratosMultiphysics.StreamSerializer()
    serialized_model.Save("ModelSerialization",simulation.model)
    serialized_parameters = KratosMultiphysics.StreamSerializer()
    serialized_parameters.Save("ParametersSerialization",simulation.project_parameters)
    # pickle dataserialized_data
    pickled_model = pickle.dumps(serialized_model, 2) # second argument is the protocol and is NECESSARY (according to pybind11 docs)
    pickled_parameters = pickle.dumps(serialized_parameters, 2)
    print("\n","#"*50," SERIALIZATION COMPLETED ","#"*50,"\n")
    return pickled_model,pickled_parameters


def splittig(qoi,Batch, MinimumBatchAdmisible = 'none'):
    #### Setting defaults
    if MinimumBatchAdmisible == 'none':
        MinimumBatchAdmisible = np.ceil(Batch/2)
    #### Creating a new list to save values
    qoi2=[]  
    #### Splitting into smaller batches for SVD computation
    Index = []
    for i in range(len(qoi)):
        Index.append(i) 
    for i,j in zip(Index[0::Batch], Index[Batch-1::Batch]):
        if j+Batch > len(Index)-1:
            if ( len(Index) - j - 1  ) >= MinimumBatchAdmisible and (j < (len(Index)-1)):   
                qoi2.append( Get_Basis_From_Basis(qoi[i:j+1]))
                qoi2.append( Get_Basis_From_Basis(qoi[j+1:]))
            else:
                qoi2.append(Get_Basis_From_Basis(qoi[i:]))
        else:
            qoi2.append(Get_Basis_From_Basis(qoi[i:j+1]))
    return qoi2



def main():
    # set the ProjectParameters.json path
    parameter_file_name = "/home/jrbravo/Desktop/PhD/Hyper-Reduction/Siemens_Project/SiemensExample/Training_Phi/ProjectParameters.json"
    # create a serialization of the model and of the project parameters
    pickled_model, pickled_parameters = SerializeModelParameters_Task(parameter_file_name)
    
    ##########################################
    #           Cases to train
    #########################################
    #Values for Imposed Temperature = []
    ImpTemp = [273.0, 278.0, 283.0, 288.0, 293.0, 298.0]
    #Values for FACE_HEAT_FLUX1 = []
    AmTemp1 = [273.0, 293.0, 313.0, 333.0]
    #Values for FACE_HEAT_FLUX2 = []
    AmTemp2 = [333.0, 353.0, 373.0, 393.0]
    
    Cases=[]
    for j in range (0,len(ImpTemp)):
        for k in range (0,len(AmTemp1)):
            for l in range (0,len(AmTemp2)):
                Cases.append([ ImpTemp[j], AmTemp1[k], AmTemp2[l] ] )

    ###############################################################################
    ######## Creating a subset of cases to send to each CPU  ######################
    qoi = []
    qoi2 = []
    TotalNumberOFCases = 33#len(Cases)
    print(TotalNumberOFCases)
    i = 0
    j= 0
    SliceOfCases = 3   ## Number of cases to send to each CPU
    MinimumSizeOfCases = np.ceil(SliceOfCases/2)
    ###############################################################################

################################################################################## Reduction Tree ##################################################################    

    #### Run Simulations and Get Basis
    while i < TotalNumberOFCases and j<TotalNumberOFCases:
        if i+SliceOfCases < TotalNumberOFCases and (TotalNumberOFCases-i) > (SliceOfCases + MinimumSizeOfCases):
            j = i+SliceOfCases
            qoi.append(Get_Basis_From_Simulations (pickled_model, pickled_parameters, Cases[i:j]))
        else:
            j = TotalNumberOFCases
            qoi.append(Get_Basis_From_Simulations (pickled_model, pickled_parameters, Cases[i:j]))
        print(i, j)
        i += SliceOfCases

    ##### Get Basis from other Basis (Setting a mini batch, to send multiple Snapshot matrices to a task)
    MiniBatchSVDs = 2    
    while MiniBatchSVDs<len(qoi):
        print('entering split by batches')
        qoi = splittig(qoi, MiniBatchSVDs)

    ##### Get singular values from Basis
    result = Get_Final_Data_From_Basis(qoi) 
    result = compss_wait_on(result)

####################################################################################################################################################################      

    return result[0],result[1]


if __name__ == '__main__':
    
    u,s = main()

    plt.plot( np.linspace(1, len(s), len(s)), s, 'bo-')
    plt.title('Singular Values')
    plt.ylabel('Log scale')
    plt.show()
    np.savetxt('LeftSingularVectors_new', u)
    np.savetxt('SingularValues_new', s)


