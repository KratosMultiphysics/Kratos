import numpy as np
import tensorflow as tf

import KratosMultiphysics as KMP

tf.keras.backend.set_floatx('float64')

class StructuralMechanics_NN_Interface():

    def __init__(self, kratos_simulation):

        self.fake_simulation = kratos_simulation
        
        self.fake_simulation.Initialize()
        self.fake_simulation.InitializeSolutionStep()

        self.space = KMP.UblasSparseSpace()
        self.strategy = self.fake_simulation._GetSolver()._GetSolutionStrategy()
        self.buildsol = self.fake_simulation._GetSolver()._GetBuilderAndSolver()
        self.scheme = self.fake_simulation._GetSolver()._GetScheme()
        self.modelpart = self.fake_simulation._GetSolver().GetComputingModelPart()
        self.var_utils = KMP.VariableUtils()

    def project_prediction_vectorial_optim_batch(self, y_pred):
        values = y_pred
        values_full=np.zeros(self.modelpart.NumberOfNodes()*2)
        values_full+=values

        dim = 2
        nodes_array=self.modelpart.Nodes
        x0_vec = self.var_utils.GetInitialPositionsVector(nodes_array,dim)
        self.var_utils.SetSolutionStepValuesVector(nodes_array, KMP.DISPLACEMENT, values_full, 0)
        x_vec=x0_vec+values_full
        self.var_utils.SetCurrentPositionsVector(nodes_array,x_vec)
    
    def get_v_loss_rdiff_batch_(self, y_pred, b_true):
        A = self.strategy.GetSystemMatrix()
        b = self.strategy.GetSystemVector()

        err_r_list=[]
        v_loss_r_list=[]

        for i in range(y_pred.shape[0]):
            self.space.SetToZeroMatrix(A)
            self.space.SetToZeroVector(b)

            self.project_prediction_vectorial_optim_batch(y_pred[i])

            self.buildsol.Build(self.scheme, self.modelpart, A, b)

            err_r=KMP.Vector(b_true[i].copy()-b)

            v_loss_r = self.space.CreateEmptyVectorPointer()
            self.space.ResizeVector(v_loss_r, self.space.Size(b))
            self.space.SetToZeroVector(v_loss_r)

            self.space.TransposeMult(A,err_r,v_loss_r)
            
            err_r_list.append(np.expand_dims(np.array(err_r, copy=False),axis=0))
            v_loss_r_list.append(np.expand_dims(np.array(v_loss_r, copy=False),axis=0))
            # The negative sign we should apply to A is compensated by the derivative of the loss

        err_r_batch = np.concatenate(err_r_list, axis = 0)
        v_loss_r_batch = np.concatenate(v_loss_r_list, axis = 0)
        
        return err_r_batch, v_loss_r_batch
    
    def get_err_rdiff_batch_(self, y_pred, b_true):
        b_pred = self.get_r_batch_(y_pred)
        err_r = b_pred - b_true
        return err_r
    
    def get_r_batch_(self, y_pred):

        A = self.strategy.GetSystemMatrix()
        b = self.strategy.GetSystemVector()

        b_list=[]

        for i in range(y_pred.shape[0]):
            self.space.SetToZeroMatrix(A)
            self.space.SetToZeroVector(b)

            self.project_prediction_vectorial_optim_batch(y_pred[i])

            self.buildsol.Build(self.scheme, self.modelpart, A, b)

            b_list.append(np.expand_dims(np.array(b, copy=True),axis=0))

        b_batch = np.concatenate(b_list, axis = 0)
        
        return b_batch
    
    @tf.function(input_signature=(tf.TensorSpec(None, tf.float64), tf.TensorSpec(None, tf.float64)))
    def get_v_loss_rdiff_batch(self, y_pred, b_true):
        y,w = tf.numpy_function(self.get_v_loss_rdiff_batch_, [y_pred, b_true], (tf.float64, tf.float64))
        return y,w
    
    @tf.function(input_signature=(tf.TensorSpec(None, tf.float64), tf.TensorSpec(None, tf.float64)))
    def get_err_rdiff_batch(self, y_pred, b_true):
        y = tf.numpy_function(self.get_err_rdiff_batch_, [y_pred, b_true], (tf.float64))
        return y
    
