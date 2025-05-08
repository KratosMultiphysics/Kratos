import numpy as np

import tensorflow as tf

class PODANN_prepost_processor():
    def __init__(self):
        super().__init__()
        self.phi_inf = None
        self.phi_inf_tf = None
        self.phi_sup = None
        self.phi_sup_tf = None
        self.sigma_inf = None
        self.sigma_inf_tf = None
        self.sigma_sup = None
        self.sigma_sup_tf = None

    def configure_processor(self, phi, sigma, S, svd_inf_size, svd_sup_size):

        self.phi=phi
        self.sigma=sigma/np.sqrt(S.shape[0])

        self.phi_inf=self.phi[:,:svd_inf_size].copy()
        self.sigma_inf=self.sigma[:svd_inf_size].copy()
        self.phi_sup=self.phi[:,svd_inf_size:svd_sup_size].copy()
        self.sigma_sup=self.sigma[svd_inf_size:svd_sup_size].copy()
        print('Phi_inf matrix shape: ', self.phi_inf.shape)
        print('Sigma_inf array shape: ', self.sigma_inf.shape)
        print('Phi_sgs matrix shape: ', self.phi_sup.shape)
        print('Sigma_sgs array shape: ', self.sigma_sup.shape)
        self.phi_inf_tf=tf.constant(self.phi_inf)
        self.sigma_inf_tf=tf.constant(self.sigma_inf)
        self.phi_sup_tf=tf.constant(self.phi_sup)
        self.sigma_sup_tf=tf.constant(self.sigma_sup)

        ## Check reconstruction error
        S_recons_aux1=self.preprocess_nn_output_data(S)
        S_recons_aux2, _ =self.preprocess_input_data(S)
        S_recons = self.postprocess_output_data(S_recons_aux1, (S_recons_aux2, None))
        print('Reconstruction error SVD (Frob): ', np.linalg.norm(S_recons-S)/np.linalg.norm(S))
        err_aux=np.linalg.norm(S-S_recons, ord=2, axis=1)/np.linalg.norm(S, ord=2, axis=1)
        print('Reconstruction error SVD (Geometric Mean of L2): ', np.exp(np.sum(np.log(err_aux))/S.shape[0]))

    def preprocess_nn_output_data(self, snapshot):
        # Returns q_sup from input snapshots
        output_data=snapshot.copy()
        output_data=np.divide(np.matmul(self.phi_sup.T,output_data.T).T,self.sigma_sup)
        # output_data=np.matmul(self.phi_sup.T,output_data.T).T
        return output_data
    
    def preprocess_nn_output_data_tf(self,snapshot_tensor):
        # Returns q_sup from input snapshots
        output_tensor=tf.transpose(tf.linalg.matmul(self.phi_sup_tf,snapshot_tensor,transpose_a=True,transpose_b=True))/self.sigma_sup_tf
        return output_tensor
    
    def preprocess_input_data(self, snapshot):
        # Returns q_inf from input snapshots
        output_data=snapshot.copy()
        output_data=np.divide(np.matmul(self.phi_inf.T,output_data.T).T,self.sigma_inf)
        # output_data=np.matmul(self.phi_inf.T,output_data.T).T
        return output_data, None
    
    def preprocess_input_data_tf(self, snapshot_tensor):
        # Returns q_inf from input snapshots
        output_tensor=tf.transpose(tf.linalg.matmul(self.phi_inf_tf,snapshot_tensor,transpose_a=True,transpose_b=True))/self.sigma_inf_tf
        return output_tensor, None

    def postprocess_output_data(self, q_sup, aux_norm_data):
        # Returns reconstructed u from given q_inf and q_sup
        q_inf, _ = aux_norm_data
        output_data_1=q_inf.copy()
        output_data_1=np.matmul(self.phi_inf,np.multiply(output_data_1, self.sigma_inf).T).T
        output_data_2=q_sup.copy()
        output_data_2=np.matmul(self.phi_sup,np.multiply(output_data_2, self.sigma_sup).T).T
        output_data = output_data_1 + output_data_2
        return output_data
    
    def postprocess_output_data_tf(self, q_sup_tensor, aux_norm_tensors):
        # Returns reconstructed u from given q_inf and q_sup
        q_inf_tensor, _ = aux_norm_tensors
        output_tensor_1=tf.transpose(tf.linalg.matmul(self.phi_inf_tf,q_inf_tensor*self.sigma_inf_tf,transpose_b=True))
        output_tensor_2=tf.transpose(tf.linalg.matmul(self.phi_sup_tf,q_sup_tensor*self.sigma_sup_tf,transpose_b=True))
        output_tensor = output_tensor_1 + output_tensor_2
        return output_tensor
    
    def get_phi_matrices(self):
        return self.phi_inf, self.phi_sup
    
    def get_sigma_vectors(self):
        return self.sigma_inf, self.sigma_sup
    
    def get_training_data(self, S_train, S_val, R_train, R_val):
        
        val_target_aux=None

        S_train_raw=S_train
        S_target_train=S_train_raw.copy()
        target_aux=R_train

        input_data, _ = self.preprocess_input_data(S_train_raw)
        target_data=(S_target_train, target_aux)

        S_val_raw=S_val
        S_target_val=S_val_raw.copy()
        val_target_aux=R_val

        val_input, _ =self.preprocess_input_data(S_val_raw)
        val_target=(S_target_val, val_target_aux)

        return input_data, target_data, val_input, val_target

