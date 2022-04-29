import os
import sys
import logging

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import json
import math

import contextlib

import h5py
import numpy as np

import keras
import tensorflow as tf

from keras import layers
from itertools import repeat

import kratos_io
import clustering
import networks.shallow_autoencoder as shallow_autoencoder
import networks.cluster_autoencoder as cluster_autoencoder
import networks.gradient_shallow_jacobi    as gradient_shallow

import utils.check_gradient as check_gradient

import matplotlib.pyplot as plt

import KratosMultiphysics as KMP
import KratosMultiphysics.RomApplication as romapp

from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

def print_gpu_info():
    gpu_devices = tf.config.experimental.list_physical_devices('GPU')
    for device in gpu_devices:
        tf.config.experimental.set_memory_growth(device, True)

def custom_loss(y_true, y_pred):
    y_diff = y_true-y_pred
    y_diff = y_diff ** 2

    return  y_diff

def custom_loss_cluster(y_true, y_pred):        
    y_diff = y_true-y_pred
    y_diff = y_diff ** 2

    return  y_diff

if __name__ == "__main__":
    
    data_inputs = [
        "hdf5_output/result_60.h5",
        "hdf5_output/result_65.h5",
        "hdf5_output/result_70.h5"
    ]

    # Save & Train Profile
    config = {
        "load_model":       False,
        "train_model":      True,
        "save_model":       True,
        "print_results":    True,
        "use_reduced":      True,
    }

    # Load Profile
    # config = {
    #     "load_model":       True,
    #     "train_model":      False,
    #     "save_model":       False,
    #     "print_results":    True,
    #     "use_reduced":      True,
    # }

    S = kratos_io.build_snapshot_grid(
        data_inputs, 
        [
            # "PRESSURE", 
            "VELOCITY"
        ]
    )

    kratos_network = gradient_shallow.GradientShallow()
    kratos_io.print_npy_snapshot(S, True)

    print("S", S.shape)

    print("=== Calculating Randomized Singular Value Decomposition ===")
    with contextlib.redirect_stdout(None):
        U,sigma,_,error = RandomizedSingularValueDecomposition().Calculate(S,0)
        SPri = U @ U.T @ S
    
    print(f"Check Norm for SVD  (All    ): {np.linalg.norm(SPri-S)/np.linalg.norm(S)}")
    print(f"Check Norm for SVD  (Col-100): {np.linalg.norm(SPri[5]-S[5])/np.linalg.norm(S[5])}")

    print("U", U.shape)

    if config["use_reduced"]:
        SReduced = U.T @ S
    else:
        SReduced = S

    print("SReduced", SReduced.shape)

    data_rows = SReduced.shape[0]
    data_cols = SReduced.shape[1]

    kratos_network.calculate_data_limits(SReduced)

    # Set the properties for the clusters
    num_clusters=1          # Number of different bases chosen
    num_cluster_col=5       # If I use num_cluster_col = num_variables result should be exact.
    num_encoding_var=num_cluster_col

    if config["load_model"]:
        autoencoder = tf.keras.models.load_model(kratos_network.model_name, custom_objects={"custom_loss":custom_loss, "custom_loss_cluster":custom_loss_cluster})
    else:
        autoencoder, autoencoder_err = kratos_network.define_network(SReduced, custom_loss_cluster, num_encoding_var)

    # Prepare the data
    SReducedNormalized = (SReduced - kratos_network.data_min) / (kratos_network.data_max - kratos_network.data_min)

    # Obtain Clusters
    print("=== Calculating Cluster Bases ===")
    with contextlib.redirect_stdout(None):
        CB, kmeans_object = clustering.calcualte_snapshots_with_columns(
            snapshot_matrix=SReducedNormalized,
            number_of_clusters=num_clusters,
            number_of_columns_in_the_basis=num_cluster_col,
            truncation_tolerance=1e-5
        )

    print(f"-> Generated {len(CB)} cluster bases with shapes:")
    for b in CB:
        print(f'-> {b=} has {CB[b].shape=}')
        SPri = CB[b] @ CB[b].T @ SReduced
        print(f"Check Norm for BASE (All    ): {np.linalg.norm(SPri-SReduced)/np.linalg.norm(SReduced)}")
        print(f"Check Norm for BASE (Col-100): {np.linalg.norm(SPri[5]-SReduced[5])/np.linalg.norm(SReduced[5])}")

    # Obtain the reduced representations (that would be our G's)
    q, s, c, b = {}, {}, {}, {}
    total_respresentations = 0

    for i in range(num_clusters):
        q[i] = CB[i] @ CB[i].T @ SReduced[:,kmeans_object.labels_==i]
        b[i] = CB[i] @ CB[i].T
        c[i] = np.array([i for _ in range(q[i].shape[1])])
        s[i] = SReduced[:,kmeans_object.labels_==i]

    print("Total number of representations:", total_respresentations)

    temp_size = num_cluster_col
    grad_size = num_cluster_col

    qc = np.empty(shape=(data_rows,0))
    sc = np.empty(shape=(data_rows,0))
    cc = np.empty(shape=(0))

    for i in range(num_clusters):
        qc = np.concatenate((qc, q[i]), axis=1)
        sc = np.concatenate((sc, s[i]), axis=1)
        cc = np.concatenate((cc, c[i]), axis=0)

    nqc = kratos_network.normalize_data(qc)
    nsc = kratos_network.normalize_data(sc)

    print(f"{cc.shape=}")

    if config["train_model"]:
        autoencoder.set_m_grad(b)
        autoencoder.set_g_weight(1)
        kratos_network.train_network(autoencoder, nsc, cc, len(data_inputs), 175)
        print("Model Initialized")

    if config["save_model"]:
        autoencoder.save(kratos_network.model_name)

    NSPredict = kratos_network.predict_snapshot(autoencoder, nsc)        # This is u', or f(q), or f(g(u)) 
    SPredict  = kratos_network.denormalize_data(NSPredict)

    print(f"{kratos_network.data_min=}")
    print(f"{kratos_network.data_max=}")
    # SPredict  = NSPredict

    # Compare a snap in the middle
    print(f'First elems of       {nsc.T[5][1:5]=}')
    print(f'First elems of {NSPredict.T[5][1:5]=}')

    # Check the gradients:
    grad_pred = np.array([NSPredict.T[5]])
    grad_expt = np.array([nsc.T[5]])

    grad_chk = check_gradient.GradientCheck(
        kratos_network.encoder_model, 
        kratos_network.get_gradients(
            kratos_network.encoder_model, 
            [kratos_network.encoder_model],
            grad_pred,
            grad_expt
        )
    )

    grad_chk.do_check(grad_expt)

    exit(0)

    print("U:", U.shape)

    if config["use_reduced"]:
        SP = U@(SPredict)
    else:
        SP = SPredict

    print("S  Shape:",  S.shape)
    print("SP Shape:", SP.shape)

    print("SP norm error", np.linalg.norm(SP-S)/np.linalg.norm(S))

    if config["print_results"]:
        current_model = KMP.Model()
        model_part = current_model.CreateModelPart("main_model_part")

        kratos_io.create_out_mdpa(model_part, "GidExampleSwaped")
        kratos_io.print_results_to_gid(model_part, S.T, SP.T)
