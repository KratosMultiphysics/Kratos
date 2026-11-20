"""Processes that run a trained model inside the Kratos solution loop.

Each one follows the same contract: gather Kratos data into tensors, call the
model, and write the prediction back onto Kratos entities - de-normalized
according to the model card when one is present
(``deployment.model_registry``).

They differ only in the *shape* of the data they gather:

=================================== ==========================================
Module                              Gathers
=================================== ==========================================
``inference_process``               flat per-node/element feature vectors
``graph_inference_process``         the element-edge graph of the ``ModelPart``
``grid_inference_process``          a regular voxel grid sampled from the mesh
``point_cloud_inference_process``   the mesh nodes as an unordered cloud
``sequence_inference_process``      a rolling window of grid states
``time_series_inference_process``   a rolling window of nodal states
``particle_inference_process``      a proximity graph over a particle cloud
``diffusion_inference_process``     a conditioning field, sampled as an ensemble
``domino_inference_process``        a preprocessed DoMINO datapipe case
``onnx_inference_process``          as ``inference_process``, via ONNX Runtime
``triton_inference_process``        as ``inference_process``, via a Triton server
``nim_inference_process``           a triangle skin, POSTed to a PhysicsNeMo NIM
``superresolution_process``         a coarse grid, writing a fine one
``rom_surrogate_process``           case parameters, writing ``u = phi q``
``hybrid_initialization_process``   nothing - it seeds the solver before it runs
``pinn_solve_process``              nothing - it *is* the solve
=================================== ==========================================
"""
