"""Taking a trained checkpoint to production.

================================== ============================================
Module                             Provides
================================== ============================================
``model_registry``                 loading checkpoints (physicsnemo ``.mdlus``,
                                   TorchScript, NGC), model cards, output
                                   de-normalization, ``torch.compile``
``onnx_utils``                     the ONNX Runtime session factory and its
                                   device parsing; ``ExportOnnxModel`` (which
                                   writes the portable artifact) lives with the
                                   other save paths in ``training_utils``
``triton_export``                  a Triton Inference Server model repository
                                   plus its generated ``config.pbtxt``
``nim_client``                     the documented HTTP contract of PhysicsNeMo
                                   NIM microservices - multipart STL in,
                                   ``.npz`` out, stdlib urllib + numpy only
``usd_export``                     time-sampled OpenUSD stages (digital twins)
                                   from plain arrays, via optional ``usd-core``
``cosim_surrogate_solver_wrapper`` a trained model as a first-class
                                   CoSimulationApplication *solver*
``surrogate_response_function``    a trained model as a Kratos *response
                                   function* - J and dJ/dX for optimization
``uncertainty_utils``              MC dropout, checkpoint ensembles and GP heads
                                   - calibrated error bars on any prediction
``ood_guard_utils``                out-of-distribution guardrails, calibrated at
                                   training time and checked per inference
``geometry_guard_utils``           the same question asked of the *shape*: a
                                   density model over surface descriptors,
                                   flagging geometries outside the training
                                   family before a prediction is trusted
================================== ============================================

The processes that *run* these artifacts are in ``processes.inference``.
"""
