"""Training loops, datasets and training-time schemes.

======================== ======================================================
Module                   Provides
======================== ======================================================
``training_utils``       ``TrainModel``/``SaveTrainedModel`` - the
                         ``Parameters``-driven loop, callbacks, warm restarts,
                         OOD calibration and FSDP2-safe checkpoint writing
    A ``performance`` block adds mixed precision, physicsnemo's static
    capture, the Muon optimizer, learning-rate schedules, profiling and
    resumable checkpoints - all off by default.
``torch_dataset``        dataset and datapipe factories over exported files
                         (``.npz``, ``.pmsh``, DoMINO/Transolver, augmentation,
                         mixing)
``streaming_dataset``    ``LiveSampleQueue`` and ``StreamingDataset`` - training
                         straight out of a running solve, no file round trip
``temporal_training``    window datasets and ``TrainAutoregressive`` (BPTT
                         through a self-fed rollout)
``diffusion_utils``      ``TrainDiffusionModel``/``TrainCorrDiffPair`` and the
                         EDM sampler wrappers
``domino_finetune``      predictor-corrector and LoRA adaptation of a pretrained
                         DoMINO checkpoint
``rom_temporal``         reduced-trajectory datasets and temporal-attention
                         training in ROM space
``rollout_utils``        ``EvaluateRollout`` - multi-step error growth of a
                         trained time-series surrogate
``sobolev_training``     ``MakeSensitivityLossTerm`` - derivative-informed
                         training against Kratos's adjoint dJ/dX
``globe_training``       ``TrainGlobe`` - the dict-in/dict-out loop GLOBE needs,
                         which the tensor-batch ``TrainModel`` cannot express
``aerojepa_pretraining`` self-supervised geometry pretraining, with the JEPA
                         objective and EMA target encoder upstream declares but
                         does not ship
======================== ======================================================

Nothing here defines a Kratos ``Process``; the process that fills the streaming
queue is ``processes.export.streaming_dataset_export_process``.
"""
