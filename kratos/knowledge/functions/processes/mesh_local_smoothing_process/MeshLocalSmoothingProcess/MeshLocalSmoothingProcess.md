---
type: C++ Method
title: MeshLocalSmoothingProcess
resource: processes/mesh_local_smoothing_process.cpp#L33-L42
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`MeshLocalSmoothingProcess::MeshLocalSmoothingProcess( ModelPart &rModelPart, double AptQuality, std::size_t IterationsNumber, const Flags& rBoundaryFlag) :mrModelPart(rModelPart), mMaxIterationsNumber(IterationsNumber), mAptQuality(AptQuality), mNumberOfLowQualityElements(0) ,mMeshMinQuality(0.00), mMeshQualityNorm(0.00), mrBoundaryFlag(rBoundaryFlag)`