---
type: C++ Method
title: MeasureMeshQualityProcess
resource: processes/measure_mesh_quality_process.cpp#L28-L39
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`MeasureMeshQualityProcess::MeasureMeshQualityProcess(ModelPart& rModelPart, std::size_t Dimension) : mrModelPart(rModelPart) , mDimension(Dimension) , mNumberOfInvertedElements(0) , mNumberOfSlivers(0) , mMinArea(std::numeric_limits<double>::max()) , mMaxArea(std::numeric_limits<double>::lowest()) , mMinH(std::numeric_limits<double>::max()) , mMinAngle(std::numeric_limits<double>::max())`