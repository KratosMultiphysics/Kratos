---
type: C++ Method
title: Communicator
resource: sources/communicator.cpp#L30-L41
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Communicator::Communicator() : mNumberOfColors(1) , mpLocalMesh(MeshType::Pointer(new MeshType)) , mpGhostMesh(MeshType::Pointer(new MeshType)) , mpInterfaceMesh(MeshType::Pointer(new MeshType)) , mrDataCommunicator(ParallelEnvironment::GetDataCommunicator("Serial"))`