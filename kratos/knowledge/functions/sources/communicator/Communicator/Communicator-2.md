---
type: C++ Method
title: Communicator
resource: sources/communicator.cpp#L44-L55
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Communicator::Communicator(const DataCommunicator& rDataCommunicator) : mNumberOfColors(1) , mpLocalMesh(MeshType::Pointer(new MeshType)) , mpGhostMesh(MeshType::Pointer(new MeshType)) , mpInterfaceMesh(MeshType::Pointer(new MeshType)) , mrDataCommunicator(rDataCommunicator)`