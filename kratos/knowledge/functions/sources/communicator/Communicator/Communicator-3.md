---
type: C++ Method
title: Communicator
resource: sources/communicator.cpp#L58-L68
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Communicator::Communicator(Communicator const& rOther) : mNumberOfColors(rOther.mNumberOfColors) , mNeighbourIndices(rOther.mNeighbourIndices) , mpLocalMesh(MeshType::Pointer(rOther.mpLocalMesh)) , mpGhostMesh(MeshType::Pointer(rOther.mpGhostMesh)) , mpInterfaceMesh(MeshType::Pointer(rOther.mpInterfaceMesh)) , mLocalMeshes(rOther.mLocalMeshes) , mGhostMeshes(rOther.mGhostMeshes) , mInterfaceMeshes(rOther.mInterfaceMeshes) , mrDataCommunicator(rOther.mrDataCommunicator)`