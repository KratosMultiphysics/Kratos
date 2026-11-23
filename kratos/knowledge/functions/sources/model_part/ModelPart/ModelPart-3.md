---
type: C++ Method
title: ModelPart
resource: sources/model_part.cpp#L124-L143
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`ModelPart::ModelPart(std::string const& NewName, IndexType NewBufferSize,VariablesList::Pointer pVariablesList, Model& rOwnerModel) : DataValueContainer() , Flags() , mBufferSize(NewBufferSize) , mpProcessInfo(new ProcessInfo()) , mpVariablesList(pVariablesList) , mpCommunicator(new Communicator) , mpParentModelPart(NULL) , mSubModelParts() , mrModel(rOwnerModel)`