---
title: PointerMapCommunicator
keywords: 
tags: [How-to-use-PointerMapCommunicator.md]
sidebar: kratos_for_developers
summary: 
---

Kratos provides a data proxy mechanism, based on the use of "GlobalPointers" to simplify communications involving non local data. In order to retrieve information on local and remote nodes, one can follow this wiki page [How to use PointerCommunicator](How-to-use-PointerCommunicator) to get more information.

The `GlobalPointerMapCommunicator` is designed to be used in a situation where, it is required to set an *entity* (could be `Node`, `Element` or `Condition`, etc...) value. This **entity* can be in a local process (which is always the case in serial runs) or local/remote processes (which is the case for distributed runs). So `GlobalPointerMapCommunicator` takes care of the communication burden therefore user does not need to bother whether simulation is run on serial or distributed.

To explain the usage, lets take the following example of doing assembly in a serial run where an element value (i.e. `TEMPERATURE`) needs to be distributed among its nodal neighbours.
```cpp
// first TEMPERATURE variable is set to zero.
VariableUtils.SetNonHistoricalVariableToZero(TEMPERATURE, model_part.Nodes);

for (auto& r_element : model_part.Elements()) {
    auto& r_geometry = r_element.GetGeometry();
    const auto& temperature = r_element.GetValue(TEMPERATURE);
    for (auto& r_node : r_geometry) {
       const auto& neighbours = node.GetValue(NODAL_NEIGHBOURS);
       for (auto& gp : neighbours) {
           gp->GetValue(TEMPERATURE) += temperature;
       }
    }
}
```
This works in serial run, but in a distributed run, `GlobalPointer`s given by `NODAL_NEIGHBOURS` will belong to either local or a remote processes. Therefore, this may require setting a value in a remote process. Above will give `Segmentaion Fault` error if a given `GlobalPointer` does not belong to the local process. These `GlobalPointer`s may also not belong to `GhostMesh` of the local process. So if the user plans to use `Assembly` methods in the `Communicator`, that will also fail to achieve the expected outcome if used without modifying existing `GhostMesh` (In there you will need to add all the remote nodes present in each of the node's `NODAL_NEIGHBOURS` list to `GhostMesh` which will hinder the performance of normal operations in distributed run). 

So in order retrofit problems in the above mentioned task, and to achieve same outcome in serial or distributed run; One can use the following block of code with `GlobalPointerMapCommunicator`
```cpp
// first TEMPERATURE variable is set to zero.
VariableUtils.SetNonHistoricalVariableToZero(TEMPERATURE, model_part.Nodes);

// create the global pointer map communicator
// will not do anything in serial case.
GlobalPointerMapCommunicator<Node, double> pointer_map_comm(r_default_comm);

// creates the proxy method for updating
// this proxy is called upon only with gps in their owning processes. No communication is done in here.
// this proxy needs to be always thread safe for parallel loops of the entitiy being passed.
// In this case the entity is of type Node
const auto& apply_proxy = pointer_map_comm.GetApplyProxy([](Node& rNode, const double& NewValue) {
    rNode.GetValue(TEMPERATURE) += NewValue;
});

// do the calculation
for (auto& r_element : model_part.Elements()) {
    auto& r_geometry = r_element.GetGeometry();
    const auto& temperature = r_element.GetValue(TEMPERATURE);
    for (auto& r_node : r_geometry) {
       const auto& neighbours = node.GetValue(NODAL_NEIGHBOURS);
       for (auto& gp : neighbours) {
           // this will directly call the apply_proxy methods lambda function
           // in the serial case or in distributed when gp is local, 
           // otherwise the value will be stored in a non-local-gp data containers. No communication is done in here.
           apply_proxy.Assign(gp, temperature);
       }
    }
}

// do mpi communication
apply_proxy.SendAndApplyRemotely();

// need to synchronize to update ghost mesh. Does nothing in the serial run
model_part.GetCommunicator().SynchronizeNonHistoricalVariable(TEMPERATURE);
```
In the above code block, If the communicator returns `IsDistributed() -> false` then followings will be the behaviour:
1. When `pointer_map_comm` is constructed, it will do nothing (except for reference assignments)
2. When `GetApplyProxy` is called with a lambda function, it will not do anything (only reference assignment)
3. When `apply_proxy.Assign` is called it will directly call the lambda function given in the `apply_proxy` (No checks performed, simply doing what is mentioned in the lambda, hence least performance impact). This can be called even in OMP parallel regions
4. `apply_proxy.SendAndApplyRemotely()` will be a blank call.

In the case where communicator returns `IsDistributed() -> true` then it will have the following behaviour
1. When `pointer_map_comm` is constructed, it will do nothing (except for reference assignments)
2. When `GetApplyProxy` is called with a lambda function, it will reserve space (the size of them will be `number_of_threads`) to store `GlobalPointers` and `DataValues` which are non-local. (Cannot be called in OMP parallel regions)
3. When `apply_proxy.Assign` is called it will directly call the lambda function given in the `apply_proxy` if the given `GlobalPointer` is local, if not it will store `GlobalPointer` and `DataValue` for future communication. This can be called even in OMP parallel regions
4. `apply_proxy.SendAndApplyRemotely()` will do communication for non-local `GlobalPointer` data values, and apply the proxy in their owning processes (This includes creating communication schedules, and the sending and receiving data to/from remote processes). This should not be called in OMP parallel regions.

**In the MPI run, it is essential to do `Synchronization` for the variables which are modified through `apply_proxy` lambda to have consistent values in the `GhostMesh` because `SendAndApplyRemotely()` will not modify any of the entities in the `GhostMesh`.**