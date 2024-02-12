---
title: PointerCommunicator
keywords: 
tags: [How-to-use-PointerCommunicator.md]
sidebar: kratos_for_developers
summary: 
---

Kratos provides a data proxy mechanism, based on the use of "GlobalPointers" to simplify communications involving non local data. To understand the mechanism, let's consider a local loop to take the average of neighbours nodes for each node. In serial it looks something like:

```cpp
for(auto& node : model_part.Nodes()) {
    double Tavg = 0;
    const auto& neighbours = node.GetValue(NODAL_NEIGHBOURS);
    for(auto& neighb : neighbours)
         Tavg += neighb.FastGetSolutionStepValue(TEMPERATURE);
    node.SetValue(AVERAGE_TEMPERATURE, Tavg/(neigbours.size());
}
```
{: data-lang="C++"}

this loop is obviously not possible in MPI. (we do not have all the neighbours on a given process unless we add a lot of ghost nodes)

Let's now consider we use global pointers instead of weak_pointers for the neighbours. **IN THE SERIAL CASE** the previous loop would look very similar, something like

```cpp
for(auto& node : model_part.Nodes()) {
    double Tavg = 0;
    const auto& global_pointers_to_neighbours = node.GetValue(GP_NODAL_NEIGHBOURS);
    for(auto& gp : global_pointers_to_neighbours)
         Tavg += gp->FastGetSolutionStepValue(TEMPERATURE);
    node.SetValue(AVERAGE_TEMPERATURE, Tavg/(global_pointers_to_neighbours.size());
}
```
{: data-lang="C++"}

Since global pointers are pointers that are only valid in the memory space of the owner processor this works ok in SMP (on a single process), but still fails in MPI, since the global neighbours may belong to different MPI processes.

In order to retrofit this problem, the former code should be modified to:

```cpp
// First of all collect all the neighbours into a list
GlobalPointersVector< Node > gp_list;

for(auto& node : model_part.Nodes())
    for(auto& gp : global_pointers_to_neighbours)
        gp_list.push_back(gp);

// Now create the pointer communicator
GlobalPointerCommunicator pointer_comm(DefaultDataCommunicator, gp_list); // Here we prepare the communication

// Define a function to access to the data
result_proxy = pointer_comm.Apply<double>(
            [](GlobalPointer<Node >& gp) -> double
            {
                  return gp->FastGetSolutionStepValue(TEMPERATURE);
            }
      );  // Here all communications happen !!

for(auto& node : model_part.Nodes()) {
    double Tavg = 0;
    const auto& global_pointers_to_neighbours = node.GetValue(GP_NODAL_NEIGHBOURS);
    for(auto& gp : global_pointers_to_neighbours)
         Tavg += result_proxy.Get(gp); // Here we get the result of the application of the function onto gp. This is so BOTH for LOCAL GPs and for REMOTE GPs
    node.SetValue(AVERAGE_TEMPERATURE, Tavg/(global_pointers_to_neighbours.size());
}
```
{: data-lang="C++"}

The code is of course a little more verbose, nevertheless all the communication is hidden in the `Apply` function (and in the constructor), and the use if still more or less similar to the original serial code.

Aside of this, the code can be made to be still efficient in serial. The point here is that the function `get` of the `result_proxy` is actually executing the function directly if the globalpointer is local. This means that no remote data at all is needed (and hence computed) if the communicator returns `IsDistributed() -> false`.
This also means that additional storage can be avoided completely in the serial case.
In the practice this means that we can further optimize the code to have no overhead in the serial case as follows:

```cpp
// Now create the pointer communicator --> NOTHING IS ACTUALLY DONE IN THE SERIAL CASE!
GlobalPointerCommunicator pointer_comm(DefaultDataCommunicator, 
    [&](){ //NOTE THAT THIS FUNCTOR WILL ONLY BE EXECUTED IN MPI MODE, it has no overhead in serial!!
        GlobalPointersVector< Node > gp_list;        
        for(auto& node : model_part.Nodes())
            for(auto& gp : global_pointers_to_neighbours)
                gp_list.push_back(gp);
    }); 

// Define a function to access to the data ----> NO COMMUNICATION HAPPENS IN THE SERIAL CASE
// we only store the access function within result_proxy and we will use it on every "Get" call
result_proxy = pointer_comm.Apply<double>(
            [](GlobalPointer<Node >& gp) -> double
            {
                  return gp->FastGetSolutionStepValue(TEMPERATURE);
            }
      ); // Here all communications happen !!

for(auto& node : model_part.Nodes()) {
    double Tavg = 0;
    const auto& global_pointers_to_neighbours = node.GetValue(GP_NODAL_NEIGHBOURS);
    for(auto& gp : global_pointers_to_neighbours)
         Tavg += result_proxy.Get(gp); // Here we get the result of the application of the function onto gp. This is so BOTH for LOCAL GPs and for REMOTE GPs
    node.SetValue(AVERAGE_TEMPERATURE, Tavg/(global_pointers_to_neighbours.size());
}
```
{: data-lang="C++"}

With such modification the code is almost as fast in serial as the original code, the only overhead being that a function call by function pointer is done on every global pointer.

