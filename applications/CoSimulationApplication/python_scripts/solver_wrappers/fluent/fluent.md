# SolverWrapperFluent

Documentation for all solver-wrappers for Fluent.



TODO: go through code, find out which stuff I should write down

TODO: go through setup of test case, write down procedure (what should be included, what shouldn't)


## Parameters

This section describes the Parameters in the JSON file.

JSON setting|type|version|description
---:|:---:|---|----
`working_directory`|string|all|absolute path to working directory or relative path w.r.t current directory
`case_file`|string|all|name of the case file; it must be present in the above defined working_directory
`dimensions`|int|all|2 for 2D and axisymmetric, 3 for 3D
`unsteady`|bool|all|true for transient FSI        
`delta_t`|double|all|fixed timestep size in flow solver
`timestep_start`|int|all|index of timestep to start transient FSI; 0 to start from case_file;
`thread_names`|list|all|list with Fluent names of the interface threads
`interface_input`|dict|all|keys are names of ModelParts for nodes, they must consist of an entry from thread_names + '_nodes'; value are (lists of) names of Variables
`interace_output`|dict|all|idem, but for faces
`cores`|int|all|number of processor cores to use (tested only on single node!)
`fluent_gui`|bool|all|true will run Fluent with graphical interface
`max_nodes_per_face`|int|all|used to get unique ID for faces, based on unique IDs of nodes; e.g. 4 for rectangular faces, 3 for triangular faces
`hybrid_initialization`|bool|all|true will run the hybrid initialization in Fluent before the first time-step; false requires that adequate reference values have been set in the case_file
`flow_iterations`|int|all|number of Fluent iterations per coupling iteration
`save_iterations`|int|all|number of timesteps between consecutive saves of the Fluent case and data files


`timestep_start` and `delta_t` are necessary parameters, but are usually defined in a higher CoSimulationObject; however, they can also be given directly as parameter of  solver wrapper; if they are defined both in higher object and in solver wrapper, then the former is used and warning printed (TODO!!)


## Version specific documentation

#### 2019R1 (19.3)

lorem ipsum

#### 2020R1 (20.1)

Not implemented yet.