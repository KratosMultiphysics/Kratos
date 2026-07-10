# UPDATE WARNING
Check:
- CMakeLists.txt -> delete examples
- source/CMakeLists.txt -> add cluster files
- gidpost_functions.h -> #include "gidpost_cluster_functions.h"
- gidpostFILES.c ->  case GiD_Cluster: && strElementType
- gidpostHDF5.c -> case GiD_Cluster: num_int=3; break;
- gidpost_types.h -> enum GiD_Cluster 
