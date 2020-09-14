
# Cluster Generator
Available for __Windows__ and __Linux__.

Based on the copyrighted __SphereTree Toolkit__ found [here](http://isg.cs.tcd.ie/spheretree/).

Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following paragraphs appear in all copies.

The __SphereTree Toolkit__ authors may be contacted at the following e-mail addresses:
- Gareth_Bradshaw@yahoo.co.uk
- isg@cs.tcd.ie



# Set up
Once downloaded and added as problemtype for GiD, the required executables must be copied to the exec folder.
- For Windows the precompiled executables can be found [here](http://isg.cs.tcd.ie/spheretree/downloads/spheretree-1.0-win32.zip)
- For Linux, the precompiled executables are already located in the exec folder.

# Use
Once the geometry has been generated and meshed, the user can choose between step-by-step process or automatically generated the final cluster file.

The step-by-step process go as follows:
- Adjust the meshing parameters and mesh the geometry.
- Execute step one in the "Sphere Cluster Creation" menu in order to generate the OBJ and MSH files
- Define the options for the spheretree algorithms or use the default options (recommended)
- Select generate SPH (see examples for execution time references)
- Select generate CLU
- Select visualize cluster to draw the generated cluster over the mesh. The cluster can also be visualized in principal axis.

# Examples
Both examples are created on a Intel i7 laptop.

- Jar, 6000 tetrahedra
branch: 20, time: 2 min
branch: 200, time: 14 min

<span>
<img align="center" src="https://github.com/KratosMultiphysics/Kratos/blob/dem-cluster_creation/applications/DEMApplication/custom_problemtype/DEMClusters.gid/images/jar1.bmp" width="288">
</span>
<br>

<span>
<img align="center" src="https://github.com/KratosMultiphysics/Kratos/blob/dem-cluster_creation/applications/DEMApplication/custom_problemtype/DEMClusters.gid/images/jar1-1.bmp" width="288">

</span>
<br>


- Candy, 268 tetrahedra
branch: 20, time: 1 min
branch: 200, time: 7 min

<span>
<img align="center" src="https://github.com/KratosMultiphysics/Kratos/blob/dem-cluster_creation/applications/DEMApplication/custom_problemtype/DEMClusters.gid/images/candy1.bmp" width="288">
</span>
<br>

<span>
<img align="center" src="https://github.com/KratosMultiphysics/Kratos/blob/dem-cluster_creation/applications/DEMApplication/custom_problemtype/DEMClusters.gid/images/candy1-1.bmp" width="288">
</span>
<br>



- Wrench, 1614 tetrahedra

<span>
<img align="center" src="https://github.com/KratosMultiphysics/Kratos/blob/dem-cluster_creation/applications/DEMApplication/custom_problemtype/DEMClusters.gid/images/wrench1.bmp" width="288">
</span>
<br>

<span>
<img align="center" src="https://github.com/KratosMultiphysics/Kratos/blob/dem-cluster_creation/applications/DEMApplication/custom_problemtype/DEMClusters.gid/images/center4.png" width="288">
</span>
<br>





# Recommendations and troubleshooting
In order to avoid typical issues when generating the cluster.
- The default values for the SPH algorithm are usually stable. Other configurations may require a more detailed mesh or calibration of other SPH options.
- Do not generate the geometry from an existing mesh (avoid bad faces definition)
- If using a copy of an existing geometry via save as, save and reload the problem before continuing.
- The mesh is not automatically generated. The user must specify the meshing parameters and generate the mesh prior to create the SPH file or the cluster.
- On Windows, GID may throw an error when finishing generating the SPH file, although the file is correctly generated in the problem folder.
- For other path related problems, try reloading the problem.


