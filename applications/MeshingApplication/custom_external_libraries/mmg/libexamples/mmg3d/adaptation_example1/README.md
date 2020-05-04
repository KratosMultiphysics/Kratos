#Example to move from mmg3d4 library to mmg3d5 library

   * the **_main.cold_** correspond to a mmg3d4 library call;
   * the **_main.c_** is the same code but with a call to the new mmg3d library (5.x.x release).

It is strongly advised to use the mmg3d5 API functions instead of hard setting your mesh as it was required by the mmg3d4 library.

**Remarks:**
  * You will find the same example (with additionnal boundary triangles) in the **_example0_** directory;
  * please, refer to the **_example0_** directory to have an example of _"clean"_ call of the new **mmg3d** library and to see how to build the example executable.