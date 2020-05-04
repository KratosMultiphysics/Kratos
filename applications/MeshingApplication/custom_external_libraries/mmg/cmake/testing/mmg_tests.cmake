## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
##
##  mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the mmg distribution only if you accept them.
## =============================================================================

FOREACH(EXEC ${LISTEXEC_MMG})

  GET_FILENAME_COMPONENT ( SHRT_EXEC ${EXEC} NAME )

  ##############################################################################
  #####
  #####         Aniso test case
  #####
  ##############################################################################
  #####

  SET ( test_names
    mmg_TorusholesAni_${SHRT_EXEC}
    mmg_TorusholesAni_chocCyl_${SHRT_EXEC}
    )

  SET ( input_files
    ${MMG_CI_TESTS}/TorusholesAni/torusholes
    ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny
    )

  SET ( args
    "-v 5 -hgrad 1.15"
    "-v 5 -hgrad 1.15"
    )

  ADD_RUN_AGAIN_TESTS ( ${EXEC} "${test_names}" "${args}" "${input_files}" )

  ADD_TEST(NAME mmg_CubeVolAni_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5
  ${MMG_CI_TESTS}/CubeVolAni/cube
  -out ${CTEST_OUTPUT_DIR}/mmg_CubeVolAni_${SHRT_EXEC}-cube.o.meshb)

  ADD_TEST(NAME mmg_CubeVolAni2_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5
  ${MMG_CI_TESTS}/CubeVolAni2/cube
  -out ${CTEST_OUTPUT_DIR}/mmg_CubeVolAni2_${SHRT_EXEC}-cube.o.meshb)

  ADD_TEST(NAME mmg_SphereVolAni_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5
  ${MMG_CI_TESTS}/SphereVolAni/sphere
  -out ${CTEST_OUTPUT_DIR}/mmg_SphereVolAni_${SHRT_EXEC}-sphere.o.meshb)

  ##############################################################################
  #####
  #####         Check Memory Leak
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME mmg_LeakCheck_AbnormalEnd2_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd2/d)
  SET(passRegex "unable to scale mesh:")
  SET_PROPERTY(TEST mmg_LeakCheck_AbnormalEnd2_${SHRT_EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
 ADD_TEST(NAME mmg_LeakCheck_AbnormalEnd7_${SHRT_EXEC}
   COMMAND ${EXEC} -v 5
   ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/d
   -out ${CTEST_OUTPUT_DIR}/AbnormalEnd7/unwrittable.meshb)
 SET(passRegex "\\*\\* UNABLE TO OPEN.*")
 SET_PROPERTY(TEST mmg_LeakCheck_AbnormalEnd7_${SHRT_EXEC}
   PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
 #####
 ADD_TEST(NAME mmg_LeakCheck_AbnormalEnd8_${SHRT_EXEC}
   COMMAND ${EXEC} -v 5
   ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/d
   -out ${CTEST_OUTPUT_DIR}/AbnormalEnd8/unwrittable.meshb)
 SET(passRegex "\\*\\* UNABLE TO OPEN.*.sol")
 SET_PROPERTY(TEST mmg_LeakCheck_AbnormalEnd8_${SHRT_EXEC}
   PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
 #####
 #####
 ADD_TEST(NAME mmg_LeakCheck_args0_${SHRT_EXEC}
   COMMAND ${EXEC} -v 5
   ${MMG_CI_TESTS}/LeakCheck_args0/d
   ${CTEST_OUTPUT_DIR}/mmg_LeakCheck_args0_${SHRT_EXEC}-d.o)
 #####
 ADD_TEST(NAME mmg_LeakCheck_args1_${SHRT_EXEC}
   COMMAND ${EXEC} -v 5
   -in ${MMG_CI_TESTS}/LeakCheck_args1/d -sol
   ${MMG_CI_TESTS}/LeakCheck_args1/dsol.sol
   -out  ${CTEST_OUTPUT_DIR}/mmg_LeakCheck_args1_${SHRT_EXEC}-dout.meshb)

 ##############################################################################
 #####
 #####         Check Local parameters
 #####
 ##############################################################################
 #####
 ADD_TEST(NAME mmg_HausdLoc_2Spheres${SHRT_EXEC}
   COMMAND ${EXEC} -v 5 -hgrad 2
   ${MMG_CI_TESTS}/HausdLoc_2Spheres/2spheres
   ${CTEST_OUTPUT_DIR}/mmg_HausdLoc_2Spheres${SHRT_EXEC}-2spheres.o.meshb
   -hgrad 2
   )
 #####
 ADD_TEST(NAME mmg_hminmaxLoc_2Spheres${SHRT_EXEC}
   COMMAND ${EXEC} -v 5 -hgrad 2
   ${MMG_CI_TESTS}/HausdLoc_2Spheres/2spheres
   ${CTEST_OUTPUT_DIR}/mmg_HausdLoc_2Spheres${SHRT_EXEC}-2spheres.o.meshb
   -hgrad 2
   )


 ADD_TEST(NAME mmg_HausdLoc_2SpheresAni${SHRT_EXEC}
   COMMAND ${EXEC} -v 5 -hgrad 2 -A
   ${MMG_CI_TESTS}/HausdLoc_2Spheres/2spheres
   ${CTEST_OUTPUT_DIR}/mmg_HausdLoc_2SpheresAni${SHRT_EXEC}-2spheres.o.meshb
   -hgrad 2
   )
 #####
 ADD_TEST(NAME mmg_hminmaxLoc_2SpheresAni${SHRT_EXEC}
   COMMAND ${EXEC} -v 5 -hgrad 2 -A
   ${MMG_CI_TESTS}/HausdLoc_2Spheres/2spheres
   ${CTEST_OUTPUT_DIR}/mmg_HausdLoc_2SpheresAni${SHRT_EXEC}-2spheres.o.meshb
   -hgrad 2
   )


 ##############################################################################
 #####
 #####         Check Precision
 #####
 ##############################################################################
 #####
 ADD_TEST(NAME mmg_MeshVersionFormatted1_${SHRT_EXEC}
   COMMAND ${EXEC} -v 5
   -in ${MMG_CI_TESTS}/MeshVersionFormatted1/d
   -sol ${MMG_CI_TESTS}/MeshVersionFormatted1/dsol.sol
   ${CTEST_OUTPUT_DIR}/mmg_MeshVersionFormatted1_${SHRT_EXEC}-d.o
   )
 #####
 ADD_TEST(NAME mmg_MeshVersionFormatted2_${SHRT_EXEC}
   COMMAND ${EXEC} -v 5
   -in ${MMG_CI_TESTS}/MeshVersionFormatted2/d
   -sol ${MMG_CI_TESTS}/MeshVersionFormatted2/dsol.sol
   ${CTEST_OUTPUT_DIR}/mmg_MeshVersionFormatted2_${SHRT_EXEC}-d.o
   )

###############################################################################
#####
#####         Options
#####
###############################################################################
ADD_TEST(NAME mmg_hsizOption_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5 -hsiz 0.1
  ${MMG_CI_TESTS}/Cube/cube
  -out ${CTEST_OUTPUT_DIR}/mmg_hsiz_${SHRT_EXEC}.o.meshb)

# hsiz Aniso
ADD_TEST(NAME mmg_hsizAni_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5 -hsiz 0.1 -sol 2 -A
  ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny
  -out ${CTEST_OUTPUT_DIR}/mmg_hsizAni_${SHRT_EXEC}.o.meshb)

ADD_TEST(NAME mmg_hsizHmax_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5 -hsiz 0.1 -hmax 0.05
  ${MMG_CI_TESTS}/Cube/cube
  -out ${CTEST_OUTPUT_DIR}/mmg_hsizHmax_${SHRT_EXEC}.o.meshb)
SET(passRegex "Mismatched options")
SET_PROPERTY(TEST mmg_hsizHmax_${SHRT_EXEC}
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

ADD_TEST(NAME mmg_hsizHmin_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5 -hsiz 0.1 -hmin 0.2
  ${MMG_CI_TESTS}/Cube/cube
  -out ${CTEST_OUTPUT_DIR}/mmg_hsizHmin_${SHRT_EXEC}.o.meshb)
SET_PROPERTY(TEST mmg_hsizHmin_${SHRT_EXEC}
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

# Required entities
ADD_TEST(NAME mmg_MultiDom_Ellipse_ReqEntities_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5 -hausd 0.002
  ${MMG_CI_TESTS}/MultiDom_Ellipse_ReqEntities/c.d
  -out ${CTEST_OUTPUT_DIR}/mmg_MultiDom_Ellipse_ReqEntities_${SHRT_EXEC}.o.meshb)

ADD_TEST(NAME mmg_MultiDom_Cube_ReqEntities_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5 -hsiz 0.02
  ${MMG_CI_TESTS}/MultiDom_Cube_ReqEntities/c
  -out ${CTEST_OUTPUT_DIR}/mmg_MultiDom_Cube_ReqEntities_${SHRT_EXEC}.o.meshb)

ADD_TEST(NAME mmg_MultiDom_Ellipse_ReqEntitiesAni_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5 -hausd 0.002 -A
  ${MMG_CI_TESTS}/MultiDom_Ellipse_ReqEntities/c.d
  -out ${CTEST_OUTPUT_DIR}/mmg_MultiDom_Ellipse_ReqEntitiesAni_${SHRT_EXEC}.o.meshb)


# -A
ADD_TEST(NAME mmg_CommandLineAni_${SHRT_EXEC}
  COMMAND ${EXEC} -v 5 -hausd 0.005 -sol 2 -A
  ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny
  -out ${CTEST_OUTPUT_DIR}/mmg_CommandLineAni_${SHRT_EXEC}.o.meshb)

  ##############################################################################
  #####
  #####         Various test cases
  #####
  ##############################################################################
  #####

  # Lot of reference edges, ridges, corners and singularities
  ADD_TEST(NAME mmg_SurfEdges_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 -hausd 0.1 -A
    ${MMG_CI_TESTS}/SurfEdges_house/housebad.meshb
    -out ${CTEST_OUTPUT_DIR}/mmg_SurfEdgesAni_${SHRT_EXEC}.o.meshb)

  ADD_TEST(NAME mmg_SurfEdgesAni_${SHRT_EXEC}
    COMMAND ${EXEC} -v 5 -hausd 0.1
    ${MMG_CI_TESTS}/SurfEdges_house/housebad.meshb
    -out ${CTEST_OUTPUT_DIR}/mmg_SurfEdges_${SHRT_EXEC}.o.meshb)


ENDFOREACH(EXEC)
