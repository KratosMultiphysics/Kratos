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

  ##############################################################################
  #####
  #####         Aniso test case
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME CubeVolAni_${EXEC}
  COMMAND ${EXEC} -v 5
  ${MMG_CI_TESTS}/CubeVolAni/cube
  -out ${MMG_CI_TESTS}/CubeVolAni/cube.o.meshb)

  ADD_TEST(NAME CubeVolAni2_${EXEC}
  COMMAND ${EXEC} -v 5
  ${MMG_CI_TESTS}/CubeVolAni2/cube
  -out ${MMG_CI_TESTS}/CubeVolAni2/cube.o.meshb)

  ADD_TEST(NAME SphereVolAni_${EXEC}
  COMMAND ${EXEC} -v 5
  ${MMG_CI_TESTS}/SphereVolAni/sphere
  -out ${MMG_CI_TESTS}/SphereVolAni/sphere.o.meshb)

  ADD_TEST(NAME TorusholesAni_${EXEC}
  COMMAND ${EXEC} -v 5 -hgrad 1.15
  ${MMG_CI_TESTS}/TorusholesAni/torusholes
  -out ${MMG_CI_TESTS}/TorusholesAni/torusholes.o.meshb)

 ADD_TEST(NAME TorusholesAni_chocCyl_${EXEC}
  COMMAND ${EXEC} -v 5 -hgrad 1.15
  ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny
  -out ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny.o.meshb)
  ##############################################################################
  #####
  #####         Check Memory Leak
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd2_${EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd2/d)
  SET(passRegex "## Unable to scale mesh.")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd2_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd7_${EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/d
    -out ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/unwrittable.meshb)
  SET(passRegex "\\*\\* UNABLE TO OPEN.*")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd7_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  ADD_TEST(NAME LeakCheck_AbnormalEnd8_${EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/d
    -out ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/unwrittable.meshb)
  SET(passRegex "\\*\\* UNABLE TO OPEN.*.sol")
  SET_PROPERTY(TEST LeakCheck_AbnormalEnd8_${EXEC}
    PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
  #####
  #####
  ADD_TEST(NAME LeakCheck_args0_${EXEC}
    COMMAND ${EXEC} -v 5
    ${MMG_CI_TESTS}/LeakCheck_args0/d)
  #####
  ADD_TEST(NAME LeakCheck_args1_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${MMG_CI_TESTS}/LeakCheck_args1/d -sol
    ${MMG_CI_TESTS}/LeakCheck_args1/dsol.sol
    -out ${MMG_CI_TESTS}/LeakCheck_args1/dout.meshb)

  ##############################################################################
  #####
  #####         Check Local parameters
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME HausdLoc_2Spheres${EXEC}
    COMMAND ${EXEC} -v 5 -hgrad 2
    ${MMG_CI_TESTS}/HausdLoc_2Spheres/2spheres
    ${MMG_CI_TESTS}/HausdLoc_2Spheres/2spheres.o.meshb
    -hgrad 2
    )
  #####
  ADD_TEST(NAME hminmaxLoc_2Spheres${EXEC}
    COMMAND ${EXEC} -v 5 -hgrad 2
    ${MMG_CI_TESTS}/HausdLoc_2Spheres/2spheres
    ${MMG_CI_TESTS}/HausdLoc_2Spheres/2spheres.o.meshb
    -hgrad 2
    )


  ##############################################################################
  #####
  #####         Check Precision
  #####
  ##############################################################################
  #####
  ADD_TEST(NAME MeshVersionFormatted1_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${MMG_CI_TESTS}/MeshVersionFormatted1/d
    -sol ${MMG_CI_TESTS}/MeshVersionFormatted1/dsol.sol)
  #####
  ADD_TEST(NAME MeshVersionFormatted2_${EXEC}
    COMMAND ${EXEC} -v 5
    -in ${MMG_CI_TESTS}/MeshVersionFormatted2/d
    -sol ${MMG_CI_TESTS}/MeshVersionFormatted2/dsol.sol)

ENDFOREACH(EXEC)
