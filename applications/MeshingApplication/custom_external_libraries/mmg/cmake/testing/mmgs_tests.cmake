## =============================================================================
##  This file is part of the Mmg software package for the tetrahedral
##  mesh modification.
##**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
##
##  Mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  Mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with Mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the Mmg distribution only if you accept them.
## =============================================================================

###############################################################################
#####
#####         Continuous Integration
#####
###############################################################################

# Simple tests: must already pass
ADD_TEST(NAME SimpleTeapot
  COMMAND ${EXECUT_MMGS}
  ${MMGS_CI_TESTS}/Teapot/teapot
  -out ${MMGS_CI_TESTS}/Teapot/teapot.d.meshb)

ADD_TEST(NAME CubeAni
  COMMAND ${EXECUT_MMGS}
  ${MMGS_CI_TESTS}/CubeAni/cube
  -out ${MMGS_CI_TESTS}/CubeAni/cube.d.meshb)

ADD_TEST(NAME SphereAni
  COMMAND ${EXECUT_MMGS}
  ${MMGS_CI_TESTS}/SphereAni/sphere
  -out ${MMGS_CI_TESTS}/SphereAni/sphere.d.meshb)


###############################################################################
#####
#####         Input/Output
#####
###############################################################################

# Binary gmsh
ADD_TEST(NAME binary_gmsh_s
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/GmshInout/cube.mshb)

# Ascii gmsh
ADD_TEST(NAME ascii_gmsh_s
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/GmshInout/cube.msh)


###############################################################################
#####
#####         Check Memory Leaks
#####
###############################################################################


###############################################################################
#####
#####         Manifold cases
#####
###############################################################################
ADD_TEST(NAME Rhino_M
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/Rhino_M/rhino -hausd 1
  -out ${MMGS_CI_TESTS}/Rhino_M/rhino.d.meshb)

###############################################################################
#####
#####         Non manifold cases
#####
###############################################################################
ADD_TEST(NAME Cow_NM_hausd10
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/Cow_NM/cow -hausd 10
  -out ${MMGS_CI_TESTS}/Cow_NM/cow.d.meshb)

###############################################################################
#####
#####         Test results
#####
###############################################################################
# Test the Ls option
ADD_TEST(NAME OptLs_teapot
  COMMAND ${EXECUT_MMGS} -v 5 -ls
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot.simple.o.meshb)

ADD_TEST(NAME OptLs_teapot_keepRef
  COMMAND ${EXECUT_MMGS} -v 5 -ls -keep-ref
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot.keep-ref.o.meshb)

ADD_TEST(NAME OptLs_teapot_0.5_keepRef
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.5 -keep-ref
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot.0.5.keep-ref.o.meshb)

ADD_TEST(NAME OptLs_teapot2
  COMMAND ${EXECUT_MMGS} -v 5 -ls -nr
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot.o.meshb)


###############################################################################
#####
#####         Detected Bugs
#####
###############################################################################
ADD_TEST(NAME Car_NM
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/Car_NM/car
  -out ${MMGS_CI_TESTS}/Car_NM/car.d.meshb)

ADD_TEST(NAME Cow_NM_hausd20
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/Cow_NM/cow -hausd 20
  -out ${MMGS_CI_TESTS}/Cow_NM/cow.d.meshb)
