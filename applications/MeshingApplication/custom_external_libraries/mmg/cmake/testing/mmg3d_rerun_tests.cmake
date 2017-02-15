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

###############################################################################
#####
#####         Continuous Integration
#####
###############################################################################

# simple test: must already pass
ADD_TEST(NAME SimpleCube2
  COMMAND ${EXECUT_MMG3D} -v 6 -d
  ${MMG3D_CI_TESTS}/Cube/cube.o
  -out ${MMG3D_CI_TESTS}/Cube/cube.o.o.meshb)


###############################################################################
#####
#####         Check Results
#####
###############################################################################
#####

# Check what happend when we refine an isotropic cube of size h with a constant
# metric (h, h/2, h/4, h/8 and h/16)
#---First with hmin=hmax
ADD_TEST(NAME CubeIso_h_hminMax2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeIso_h_hminMax/CubeIso0.1.o -hmax 0.1 -hmin 0.1 -v 5
  -out ${MMG3D_CI_TESTS}/CubeIso_h_hminMax/CubeIso0.1.o.o.meshb)
ADD_TEST(NAME CubeIso_0.5h_hminMax2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeIso_0.5h_hminMax/CubeIso0.1.o -hmax 0.05 -hmin 0.05 -v 5
  -out ${MMG3D_CI_TESTS}/CubeIso_0.5h_hminMax/CubeIso0.1.o.o.meshb)
ADD_TEST(NAME CubeIso_0.25h_hminMax2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeIso_0.25h_hminMax/CubeIso0.1.o -hmax 0.025 -hmin 0.025 -v 5
  -out ${MMG3D_CI_TESTS}/CubeIso_0.25h_hminMax/CubeIso0.1.o.o.meshb)
#ADD_TEST(NAME CubeIso_0.125h_hminMax
#  COMMAND ${EXECUT_MMG3D}
#  ${MMG3D_CI_TESTS}/CubeIso_0.125h_hminMax/CubeIso0.1 -hmax 0.0125 -hmin 0.0125 -v 5)

#---Second with sol file
ADD_TEST(NAME CubeIso_h_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeIso_h_met/CubeIso0.1.o -v 5
  -out ${MMG3D_CI_TESTS}/CubeIso_h_met/CubeIso0.1.o.o.meshb)
ADD_TEST(NAME CubeIso_0.5h_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeIso_0.5h_met/CubeIso0.1.o -v 5
  -out ${MMG3D_CI_TESTS}/CubeIso_0.5h_met/CubeIso0.1.o.o.meshb)
ADD_TEST(NAME CubeIso_0.25h_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeIso_0.25h_met/CubeIso0.1.o -v 5
  -out ${MMG3D_CI_TESTS}/CubeIso_0.25h_met/CubeIso0.1.o.o.meshb)
ADD_TEST(NAME CubeIso_0.125h_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeIso_0.125h_met/CubeIso0.1.o -v 5
  -out ${MMG3D_CI_TESTS}/CubeIso_0.125h_met/CubeIso0.1.o.o.meshb)

#####

# Check what happend when we refine a sphere of size h with a constant metric
# (h, h/2, h/4 and h/8)
##---First with hmin=hmax
#ADD_TEST(NAME SphereIso_h_hminMax
#  COMMAND ${EXECUT_MMG3D}
#  ${MMG3D_CI_TESTS}/SphereIso_h_hminMax/SphereIso0.5
#  -hmax 0.5 -hmin 0.5 -hausd 0.1 -v 5)
#ADD_TEST(NAME SphereIso_0.5h_hminMax
#  COMMAND ${EXECUT_MMG3D}
#  ${MMG3D_CI_TESTS}/SphereIso_0.5h_hminMax/SphereIso0.5
#  -hmax 0.25 -hmin 0.25 -hausd 0.1 -v 5)
#ADD_TEST(NAME SphereIso_0.25h_hminMax
#  COMMAND ${EXECUT_MMG3D}
#  ${MMG3D_CI_TESTS}/SphereIso_0.25h_hminMax/SphereIso0.5
#  -hmax 0.125 -hmin 0.125 -hausd 0.1 -v 5)
#ADD_TEST(NAME SphereIso_0.125h_hminMax
#  COMMAND ${EXECUT_MMG3D}
#  ${MMG3D_CI_TESTS}/SphereIso_0.125h_hminMax/SphereIso0.5
#  -hmax 0.0625 -hmin 0.0625 -hausd 0.1 -v 5)
#---Second with sol file
ADD_TEST(NAME SphereIso_h_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_h_met/SphereIso0.5.o -hausd 0.1 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_h_met/SphereIso0.5.o.o.meshb)
ADD_TEST(NAME SphereIso_0.5h_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_0.5h_met/SphereIso0.5.o -hausd 0.1 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_0.5h_met/SphereIso0.5.o.o.meshb)
ADD_TEST(NAME SphereIso_0.25h_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_0.25h_met/SphereIso0.5.o -hausd 0.1 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_0.25h_met/SphereIso0.5.o.o.meshb)
ADD_TEST(NAME SphereIso_0.125h_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_0.125h_met/SphereIso0.5.o -hausd 0.1 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_0.125h_met/SphereIso0.5.o.o.meshb)

# Check what happend when we unrefine a sphere of size smallh with a constant metric
# (2*smallh, 4*smallh and 8*smallh)
ADD_TEST(NAME SphereIso_2smallh_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_2smallh_met/SphereIso0.0625.o -hausd 0.1 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_2smallh_met/SphereIso0.0625.o.o.meshb)
ADD_TEST(NAME SphereIso_4smallh_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_4smallh_met/SphereIso0.0625.o -hausd 0.1 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_4smallh_met/SphereIso0.0625.o.o.meshb)
ADD_TEST(NAME SphereIso_8smallh_met2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_8smallh_met/SphereIso0.0625.o -hausd 0.1 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_8smallh_met/SphereIso0.0625.o.o.meshb)

# Check what happend when we use hausdorff number to refine the skin and a big hgrad
# to have an inside of the initial size (0.5)
ADD_TEST(NAME SphereIso_h_hausd0.001_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_h_hausd0.001/SphereIso0.5.o -hausd 0.001 -hgrad 500 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_h_hausd0.001/SphereIso0.5.o.o.meshb)
ADD_TEST(NAME SphereIso_h_hausd0.005_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/SphereIso_h_hausd0.005/SphereIso0.5.o -hausd 0.005 -hgrad 100 -v 5
  -out ${MMG3D_CI_TESTS}/SphereIso_h_hausd0.005/SphereIso0.5.o.o.meshb)

# Check what happend when we refine a cube whose skin has already the good size
ADD_TEST(NAME CubeSkin0.05_Inside0.4_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeSkin0.05_Inside0.4/CubeSkin0.05.o -v 5
  ${MMG3D_CI_TESTS}/CubeSkin0.05_Inside0.4/CubeSkin0.05.o.o.meshb)
ADD_TEST(NAME CubeSkin0.1_Inside0.4_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeSkin0.1_Inside0.4/CubeSkin0.1.o -v 5
  ${MMG3D_CI_TESTS}/CubeSkin0.1_Inside0.4/CubeSkin0.1.o.o.meshb)
ADD_TEST(NAME CubeSkin0.2_Inside0.4_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeSkin0.2_Inside0.4/CubeSkin0.2.o -v 5
  ${MMG3D_CI_TESTS}/CubeSkin0.2_Inside0.4/CubeSkin0.2.o.o.meshb)
ADD_TEST(NAME CubeSkin0.0125_Inside0.125_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.125/CubeSkin0.125.o -v 5
  -out ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.125/CubeSkin0.125.o.o.meshb)
ADD_TEST(NAME CubeSkin0.0125_Inside0.25_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.25/CubeSkin0.25.o -v 5
  -out ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.25/CubeSkin0.25.o.o.meshb)
ADD_TEST(NAME CubeSkin0.0125_Inside0.5_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.5/CubeSkin0.5.o -v 5
  ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.5/CubeSkin0.5.o.o.meshb)


# Check results on various meshes
# First: Meshes that we want unrefined
ADD_TEST(NAME Various_unref_Linkrods_met0.2_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/Various_unref_Linkrods_met0.2/linkrods.o -hausd 0.1 -v 5
  ${MMG3D_CI_TESTS}/Various_unref_Linkrods_met0.2/linkrods.o.o.meshb)
ADD_TEST(NAME Various_unref_Linkrods_met0.2_hausd0.01_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/Various_unref_Linkrods_met0.2_hausd0.01/linkrods.o
  -hausd 0.01 -v 5
  ${MMG3D_CI_TESTS}/Various_unref_Linkrods_met0.2_hausd0.01/linkrods.o.o.meshb)



# Second: Meshes that we want refined
ADD_TEST(NAME Various_ref_Linkrods_met0.05_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05/linkrods.o -hausd 0.1 -v 5
  ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05/linkrods.o.o.meshb)
ADD_TEST(NAME Various_ref_Linkrods_met0.05_hausd0.01_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05_hausd0.01/linkrods.o
  -hausd 0.01 -v 5
  ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05_hausd0.01/linkrods.o.o.meshb)
ADD_TEST(NAME Various_ref_Linkrods_met0.05_hausd0.001_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05_hausd0.001/linkrods.o
  -hausd 0.001 -v 5
  ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05_hausd0.001/linkrods.o.o.meshb)

# Third: We refine some parts and unrefined others
ADD_TEST(NAME Various_refunref_Santa_met0.05_hausd0.001_ar90_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/Various_refunref_Santa_met0.05_hausd0.001_ar90/santa.o
  -hausd 0.001 -ar 90 -v 5
  ${MMG3D_CI_TESTS}/Various_refunref_Santa_met0.05_hausd0.001_ar90/santa.o.o.meshb)
ADD_TEST(NAME Various_refunref_Santa_met0.05_hausd0.0001_ar90_2
  COMMAND ${EXECUT_MMG3D}
  ${MMG3D_CI_TESTS}/Various_refunref_Santa_met0.05_hausd0.0001_ar90/santa.o
  -hausd 0.0001 -ar 90 -v 5
  ${MMG3D_CI_TESTS}/Various_refunref_Santa_met0.05_hausd0.0001_ar90/santa.o.o.meshb)

# 4: Refinment on a solution
IF ( PATTERN )
  ADD_TEST(NAME Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2_2
    COMMAND ${EXECUT_MMG3D}
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6.o
    -sol ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6.o.sol
    -hausd 0.1 -ar 60 -hgrad 1 -v 5
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6.o.o.meshb)
  ADD_TEST(NAME Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3_2
    COMMAND ${EXECUT_MMG3D}
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6.o
    -sol
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6.o.sol
    -hausd 0.1 -ar 60 -hgrad 1 -v 5
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6.o.o.meshb)
ELSE ()
  ADD_TEST(NAME Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2_2
    COMMAND ${EXECUT_MMG3D}
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6.o
    -sol ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6.o.sol
    -hausd 0.1 -ar 60 -hgrad 1
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6.o.o.meshb)
  ADD_TEST(NAME Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3_2
    COMMAND ${EXECUT_MMG3D}
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6.o -v 5
    -sol
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6.o.sol
    -hausd 0.1 -ar 60 -hgrad 1
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6.o.o.meshb)
ENDIF ()

# Test multi-domain remeshing
ADD_TEST(NAME MultiDom_Cube2
  COMMAND ${EXECUT_MMG3D} -v 6 -hmax 0.02 ${MMG3D_CI_TESTS}/MultiDom_Cube/c.o
  -out ${MMG3D_CI_TESTS}/MultiDom_Cube/c.o.o.meshb)

ADD_TEST(NAME MultiDom_Ellipse_2
  COMMAND ${EXECUT_MMG3D} -v 6 -hausd 0.0003 ${MMG3D_CI_TESTS}/MultiDom_Ellipse/c.d.o
  -out ${MMG3D_CI_TESTS}/MultiDom_Ellipse/c.d.o.o.meshb)

# Non-manifold test case
ADD_TEST(NAME NM_Cube2
  COMMAND ${EXECUT_MMG3D} -v 6 -d -hmax 0.05 ${MMG3D_CI_TESTS}/NM_Cube/nm.o
  -out ${MMG3D_CI_TESTS}/NM_Cube/nm.o.o.meshb)


ADD_TEST(NAME TorusholesAni_${EXECUT_MMG3D}_2
 COMMAND ${EXECUT_MMG3D} -v 5 -hgrad 1.15
 ${MMG_CI_TESTS}/TorusholesAni/torusholes.o.meshb
 -out ${MMG_CI_TESTS}/TorusholesAni/torusholes.o.o.meshb)

SET_TESTS_PROPERTIES(TorusholesAni_${EXECUT_MMG3D}_2 PROPERTIES DEPENDS TorusholesAni_${EXECUT_MMG3D})

ADD_TEST(NAME TorusholesAni_chocCyl_${EXECUT_MMG3D}_2
 COMMAND ${EXECUT_MMG3D} -v 5 -hgrad 1.15
 ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny.o.meshb
 -out ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny.o.o.meshb)

SET_TESTS_PROPERTIES(TorusholesAni_chocCyl_${EXECUT_MMG3D}_2 PROPERTIES DEPENDS TorusholesAni_chocCyl_${EXECUT_MMG3D})
