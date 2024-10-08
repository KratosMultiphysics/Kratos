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

##############################################################################
#####
#####         Tests that may be run twice
#####
##############################################################################

SET ( test_names
  # Simple test: must already pass
  mmg3d_SimpleCube_fast
  # MultiDomain
  mmg3d_MultiDom_Ellipse_fast
  # Non-manifold test case
  mmg3d_NM_Cube_fast
  mmg3d_NM_Complex_fast
  # test case with non-manifold, ridges, ref edges and a curve surface
  mmg3d_NM_cone_fast
  # mmg3d_NM_cone_ani_fast #Fail because at second run a tetra we have a tet with 4 ridge vertices
  )

SET ( input_files
  ${MMG3D_CI_TESTS}/Cube/cube
  ### Multidom
  ${MMG3D_CI_TESTS}/MultiDom_Ellipse/c.d
   ### non-manifold
  ${MMG3D_CI_TESTS}/NM_Cube/nm
  ${MMG3D_CI_TESTS}/NM_Complex/nm4
  ${MMG3D_CI_TESTS}/cone-nm.mesh
  #${MMG3D_CI_TESTS}/cone-nm.mesh
  )

SET ( args
  "-v 5"
  ### MultiDomain
  "-v 5 -hausd 0.002"
  ### non-manifold
  "-v 5 -hmax 0.1"
  "-v 5"
  "-v 5"
  #"-v 5 -A"
  )

IF ( LONG_TESTS )
  SET ( test_names ${test_names}
    # Check what happend when we refine an isotropic cube of size h with a
    # constant metric (h, h/2, h/4, h/8 and h/16) ---First with hmin=hmax
    mmg3d_CubeIso_h_hminMax
    mmg3d_CubeIso_0.5h_hminMax
    mmg3d_CubeIso_0.25h_hminMax
    #---Second with sol file
    mmg3d_CubeIso_h_met
    mmg3d_CubeIso_0.5h_met
    mmg3d_CubeIso_0.25h_met
    mmg3d_CubeIso_0.125h_met
    mmg3d_CubeAniIso_0.125h_met
    #####
    mmg3d_SphereIso_h_met
    mmg3d_SphereIso_0.5h_met
    mmg3d_SphereIso_0.25h_met
    mmg3d_SphereIso_0.125h_met
    mmg3d_SphereIso_0.020_met
    # mmg3d_SphereIso_0.020-0.015_met # not enough mem on windows 4G
    mmg3d_SphereAni_0.02
    # Check what happend when we unrefine a sphere of size smallh with a
    # constant metric (2*smallh, 4*smallh and 8*smallh)
    mmg3d_SphereIso_2smallh_met
    mmg3d_SphereIso_4smallh_met
    mmg3d_SphereIso_8smallh_met
    # Check what happend when we use hausdorff number to refine the skin and a
    # big hgrad to have an inside of the initial size (0.5)
    mmg3d_SphereIso_h_hausd0.001
    mmg3d_SphereIso_h_hausd0.005
    # Check what happend when we refine a cube whose skin has already the good size
    mmg3d_CubeSkin0.05_Inside0.4
    mmg3d_CubeSkin0.1_Inside0.4
    mmg3d_CubeSkin0.2_Inside0.4
    mmg3d_CubeSkin0.0125_Inside0.125
    # mmg3d_CubeSkin0.0125_Inside0.25 # too long on OSX
    # mmg3d_CubeSkin0.0125_Inside0.5 # too long on all machine
    # Check results on various meshes
    # First: Meshes that we want unrefined
    mmg3d_Various_unref_Linkrods_met0.2
    mmg3d_Various_unref_Linkrods_met0.2_hausd0.01
    # Second: Meshes that we want refined
    mmg3d_Various_ref_Linkrods_met0.05
    mmg3d_Various_ref_Linkrods_met0.05_hausd0.01
    mmg3d_Various_ref_Linkrods_met0.05_hausd0.001
    # Third: We refine some parts and unrefined others
    mmg3d_Various_refunref_Santa_met0.05_hausd0.001_ar90
    mmg3d_Various_refunref_Santa_met0.05_hausd0.0001_ar90
    # 5: MultiDomain
    mmg3d_MultiDom_Cube
    mmg3d_MultiDom_Ellipse
    # Non-manifold test case
    mmg3d_NM_Cube
    mmg3d_NM_Complex
    )

  SET ( input_files  ${input_files}
    ### Cube
    ${MMG3D_CI_TESTS}/CubeIso_h_hminMax/CubeIso0.1
    ${MMG3D_CI_TESTS}/CubeIso_0.5h_hminMax/CubeIso0.1
    ${MMG3D_CI_TESTS}/CubeIso_0.25h_hminMax/CubeIso0.1
    ###
    ${MMG3D_CI_TESTS}/CubeIso_h_met/CubeIso0.1
    ${MMG3D_CI_TESTS}/CubeIso_0.5h_met/CubeIso0.1
    ${MMG3D_CI_TESTS}/CubeIso_0.25h_met/CubeIso0.1
    ${MMG3D_CI_TESTS}/CubeIso_0.125h_met/CubeIso0.1
    ${MMG3D_CI_TESTS}/CubeAniIso_0.125h_met/CubeIso0.1
    ### Sphere
    ${MMG3D_CI_TESTS}/SphereIso_h_met/SphereIso0.5
    ${MMG3D_CI_TESTS}/SphereIso_0.5h_met/SphereIso0.5
    ${MMG3D_CI_TESTS}/SphereIso_0.25h_met/SphereIso0.5
    ${MMG3D_CI_TESTS}/SphereIso_0.125h_met/SphereIso0.5
    ${MMG3D_CI_TESTS}/SphereIso_0.020_met/SphereIso0.5
    # ${MMG3D_CI_TESTS}/SphereIso_0.020-0.015_met/SphereIso0.020
    ${MMG3D_CI_TESTS}/SphereAni_0.02/sphere
    ###
    ${MMG3D_CI_TESTS}/SphereIso_2smallh_met/SphereIso0.0625
    ${MMG3D_CI_TESTS}/SphereIso_4smallh_met/SphereIso0.0625
    ${MMG3D_CI_TESTS}/SphereIso_8smallh_met/SphereIso0.0625
    ###
    ${MMG3D_CI_TESTS}/SphereIso_h_hausd0.001/SphereIso0.5
    ${MMG3D_CI_TESTS}/SphereIso_h_hausd0.005/SphereIso0.5
    ### CubeSkin
    ${MMG3D_CI_TESTS}/CubeSkin0.05_Inside0.4/CubeSkin0.05
    ${MMG3D_CI_TESTS}/CubeSkin0.1_Inside0.4/CubeSkin0.1
    ${MMG3D_CI_TESTS}/CubeSkin0.2_Inside0.4/CubeSkin0.2
    ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.125/CubeSkin0.125
    # ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.25/CubeSkin0.25
    # ${MMG3D_CI_TESTS}/CubeSkin0.0125_Inside0.5/CubeSkin0.5
    ### Linkrods
    ${MMG3D_CI_TESTS}/Various_unref_Linkrods_met0.2/linkrods
    ${MMG3D_CI_TESTS}/Various_unref_Linkrods_met0.2_hausd0.01/linkrods
    ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05/linkrods
    ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05_hausd0.01/linkrods
    ${MMG3D_CI_TESTS}/Various_ref_Linkrods_met0.05_hausd0.001/linkrods
    ### Santa
    ${MMG3D_CI_TESTS}/Various_refunref_Santa_met0.05_hausd0.001_ar90/santa
    ${MMG3D_CI_TESTS}/Various_refunref_Santa_met0.05_hausd0.0001_ar90/santa
    ### MultiDomain
    ${MMG3D_CI_TESTS}/MultiDom_Cube/c
    ${MMG3D_CI_TESTS}/MultiDom_Ellipse/c.d
    ${MMG3D_CI_TESTS}/NM_Cube/nm
    ${MMG3D_CI_TESTS}/NM_Complex/nm4
    )

  SET ( args ${args}
    ### Cube
    "-v 5 -hmax 0.1 -hmin 0.1"
    "-v 5 -hmax 0.05 -hmin 0.05"
    "-v 5 -hmax 0.025 -hmin 0.025"
    ###
    "-v 5"
    "-v 5"
    "-v 5"
    "-v 5"
    "-v 5"
    ### Sphere
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    # "-v 5 -hausd 0.1"
    "-v 5"
    ###
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.1"
    ###
    "-v 5 -hausd 0.001 -hgrad -1"
    "-v 5 -hausd 0.005 -hgrad -1"
    ### CubeSkin
    "-v 5"
    "-v 5"
    "-v 5"
    "-v 5"
    # "-v 5"
    # "-v 5"
    ### Linkrods
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.01"
    "-v 5 -hausd 0.1"
    "-v 5 -hausd 0.01"
    "-v 5 -hausd 0.001"
    ### Santa
    "-v 5 -hausd 0.001 -ar 90"
    "-v 5 -hausd 0.0001 -ar 90"
    ### MultiDomain
    "-v 5 -hmax 0.02"
    "-v 5 -hausd 0.0003"
    "-v 5 -hmax 0.05"
    "-v 5"
    )

ENDIF ( )

ADD_RUN_AGAIN_TESTS ( ${EXECUT_MMG3D} "${test_names}" "${args}" "${input_files}" )

IF ( LONG_TESTS )
  ### M6
  SET ( test_name
    # 4: Refinment on a solution
    mmg3d_Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2
    )
  SET ( input_file
    ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/M6
    )

  ADD_TEST(NAME ${test_name}
    COMMAND ${EXECUT_MMG3D}
    ### M6
    ${input_file}
    -v 5 -sol ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.001_hmin0.001_hmax2/metM6.sol -hausd 0.1 -ar 60
    -out ${CTEST_OUTPUT_DIR}/${test_name}-out.o.meshb )

  SET_TESTS_PROPERTIES ( ${test_name}
    PROPERTIES FIXTURES_SETUP ${test_name} )

  IF ( RUN_AGAIN )
    ADD_TEST(NAME ${test_name}_2
      COMMAND ${EXECUT_MMG3D}
      -v 5 -hausd 0.1 -ar 60 -hgrad -1
      ${CTEST_OUTPUT_DIR}/${test_name}-out.o.meshb
      -out ${CTEST_OUTPUT_DIR}/${test_name}_2-out.o.meshb
      )

    SET_TESTS_PROPERTIES ( ${test_name}_2
      PROPERTIES FIXTURES_REQUIRED ${test_name} )

  ENDIF ( RUN_AGAIN )

  SET ( test_name
      # 4: Refinment on a solution
      mmg3d_Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3
      )
  SET ( input_file
      ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/M6
      )

  ADD_TEST(NAME ${test_name}
    COMMAND ${EXECUT_MMG3D}
    ${input_file}
    ### M6
    ${input_file}
    -v 5 -sol ${MMG3D_CI_TESTS}/Various_adpsol_hgrad1_M6Mach_Eps0.0005_hmin0.0001_hmax3/metM6.sol -hausd 0.1 -ar 60
    -out ${CTEST_OUTPUT_DIR}/${test_name}-out.o.meshb )

  SET_TESTS_PROPERTIES ( ${test_name}
    PROPERTIES FIXTURES_SETUP ${test_name} )

    IF ( RUN_AGAIN )
      ADD_TEST(NAME ${test_name}_2
        COMMAND ${EXECUT_MMG3D}
        -v 5 -hausd 0.1 -ar 60 -hgrad -1
        ${CTEST_OUTPUT_DIR}/${test_name}-out.o.meshb
        -out ${CTEST_OUTPUT_DIR}/${test_name}_2-out.o.meshb
        )

    SET_TESTS_PROPERTIES ( ${test_name}_2
      PROPERTIES FIXTURES_REQUIRED ${test_name} )

  ENDIF ( )
ENDIF ( LONG_TESTS )


###############################################################################
#####
#####         Input/Output
#####
###############################################################################

# Binary gmsh
ADD_TEST(NAME mmg3d_binary_gmsh_3d
  COMMAND ${EXECUT_MMG3D} -v 5
  ${MMG3D_CI_TESTS}/GmshInout/cube.mshb
  ${CTEST_OUTPUT_DIR}/mmg3d_binary_gmsh_3d-cube.o
  )

# Ascii gmsh
ADD_TEST(NAME mmg3d_ascii_gmsh_3d
  COMMAND ${EXECUT_MMG3D} -v 5
  ${MMG3D_CI_TESTS}/GmshInout/cube.msh
  ${CTEST_OUTPUT_DIR}/mmg3d_ascii_gmsh_3d-cube.o
)

# Tetgen
# Default Tetgen behaviour saves only boundary tria (resp. edges) in
# .face (resp. .edge) file.
ADD_TEST ( NAME mmg3d_cube-tetgen
  COMMAND ${EXECUT_MMG3D} -v 5
  ${MMG3D_CI_TESTS}/Cube/cube
  ${CTEST_OUTPUT_DIR}/mmg3d_cube-tetgen.o.node
 )

##############################################################################
#####
#####         Check Memory Leak
#####
##############################################################################
#####
ADD_TEST(NAME mmg3d_LeakCheck_AbnormalEnd3
  COMMAND ${EXECUT_MMG3D} -v 5
  ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd3/d -sol
  ${MMG3D_CI_TESTS}/LeakCheck_AbnormalEnd3/dsol.sol -ls
  -out ${CTEST_OUTPUT_DIR}/mmg3d_LeakCheck_AbnormalEnd3-d.o.meshb)
SET(passRegex "## ERROR: UNABLE TO LOAD SOLUTION")
SET_PROPERTY(TEST mmg3d_LeakCheck_AbnormalEnd3
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
#####
ADD_TEST(NAME mmg3d_LeakCheck_optLevelSet
  COMMAND ${EXECUT_MMG3D} -v 5  -ls -hgrad 1.5
  ${MMG3D_CI_TESTS}/LeakCheck_optLevelSet/rect03d
  -out ${CTEST_OUTPUT_DIR}/mmg3d_LeakCheck_optLevelSet-rect03d.o.meshb)

##############################################################################
#####
#####         Check Options
#####
##############################################################################
#####
ADD_TEST(NAME mmg3d_memOption
  COMMAND ${EXECUT_MMG3D} -v 5 -m 100
  ${MMG3D_CI_TESTS}/Cube/cube
  -out ${CTEST_OUTPUT_DIR}/mmg3d_memOption.o.meshb)

ADD_TEST(NAME mmg3d_hsizAndNosurfOption
  COMMAND ${EXECUT_MMG3D} -v 5 -hsiz 0.1 -nosurf
  ${MMG3D_CI_TESTS}/Cube/cube
  -out ${CTEST_OUTPUT_DIR}/mmg3d_hsizNosurf.o.meshb)

ADD_TEST(NAME mmg3d_hsizAndNosurfAni
  COMMAND ${EXECUT_MMG3D} -v 5 -hsiz 0.1 -nosurf -A
  ${MMG3D_CI_TESTS}/Cube/cube
  -out ${CTEST_OUTPUT_DIR}/mmg3d_hsizNosurfAni.o.meshb)

ADD_TEST(NAME mmg3d_val
  COMMAND ${EXECUT_MMG3D} -v 5 -val
  ${MMG3D_CI_TESTS}/Cube/cube
  ${CTEST_OUTPUT_DIR}/mmg3d_cube-val.o.meshb
  )
SET_PROPERTY(TEST mmg3d_val
  PROPERTY WILL_FAIL TRUE)

ADD_TEST(NAME mmg3d_locParamCrea
  COMMAND ${EXECUT_MMG3D} -v 5 -default
  ${MMG3D_CI_TESTS}/LocParamsCrea/step.0)

SET_TESTS_PROPERTIES ( mmg3d_locParamCrea
  PROPERTIES FIXTURES_SETUP mmg3d_locParamCrea )
ADD_TEST(NAME mmg3d_locParamClean
  COMMAND ${CMAKE_COMMAND} -E remove -f
  ${MMG3D_CI_TESTS}/LocParamsCrea/step.mmg3d)
SET_TESTS_PROPERTIES ( mmg3d_locParamClean
  PROPERTIES FIXTURES_REQUIRED mmg3d_locParamCrea )

# default hybrid
ADD_TEST(NAME mmg3d_hybrid_3d
  COMMAND ${EXECUT_MMG3D} -v 5
  ${MMG3D_CI_TESTS}/Hybrid/prism.mesh
  ${CTEST_OUTPUT_DIR}/mmg3d_hybrid_3d-default.msh)

# nsd + hybrid
ADD_TEST(NAME mmg3d_hybrid-nsd1
  COMMAND ${EXECUT_MMG3D} -v 5 -nsd 1
  ${MMG3D_CI_TESTS}/Hybrid/prism.mesh
  ${CTEST_OUTPUT_DIR}/mmg3d_hybrid-nsd.mesh)

###############################################################################
#####
#####         Check Boundaries
#####
###############################################################################
#####
ADD_TEST(NAME mmg3d_ChkBdry_optls_test4
  COMMAND ${EXECUT_MMG3D} -v 5  -ls -hgrad 1.5
  -in ${MMG3D_CI_TESTS}/ChkBdry_optls_test4/test4
  -sol ${MMG3D_CI_TESTS}/ChkBdry_optls_test4/test4.sol
  -out ${CTEST_OUTPUT_DIR}/mmg3d_ChkBdry_optls_test4-test4.o.meshb)
#####
ADD_TEST(NAME mmg3d_ChkBdry_optls_temp
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -hmin 5 -hmax 6
  -nr -hausd 0.5 -hgrad 1.2
  -in ${MMG3D_CI_TESTS}/ChkBdry_optls_temp/temp
  -sol ${MMG3D_CI_TESTS}/ChkBdry_optls_temp/temp.sol
  -out ${CTEST_OUTPUT_DIR}/mmg3d_ChkBdry_optls_temp-temp.o.meshb)
####
ADD_TEST(NAME mmg3d_ChkBdry_optls_temp2
  COMMAND ${EXECUT_MMG3D} -v 5  -ls -hmin 5 -hmax 6
  -nr -hausd 0.5 -hgrad 1.2
  -in ${MMG3D_CI_TESTS}/ChkBdry_optls_temp/temp
  -sol ${MMG3D_CI_TESTS}/ChkBdry_optls_temp/temp.sol
  -out ${CTEST_OUTPUT_DIR}/mmg3d_ChkBdry_optls_temp2-temp.o.meshb)
#####
ADD_TEST(NAME mmg3d_ChkBdry_cube
  COMMAND ${EXECUT_MMG3D} -v 5
  ${MMG3D_CI_TESTS}/ChkBdry_cube/cube
  ${CTEST_OUTPUT_DIR}/mmg3d_ChkBdry_cube-cube.o
  )
#####
ADD_TEST(NAME mmg3d_ChkBdry_multidomCube
  COMMAND ${EXECUT_MMG3D} -v 5 -hmax 0.1
  ${MMG3D_CI_TESTS}/ChkBdry_multidomCube/c
  ${CTEST_OUTPUT_DIR}/mmg3d_ChkBdry_multidomCube-cube.o
  )
#####
ADD_TEST(NAME mmg3d_ChkBdry_multidomCube2
  COMMAND ${EXECUT_MMG3D} -v 5 -hmax 0.1
  ${MMG3D_CI_TESTS}/ChkBdry_multidomCube2/c
  ${CTEST_OUTPUT_DIR}/mmg3d_ChkBdry_multidomCube2-c.o
  )
#####
ADD_TEST(NAME mmg3d_ChkBdry_multidomCube3
  COMMAND ${EXECUT_MMG3D} -v 5 -hmax 0.1
  ${MMG3D_CI_TESTS}/ChkBdry_multidomCube3/c
  ${CTEST_OUTPUT_DIR}/mmg3d_ChkBdry_multidomCube2-cube.o
  )

ADD_TEST(NAME mmg3d_opnbdy_unref_peninsula
  COMMAND ${EXECUT_MMG3D} -v 5 -opnbdy
  -in ${MMG3D_CI_TESTS}/OpnBdy_peninsula/peninsula
  -out ${CTEST_OUTPUT_DIR}/mmg3d_OpnBdy_peninsula.o.meshb)

ADD_TEST(NAME mmg3d_opnbdy_ls_peninsula
  COMMAND ${EXECUT_MMG3D} -v 5 -opnbdy -ls
  -in ${MMG3D_CI_TESTS}/OpnBdy_peninsula/peninsula
  -sol  ${MMG3D_CI_TESTS}/OpnBdy_peninsula/ls.sol
  -out ${CTEST_OUTPUT_DIR}/mmg3d_OpnBdy_ls_peninsula.o.meshb)

ADD_TEST(NAME mmg3d_opnbdy_lssurf-nofile_peninsula
  COMMAND ${EXECUT_MMG3D} -v 5 -opnbdy -lssurf 0.6 -nr -hgrad 1.5 -hausd 0.02
  -in ${MMG3D_CI_TESTS}/OpnBdy_peninsula/peninsula
  -sol  ${MMG3D_CI_TESTS}/OpnBdy_peninsula/ls.sol
  -out ${CTEST_OUTPUT_DIR}/mmg3d_OpnBdy_lssurf_peninsula.o.meshb)

# ls + nsd
ADD_TEST(NAME mmg3d_opnbdy_ls_peninsula-nsd3
  COMMAND ${EXECUT_MMG3D} -v 5 -opnbdy -ls -nsd 3
  -in ${MMG3D_CI_TESTS}/OpnBdy_peninsula/peninsula
  -sol  ${MMG3D_CI_TESTS}/OpnBdy_peninsula/ls.sol
  -out ${CTEST_OUTPUT_DIR}/mmg3d_OpnBdy_ls_peninsula-nsd3.o.meshb)

ADD_TEST(NAME mmg3d_opnbdy_ref_peninsula
  COMMAND ${EXECUT_MMG3D} -v 5 -hmax 0.06 -opnbdy
  -in ${MMG3D_CI_TESTS}/OpnBdy_peninsula/peninsula
  -out ${CTEST_OUTPUT_DIR}/mmg3d_OpnBdy_peninsula.o.meshb)

ADD_TEST(NAME mmg3d_opnbdy_unref_island
  COMMAND ${EXECUT_MMG3D} -v 5 -opnbdy
  -in ${MMG3D_CI_TESTS}/OpnBdy_island/island
  -out ${CTEST_OUTPUT_DIR}/mmg3d_OpnBdy_island.o.meshb)

ADD_TEST(NAME mmg3d_opnbdy_ref_island
  COMMAND ${EXECUT_MMG3D} -v 5 -hmax 0.06 -opnbdy
  -in ${MMG3D_CI_TESTS}/OpnBdy_island/island
  -out ${CTEST_OUTPUT_DIR}/mmg3d_OpnBdy_island.o.meshb)

###############################################################################
#####
#####         Check Lagrangian motion option
#####
###############################################################################
#####
IF ( ELAS_FOUND AND NOT USE_ELAS MATCHES OFF )
  ADD_TEST(NAME mmg3d_LagMotion0_tinyBoxt
    COMMAND ${EXECUT_MMG3D} -v 5  -lag 0
    -in ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt
    -sol ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt.sol
    -out ${CTEST_OUTPUT_DIR}/mmg3d_LagMotion0_tinyBoxt-tinyBoxt.o.meshb
    )
  ADD_TEST(NAME mmg3d_LagMotion1_tinyBoxt
    COMMAND ${EXECUT_MMG3D} -v 5  -lag 1
    -in ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt
    -sol ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt.sol
    -out ${CTEST_OUTPUT_DIR}/mmg3d_LagMotion1_tinyBoxt-tinyBoxt.o.meshb
    )
  ADD_TEST(NAME mmg3d_LagMotion2_tinyBoxt
    COMMAND ${EXECUT_MMG3D} -v 5  -lag 2
    -in ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt
    -sol ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt.sol
    -out ${CTEST_OUTPUT_DIR}/mmg3d_LagMotion2_tinyBoxt-tinyBoxt.o.meshb
    )
  # nsd
  ADD_TEST(NAME mmg3d_LagMotion2_tinyBoxt-nsd3
    COMMAND ${EXECUT_MMG3D} -v 5  -lag 2 -nsd 3
    -in ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt
    -sol ${MMG3D_CI_TESTS}/LagMotion1_tinyBoxt/tinyBoxt.sol
    -out ${CTEST_OUTPUT_DIR}/mmg3d_LagMotion2_tinyBoxt-nsd3.o.meshb
    )

ENDIF()

##############################################################################
#####
#####         Check Local parameters at tetra
#####
##############################################################################
#####
ADD_TEST(NAME mmg3d_TetLoc_Ellipse
  COMMAND ${EXECUT_MMG3D} -v 5 -hgrad -1
  ${MMG3D_CI_TESTS}/TetLoc_Ellipse/c
  ${CTEST_OUTPUT_DIR}/mmg3d_TetLoc_Ellipse-c.o.meshb
  -hgrad 2
  )

##############################################################################
#####
#####         Check optim + aniso option
#####
##############################################################################
#####
ADD_TEST(NAME mmg3d_OptimAni_Sphere
  COMMAND ${EXECUT_MMG3D} -v 5 -optim -A
  ${MMG3D_CI_TESTS}/SphereIso_h_met/SphereIso0.5.meshb -sol 2
  ${CTEST_OUTPUT_DIR}/mmg3d_OptimAni_Sphere.o.mesh
  )

ADD_TEST(NAME mmg3d_OptimAni_Cube
  COMMAND ${EXECUT_MMG3D} -v 5 -optim -A -hgrad -1
  ${MMG3D_CI_TESTS}/Cube/cube-ani
  -out ${CTEST_OUTPUT_DIR}/mmg3d_OptimAni_cube.o.meshb)


##############################################################################
#####
#####         Check optimLES
#####
##############################################################################
#####
ADD_TEST(NAME mmg3d_OptimLES_sphere
  COMMAND ${EXECUT_MMG3D} -v 5 -optimLES
  ${MMG3D_CI_TESTS}/SphereIso_0.25h_met/SphereIso0.5
  ${CTEST_OUTPUT_DIR}/mmg3d_OptimLES_Sphere.o.mesh
  )

###############################################################################
#####
#####         Check Results
#####
###############################################################################
#####

# lssurf: discretization of boundaries only
ADD_TEST(NAME mmg3d_OptLsSurf_box
  COMMAND ${EXECUT_MMG3D} -v 5 -lssurf
  -sol ${MMG3D_CI_TESTS}/OptLsSurf_box/box.sol
  ${MMG3D_CI_TESTS}/OptLsSurf_box/box.mesh
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLsSurf_box.o.meshb
  )

# multi-mat
ADD_TEST(NAME mmg3d_LSMultiMat
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -nr
  ${MMG3D_CI_TESTS}/LSMultiMat/step.0.mesh
  -sol ${MMG3D_CI_TESTS}/LSMultiMat/step.0.phi.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_LSMultiMat.o.meshb)

#multi-mat + opnbdy + non-manifold check
ADD_TEST(NAME mmg3d_LSMultiMat_nm
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -0.1 -hausd 0.05 -hgrad 1.8 -nr -opnbdy
  ${MMG3D_CI_TESTS}/LSMultiMat/3d-opn
  ${CTEST_OUTPUT_DIR}/mmg3d_3d-opn.o.meshb)

ADD_TEST(NAME mmg3d_OptLs_plane_val
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -val
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  -met ${MMG3D_CI_TESTS}/OptLs_plane/met.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

#ADD_TEST(NAME mmg3d_OptLs_plane_default
#  COMMAND ${EXECUT_MMG3D} -v 5 -ls -default
#  ${MMG3D_CI_TESTS}/OptLs_plane/plane
#  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
#  -met ${MMG3D_CI_TESTS}/OptLs_plane/met.sol
#  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

SET_PROPERTY(TEST  mmg3d_OptLs_plane_val #mmg3d_OptLs_plane_default
  PROPERTY WILL_FAIL TRUE)

# ls oritentation
ADD_TEST(NAME mmg3d_OptLs_plane_p
  COMMAND ${EXECUT_MMG3D} -v 5 -ls
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/p.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-p.o.meshb)

ADD_TEST(NAME mmg3d_OptLs_plane_m
  COMMAND ${EXECUT_MMG3D} -v 5 -ls
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-m.o.meshb)

# ridge preservation
IF ( (NOT SCOTCH_FOUND) OR USE_SCOTCH MATCHES OFF )
  SET ( DISABLE_RENUM "" )
ELSE()
  SET ( DISABLE_RENUM -rn 0 )
ENDIF()

ADD_TEST(NAME mmg3d_OptLs_NM_ridge
  COMMAND ${EXECUT_MMG3D} -v 5 -ls 0.5 -noinsert -noswap -nomove -nr ${DISABLE_RENUM}
  ${MMG3D_CI_TESTS}/OptLs_NM_ridge/cube-it2.mesh
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_NM_cube-it2.o.mesh)

SET_TESTS_PROPERTIES ( mmg3d_OptLs_NM_ridge
  PROPERTIES FIXTURES_SETUP mmg3d_OptLs_NM_ridge )

# non-zero ls
ADD_TEST(NAME mmg3d_OptLs_plane_nonzero
  COMMAND ${EXECUT_MMG3D} -v 5 -ls 0.1
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

# ls discretization + optim
ADD_TEST(NAME mmg3d_OptLs_plane_optim
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -optim
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

# ls discretization + optim + aniso
ADD_TEST(NAME mmg3d_OptLs_plane_optimAni
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -optim -A
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

# ls discretization + hsiz
ADD_TEST(NAME mmg3d_OptLs_plane_hsiz
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -hsiz 0.2
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

# ls discretization + hsiz
ADD_TEST(NAME mmg3d_OptLs_plane_hsizAni
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -hsiz 0.2 -A
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

# ls discretization + metric
ADD_TEST(NAME mmg3d_OptLs_plane_withMetAndLs
  COMMAND ${EXECUT_MMG3D} -v 5 -ls
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/m.sol
  -met ${MMG3D_CI_TESTS}/OptLs_plane/met.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-nonzero.o.meshb)

# ls + rmc + LSBaseReference
ADD_TEST(NAME mmg3d_OptLs_LSBaseReferences-rmc
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -rmc -nr
  ${MMG3D_CI_TESTS}/LSBaseReferences/box
  -sol ${MMG3D_CI_TESTS}/LSBaseReferences/box.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_LSBaseReferences-rmc.o.meshb)

ADD_TEST(NAME mmg3d_OptLs_LSBaseReferences-normc
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -nr
  ${MMG3D_CI_TESTS}/LSBaseReferences/box
  -sol ${MMG3D_CI_TESTS}/LSBaseReferences/box.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_LSBaseReferences-normc.o.meshb)

# ls + rmc
ADD_TEST(NAME mmg3d_OptLs_plane_withbub
  COMMAND ${EXECUT_MMG3D} -v 5 -ls
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/bub.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-withbub.o.meshb)

# ls + rmc: max pile bug
ADD_TEST(NAME mmg3d_OptLs_plane_rmcmaxpile
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -rmc
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/whole.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-rmcmaxpile.o.meshb)

ADD_TEST(NAME mmg3d_OptLs_plane_rembub
  COMMAND ${EXECUT_MMG3D} -v 5 -ls
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/bub.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-rembub.o.meshb -rmc)

ADD_TEST(NAME mmg3d_OptLs_plane_rembub2
  COMMAND ${EXECUT_MMG3D} -v 5 -ls -rmc 0.1
  ${MMG3D_CI_TESTS}/OptLs_plane/plane
  -sol ${MMG3D_CI_TESTS}/OptLs_plane/bub.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_plane-rembub2.o.meshb)

# Preservation of orphan points
ADD_TEST(NAME mmg3d_OptLs_temp_orphan
  COMMAND ${EXECUT_MMG3D} -v 5 -ls
  ${MMG3D_CI_TESTS}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp
  -sol ${MMG3D_CI_TESTS}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp.sol
  -hausd 0.5 -nr -hgrad -1 -nsd 3
  ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_temp_orphan.o.meshb)

# OptLs and isoref option: compare the result of ls discretization with ref 10
# and results of the same case with ref 5
#include(FindUnixCommands)

add_test(
  NAME mmg3d_OptLs_isoref_defaut
  COMMAND ${EXECUT_MMG3D} -v 5 -ls ${MMG3D_CI_TESTS}/OptLs_isoref/3d-mesh.mesh
  -sol ${MMG3D_CI_TESTS}/OptLs_isoref/3d-mesh.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_isoref.o.mesh
  )
add_test(
  NAME mmg3d_OptLs_isoref_5
  COMMAND ${EXECUT_MMG3D} -v 5 -isoref 5 -ls
  ${MMG3D_CI_TESTS}/OptLs_isoref/3d-mesh-isoref5.mesh
  -sol ${MMG3D_CI_TESTS}/OptLs_isoref/3d-mesh.sol
  ${CTEST_OUTPUT_DIR}/mmg3d_isoref5.o.mesh
  )

if (BASH)
  add_test(
    NAME mmg3d_optLs_isoref
    COMMAND ${BASH} -c "diff <(wc -wl ${CTEST_OUTPUT_DIR}/mmg3d_isoref.o.mesh  | awk '{print $1 $2}') <(wc -wl ${CTEST_OUTPUT_DIR}/mmg3d_isoref5.o.mesh | awk '{print $1 $2}')"
    )
endif()

ADD_TEST(NAME test_para_tria
  COMMAND ${EXECUT_MMG3D}
  -ar 0.02 -nofem -nosizreq -hgradreq -1 -hgrad -1
  ${MMG3D_CI_TESTS}/test_para_tria/proc0.mesh
  -sol ${MMG3D_CI_TESTS}/test_para_tria/proc0.sol
  ${CTEST_OUTPUT_DIR}/proc0.o.mesh
  )

SET_TESTS_PROPERTIES ( test_para_tria
  PROPERTIES FIXTURES_SETUP test_para_tria )


IF ( LONG_TESTS )
  # Test the Ls option
  ADD_TEST(NAME mmg3d_OptLs_cube303d_hminMax_hgrad1.2_hausd0.005
    COMMAND ${EXECUT_MMG3D} -v 5 -ls
    ${MMG3D_CI_TESTS}/OptLs_cube303d_hminMax_hgrad1.2_hausd0.005/cube303d
    -sol ${MMG3D_CI_TESTS}/OptLs_cube303d_hminMax_hgrad1.2_hausd0.005/cube303d.sol
    -hausd 0.005 -nr -hgrad 1.2 -hmin 0.001 -hmax 0.1
    ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_cube303d_hminMax_hgrad1.2_hausd0.005-cube303d.o.meshb)
  ADD_TEST(NAME mmg3d_OptLs_temp_hminMax_hgrad1.2_hausd0.1
    COMMAND ${EXECUT_MMG3D} -v 5 -ls
    ${MMG3D_CI_TESTS}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp
    -sol ${MMG3D_CI_TESTS}/OptLs_temp_hminMax_hgrad1.2_hausd0.1/temp.sol
    -hausd 0.1 -nr -hgrad 1.2 -hmin 3 -hmax 4
    ${CTEST_OUTPUT_DIR}/mmg3d_OptLs_temp_hminMax_hgrad1.2_hausd0.1-temp.o.meshb)

  ###############################################################################
  #####
  #####         Check Lagrangian motion option
  #####
  ###############################################################################
  #####
  IF ( ELAS_FOUND AND NOT USE_ELAS MATCHES OFF )
    ADD_TEST(NAME mmg3d_LagMotion0_boxt
      COMMAND ${EXECUT_MMG3D} -v 5  -lag 0
      -in ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt
      -sol ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt.sol
      -out ${CTEST_OUTPUT_DIR}/mmg3d_LagMotion0_boxt-boxt.o.meshb
      )
    ADD_TEST(NAME mmg3d_LagMotion1_boxt
      COMMAND ${EXECUT_MMG3D} -v 5  -lag 1
      -in ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt
      -sol ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt.sol
      -out ${CTEST_OUTPUT_DIR}/mmg3d_LagMotion1_boxt-boxt.o.meshb
      )
    ADD_TEST(NAME mmg3d_LagMotion2_boxt
      COMMAND ${EXECUT_MMG3D} -v 5  -lag 2
      -in ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt
      -sol ${MMG3D_CI_TESTS}/LagMotion1_boxt/boxt.sol
      -out ${CTEST_OUTPUT_DIR}/mmg3d_LagMotion2_boxt-boxt.o.meshb
      )
  ENDIF()

ENDIF()

###############################################################################
#####
#####         Bug Fix
#####
###############################################################################
#####
#ADD_TEST(NAME mmg3d_BUG_OptLsSingularities
# COMMAND ${EXECUT_MMG3D} -v 5  -ls
# ${MMG3D_CI_TESTS}/BUG_OptLsSingularities/test4
# ${CTEST_OUTPUT_DIR}/mmg3d_BUG_OptLsSingularities-test4.o.meshb)
#
#ADD_TEST(NAME mmg3d_TestDoSol_1
# COMMAND ${EXECUT_MMG3D} -v 5  -hgrad -1 -hausd 1 -m 100
# ${MMG3D_CI_TESTS}/TestDoSol_1/66_shaver3.mesh
# ${CTEST_OUTPUT_DIR}/mmg3d_TestDoSol_1-66_shaver3.o.meshb)
