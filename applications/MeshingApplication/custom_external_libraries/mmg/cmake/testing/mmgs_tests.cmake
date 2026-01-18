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
SET ( test_names mmgs_SimpleTeapot )
SET ( input_files ${MMGS_CI_TESTS}/Teapot/teapot )
SET ( args  "-v 5" )
SET ( common_args "" )

ADD_RUN_AGAIN_TESTS ( ${EXECUT_MMGS} "${test_names}" "${args}" "${input_files}" )

ADD_TEST(NAME mmgs_CubeAni
  COMMAND ${EXECUT_MMGS}
  ${MMGS_CI_TESTS}/CubeAni/cube
  -out ${CTEST_OUTPUT_DIR}/mmgs_CubeAni-cube.d.meshb)

ADD_TEST(NAME mmgs_SphereAni
  COMMAND ${EXECUT_MMGS}
  ${MMGS_CI_TESTS}/SphereAni/sphere
  -out ${CTEST_OUTPUT_DIR}/mmgs_SphereAni-sphere.d.meshb)

###############################################################################
#####
#####         Options
#####
###############################################################################

ADD_TEST(NAME mmgs_memOption
  COMMAND ${EXECUT_MMGS} -v 5 -m 100 ${common_args}
  ${MMGS_CI_TESTS}/Teapot/teapot
  -out ${CTEST_OUTPUT_DIR}/mmgs_memOption.o.meshb)

ADD_TEST(NAME mmgs_val
  COMMAND ${EXECUT_MMGS} -val
  ${MMGS_CI_TESTS}/Teapot/teapot
  )
SET_PROPERTY(TEST mmgs_val
  PROPERTY WILL_FAIL TRUE)

# nsd
ADD_TEST(NAME mmgs_nsd24
  COMMAND ${EXECUT_MMGS} -v 5 -nsd 24 ${common_args}
  ${MMGS_CI_TESTS}/Teapot/teapot
  -out ${CTEST_OUTPUT_DIR}/mmgs_nsd24.o.meshb)

ADD_TEST(NAME mmgs_locParamCrea
  COMMAND ${EXECUT_MMGS} -v 5 -default
  ${MMGS_CI_TESTS}/LocParamsCrea/circle2refs.mesh)
SET_TESTS_PROPERTIES ( mmgs_locParamCrea
  PROPERTIES FIXTURES_SETUP mmgs_locParamCrea )
ADD_TEST(NAME mmgs_locParamClean
  COMMAND ${CMAKE_COMMAND} -E remove -f
  ${MMGS_CI_TESTS}/LocParamsCrea/circle2refs.mmgs)
SET_TESTS_PROPERTIES ( mmgs_locParamClean
  PROPERTIES FIXTURES_REQUIRED mmgs_locParamCrea )

###############################################################################
#####
#####         Input/Output
#####
###############################################################################

# Binary gmsh
ADD_TEST(NAME mmgs_binary_gmsh_s
  COMMAND ${EXECUT_MMGS} -v 5 ${common_args}
  ${MMGS_CI_TESTS}/GmshInout/cube.mshb
  ${CTEST_OUTPUT_DIR}/)

# Ascii gmsh
ADD_TEST(NAME mmgs_ascii_gmsh_s
  COMMAND ${EXECUT_MMGS} -v 5 ${common_args}
  ${MMGS_CI_TESTS}/GmshInout/cube.msh
  ${CTEST_OUTPUT_DIR}/mmgs-cube-gmsh.o.msh)

# VTK .vtk no metric
ADD_TEST(NAME mmgs_vtkvtk
COMMAND ${EXECUT_MMGS} -v 5
${MMGS_CI_TESTS}/VtkInout/c1.vtk
${CTEST_OUTPUT_DIR}/mmgs_vtkvtk)

# VTK .vtu no metric
ADD_TEST(NAME mmgs_vtkvtu
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/VtkInout/c1.vtu
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtu)

# VTK .vtp no metric
ADD_TEST(NAME mmgs_vtkvtp
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/VtkInout/c1.vtp
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtp)

# VTK .vtk with iso metric
ADD_TEST(NAME mmgs_vtkvtk_iso
COMMAND ${EXECUT_MMGS} -v 5
${MMGS_CI_TESTS}/VtkInout/iso.vtk
${CTEST_OUTPUT_DIR}/mmgs_vtkvtk_iso)

# VTK .vtu with iso metric
ADD_TEST(NAME mmgs_vtkvtu_iso
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/VtkInout/iso.vtu
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtu_iso)

  # VTK .vtp with iso metric
ADD_TEST(NAME mmgs_vtkvtp_iso
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/VtkInout/iso.vtp
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtp_iso)

# VTK .vtk with aniso metric
ADD_TEST(NAME mmgs_vtkvtk_ani
COMMAND ${EXECUT_MMGS} -v 5
${MMGS_CI_TESTS}/VtkInout/ani.vtk
${CTEST_OUTPUT_DIR}/mmgs_vtkvtk_ani)

# VTK .vtu with aniso metric
ADD_TEST(NAME mmgs_vtkvtu_ani
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/VtkInout/ani.vtu
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtu_ani)

# VTK .vtp with aniso metric
ADD_TEST(NAME mmgs_vtkvtp_ani
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMGS_CI_TESTS}/VtkInout/ani.vtp
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtp_ani)

# VTK .vtk with ls
ADD_TEST(NAME mmgs_vtkvtk_ls
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.1
  ${MMGS_CI_TESTS}/VtkInout/surf_ls.vtk
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtk_ls)

# VTK .vtu with ls
ADD_TEST(NAME mmgs_vtkvtu_ls
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.1
  ${MMGS_CI_TESTS}/VtkInout/surf_ls.vtu
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtu_ls)

# VTK .vtp with ls
ADD_TEST(NAME mmgs_vtkvtp_ls
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.1
  ${MMGS_CI_TESTS}/VtkInout/surf_ls.vtp
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtp_ls)

# VTK .vtk with ls and metric
ADD_TEST(NAME mmgs_vtkvtk_ls_metric
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.1
  ${MMGS_CI_TESTS}/VtkInout/surf_ls_metric.vtk
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtk_ls_metric)

# VTK .vtu with ls and metric
ADD_TEST(NAME mmgs_vtkvtu_ls_metric
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.1
  ${MMGS_CI_TESTS}/VtkInout/surf_ls_metric.vtu
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtu_ls_metric)

# VTK .vtp with ls and metric
ADD_TEST(NAME mmgs_vtkvtp_ls_metric
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.1
  ${MMGS_CI_TESTS}/VtkInout/surf_ls_metric.vtp
  ${CTEST_OUTPUT_DIR}/mmgs_vtkvtp_ls_metric)

# VTK .vtk with metric and ls
ADD_TEST(NAME mmgs_vtkvtk_metric_ls
COMMAND ${EXECUT_MMGS} -v 5 -ls 0.1
${MMGS_CI_TESTS}/VtkInout/surf_metric_ls.vtk
${CTEST_OUTPUT_DIR}/mmgs_vtkvtk_metric_ls)

IF ( (NOT VTK_FOUND) OR USE_VTK MATCHES OFF )
  SET(expr "VTK library not founded")
  SET_PROPERTY(TEST mmgs_vtkvtk
    mmgs_vtkvtu
    mmgs_vtkvtp
    mmgs_vtkvtk_iso
    mmgs_vtkvtu_iso
    mmgs_vtkvtp_iso
    mmgs_vtkvtk_ani
    mmgs_vtkvtu_ani
    mmgs_vtkvtp_ani
    mmgs_vtkvtk_ls
    mmgs_vtkvtu_ls
    mmgs_vtkvtp_ls
    mmgs_vtkvtk_ls_metric
    mmgs_vtkvtu_ls_metric
    mmgs_vtkvtp_ls_metric
    mmgs_vtkvtk_metric_ls
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
  ENDIF ( )

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
ADD_TEST(NAME mmgs_Rhino_M
  COMMAND ${EXECUT_MMGS} -v 5 ${common_args}
  ${MMGS_CI_TESTS}/Rhino_M/rhino -hausd 1
  -out ${CTEST_OUTPUT_DIR}/mmgs_Rhino_M-rhino.d.meshb)

ADD_TEST(NAME mmgs_moebius
  COMMAND ${EXECUT_MMGS} -v 5 ${common_args} -d
  ${MMGS_CI_TESTS}/Moebius-strip/moebius-strip.mesh -nr
  -out ${CTEST_OUTPUT_DIR}/mmgs_moebius-strip.d.mesh)

###############################################################################
#####
#####         Non manifold cases
#####
###############################################################################
ADD_TEST(NAME mmgs_Cow_NM_hausd10
  COMMAND ${EXECUT_MMGS} -v 5 ${common_args}
  ${MMGS_CI_TESTS}/Cow_NM/cow -hausd 10
  -out ${CTEST_OUTPUT_DIR}/mmgs_Cow_NM_hausd10-cow.d.meshb)

###############################################################################
#####
#####         Test results
#####
###############################################################################
# Test the Ls option
ADD_TEST(NAME mmgs_OptLs_val
  COMMAND ${EXECUT_MMGS} -v 5 -ls -val
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${CTEST_OUTPUT_DIR}/mmgs_teapot-val.o.meshb)
#ADD_TEST(NAME mmgs_OptLs_default
#  COMMAND ${EXECUT_MMGS} -v 5 -ls -default
#  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
#  ${CTEST_OUTPUT_DIR}/mmgs_teapot-val.o.meshb)

SET_PROPERTY(TEST mmgs_OptLs_val #mmgs_OptLs_default
  PROPERTY WILL_FAIL TRUE)


ADD_TEST(NAME mmgs_OptLs_teapot
  COMMAND ${EXECUT_MMGS} -v 5 -ls ${common_args}
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${CTEST_OUTPUT_DIR}/mmgs_OptLs_teapot-teapot.simple.o.meshb)

ADD_TEST(NAME mmgs_OptLs_teapot_keepRef
  COMMAND ${EXECUT_MMGS} -v 5 -ls -keep-ref ${common_args}
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${CTEST_OUTPUT_DIR}/mmgs_OptLs_teapot_keepRef-teapot.keep-ref.o.meshb)

ADD_TEST(NAME mmgs_OptLs_teapot_0.5_keepRef
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.5 -keep-ref ${common_args}
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${CTEST_OUTPUT_DIR}/mmgs_OptLs_teapot_0.5_keepRef-teapot.0.5.keep-ref.o.meshb)

ADD_TEST(NAME mmgs_OptLs_teapot2
  COMMAND ${EXECUT_MMGS} -v 5 -ls -nr ${common_args}
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${CTEST_OUTPUT_DIR}/mmgs_OptLs_teapot2-teapot.o.meshb)

add_test(
  NAME mmgs_OptLs_isoref_defaut
  COMMAND ${EXECUT_MMGS} -v 5 -ls ${MMGS_CI_TESTS}/OptLs_isoref/surf-mesh.mesh
  -sol ${MMGS_CI_TESTS}/OptLs_isoref/surf-mesh.sol
  ${CTEST_OUTPUT_DIR}/mmgs_isoref.o.mesh
  )
add_test(
  NAME mmgs_OptLs_isoref_5
  COMMAND ${EXECUT_MMGS} -v 5 -isoref 5 -ls
  ${MMGS_CI_TESTS}/OptLs_isoref/surf-mesh-isoref5.mesh
  -sol ${MMGS_CI_TESTS}/OptLs_isoref/surf-mesh.sol
  ${CTEST_OUTPUT_DIR}/mmgs_isoref5.o.mesh
  )

if (BASH)
  add_test(
    NAME mmgs_optLs_isoref
    COMMAND ${BASH} -c "diff <(wc -wl ${CTEST_OUTPUT_DIR}/mmgs_isoref.o.mesh  | awk '{print $1 $2}') <(wc -wl ${CTEST_OUTPUT_DIR}/mmgs_isoref5.o.mesh | awk '{print $1 $2}')"
    )
endif()


####### -met option
ADD_TEST(NAME mmgs_2squares-withMet
  COMMAND ${EXECUT_MMGS} -v 5 -d
  ${MMG2D_CI_TESTS}/2squares/2squares -met ${MMG2D_CI_TESTS}/2squares/2s.sol
  -out ${CTEST_OUTPUT_DIR}/mmgs_2squares-met.o.meshb)

####### -sol option
ADD_TEST(NAME mmgs_2squares-withSol
  COMMAND ${EXECUT_MMGS} -v 5
  ${MMG2D_CI_TESTS}/2squares/2squares -sol ${MMG2D_CI_TESTS}/2squares/2s.sol
  -out ${CTEST_OUTPUT_DIR}/mmgs_2squares-sol.o.meshb)

####### orphan points
ADD_TEST(NAME mmgs_2squares-orphan
  COMMAND ${EXECUT_MMGS} -v 5 -nsd 10
  ${MMG2D_CI_TESTS}/2squares/2squares
  -out ${CTEST_OUTPUT_DIR}/mmgs_2squares-orphan.o.meshb)


# nsd + ls
ADD_TEST(NAME mmgs_OptLs_teapot-nsd3
  COMMAND ${EXECUT_MMGS} -v 5 -ls -nsd 3 ${common_args}
  ${MMGS_CI_TESTS}/OptLs_teapot/teapot
  ${CTEST_OUTPUT_DIR}/mmgs_OptLs_teapot-ls-nsd3.o.meshb)

###############################################################################
#####
#####         Detected Bugs
#####
###############################################################################
ADD_TEST(NAME mmgs_Car_NM
  COMMAND ${EXECUT_MMGS} -v 5 ${common_args}
  ${MMGS_CI_TESTS}/Car_NM/car
  -out ${CTEST_OUTPUT_DIR}/mmgs_Car_NM-car.d.meshb)

ADD_TEST(NAME mmgs_Cow_NM_hausd20
  COMMAND ${EXECUT_MMGS} -v 5 ${common_args}
  ${MMGS_CI_TESTS}/Cow_NM/cow -hausd 20
  -out ${CTEST_OUTPUT_DIR}/mmgs_Cow_NM_hausd20-cow.d.meshb)

###############################################################################
#####
#####         Implicit domain discretization
#####
###############################################################################
ADD_TEST(NAME mmgs_LSMultiMat_val
  COMMAND ${EXECUT_MMGS} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3 -val
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  -met ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-val.o.meshb
  )
#ADD_TEST(NAME mmgs_LSMultiMat_default
#  COMMAND ${EXECUT_MMGS} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3 -default
#  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
#  -met ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-met.sol
#  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
#  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-default.o.meshb)
SET_PROPERTY(TEST mmgs_LSMultiMat_val #mmgs_LSMultiMat_default
  PROPERTY WILL_FAIL TRUE)


ADD_TEST(NAME mmgs_LSMultiMat
  COMMAND ${EXECUT_MMGS} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat.o.meshb)

ADD_TEST(NAME mmgs_LSMultiMat-rmc
  COMMAND ${EXECUT_MMGS} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3 -rmc
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat.o.meshb)

# ls + rmc + LSBaseReference
ADD_TEST(NAME mmgs_OptLs_LSBaseReferences-rmc
  COMMAND ${EXECUT_MMGS} -v 5 -ls -rmc
  ${MMGS_CI_TESTS}/LSBaseReferences/box
  -sol ${MMGS_CI_TESTS}/LSBaseReferences/box.sol
  ${CTEST_OUTPUT_DIR}/mmgs_OptLs_LSBaseReferences-rmc.o.meshb)

ADD_TEST(NAME mmgs_OptLs_LSBaseReferences-normc
  COMMAND ${EXECUT_MMGS} -v 5 -ls
  ${MMGS_CI_TESTS}/LSBaseReferences/box
  -sol ${MMGS_CI_TESTS}/LSBaseReferences/box.sol
  ${CTEST_OUTPUT_DIR}/mmgs_OptLs_LSBaseReferences-normc.o.meshb)

# non 0 ls
ADD_TEST(NAME mmgs_LSMultiMat_nonzero
  COMMAND ${EXECUT_MMGS} -v 5 -ls 0.01 -hausd 0.001
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-nonzero.o.meshb)

# optim + aniso option with corners
ADD_TEST(NAME mmgs_CubeOptimAni
  COMMAND ${EXECUT_MMGS} -v 5 -optim -A -hgrad -1
  ${MMGS_CI_TESTS}/CubeOptimAni/cube-ani
  -out ${CTEST_OUTPUT_DIR}/mmgs_CubeOptimAni-cube.o.meshb)

# optim + aniso option for open surface
ADD_TEST(NAME mmgs_OpnbdyOptimAni-circle
  COMMAND ${EXECUT_MMGS} -v 5 -optim -A -hgrad -1
  ${MMGS_CI_TESTS}/OpnbdyOptimAni/cercle-3D.mesh
  -out ${CTEST_OUTPUT_DIR}/mmgs_OpnbdyOptimAni-circle.o.meshb)

ADD_TEST(NAME mmgs_OpnbdyOptimAni-adap1
  COMMAND ${EXECUT_MMGS} -v 5 -optim -A -hgrad -1
  ${MMGS_CI_TESTS}/OpnbdyOptimAni/adap1-3D.mesh
  -out ${CTEST_OUTPUT_DIR}/mmgs_OpnbdyOptimAni-adap1.o.meshb)

# ls discretisation + parameter file
ADD_TEST(NAME mmgs_ParsOpName
COMMAND ${EXECUT_MMGS} -v 5 -ls
-sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
-f ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-refs.mmgs
${MMGS_CI_TESTS}/LSMultiMat/multi-mat.mesh
${CTEST_OUTPUT_DIR}/mmgs_ParsOpName.o.meshb)

SET(parsopName "multi-mat-refs.mmgs OPENED")
SET_PROPERTY(TEST mmgs_ParsOpName
PROPERTY PASS_REGULAR_EXPRESSION "${parsopName}")

# ls discretisation + wrong name of parameter file
ADD_TEST(NAME mmgs_ParsOpName_wrongFile
COMMAND ${EXECUT_MMGS} -v 5 -ls
-sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
-f ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-false.mmg
${MMGS_CI_TESTS}/LSMultiMat/multi-mat.mesh
${CTEST_OUTPUT_DIR}/mmgs_ParsOpName_wrongFile.o.meshb)

SET(parsopNameWrong "multi-mat-false.mmgs file NOT FOUND.")
SET_PROPERTY(TEST mmgs_ParsOpName_wrongFile
PROPERTY PASS_REGULAR_EXPRESSION "${parsopNameWrong}")

# ls discretisation + no name of parameter file
ADD_TEST(NAME mmgs_ParsOpName_NoFileName
COMMAND ${EXECUT_MMGS} -v 5 -f -ls
-sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
${MMGS_CI_TESTS}/LSMultiMat/multi-mat.mesh
${CTEST_OUTPUT_DIR}/mmgs_ParsOpName_NoFileName.o.meshb)

SET(parsopNameNo "Missing filename for -f")
SET_PROPERTY(TEST mmgs_ParsOpName_NoFileName
PROPERTY PASS_REGULAR_EXPRESSION "${parsopNameNo}")

# ls discretisation + optim option
ADD_TEST(NAME mmgs_LSMultiMat_optim
  COMMAND ${EXECUT_MMGS} -v 5 -ls -optim -hausd 0.001
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-optim.o.meshb)

# ls discretisation + optim + aniso option
ADD_TEST(NAME mmgs_LSMultiMat_optimAni
  COMMAND ${EXECUT_MMGS} -v 5 -ls -optim -A -hausd 0.001
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-optimAni.o.meshb)

SET(passRegex "## ERROR: MISMATCH OPTIONS: OPTIM OPTION IS NOT AVAILABLE")

SET_PROPERTY(TEST mmgs_LSMultiMat_optim
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
SET_PROPERTY(TEST mmgs_LSMultiMat_optimAni
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

# ls discretisation + hsiz option
ADD_TEST(NAME mmgs_LSMultiMat_hsiz
  COMMAND ${EXECUT_MMGS} -v 5 -ls -hsiz 0.05 -hausd 0.001
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-hsiz.o.meshb)

# ls discretisation + hsiz Ani option
ADD_TEST(NAME mmgs_LSMultiMat_hsizAni
  COMMAND ${EXECUT_MMGS} -v 5 -ls -hsiz 0.05 -A -hausd 0.001
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-hsizAni.o.meshb)

# ls discretisation + metric
ADD_TEST(NAME mmgs_LSMultiMat_withMet
  COMMAND ${EXECUT_MMGS} -v 5 -ls -hausd 0.001
  -met ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-withMet.o.meshb)

# ls discretisation + metric + ls
ADD_TEST(NAME mmgs_LSMultiMat_withMetAndLs
  COMMAND ${EXECUT_MMGS} -v 5 -ls -hausd 0.001
  -met ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMGS_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${MMGS_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmgs_LSMultiMat-withMetAndLs.o.meshb)

# ls discretization with wrong orientation of input triangles
ADD_TEST(NAME mmgs_LSTriaOri
  COMMAND ${EXECUT_MMGS} -v 5 -ls -hausd 0.001
  ${MMGS_CI_TESTS}/LSTriaOri/fault.mesh
  ${CTEST_OUTPUT_DIR}/mmgs_LSTriaOri.o.meshb)

# lssurf: discretization of boundaries only
ADD_TEST(NAME mmgs_OptLsSurf_box
  COMMAND ${EXECUT_MMGS} -v 5 -lssurf
  -sol ${MMGS_CI_TESTS}/OptLsSurf_box/box.sol
  ${MMGS_CI_TESTS}/OptLsSurf_box/box-3D.mesh
  ${CTEST_OUTPUT_DIR}/mmgs_OptLsSurf_box.o.meshb
  )

# lssurf + multimat: discretization of boundaries only
ADD_TEST(NAME mmgs_OptLsSurf_multiMat_box
  COMMAND ${EXECUT_MMGS} -v 5 -lssurf
  -sol ${MMGS_CI_TESTS}/OptLsSurf_box/box.sol
  ${MMGS_CI_TESTS}/OptLsSurf_box/box_multiMat-3D.mesh
  ${CTEST_OUTPUT_DIR}/mmgs_OptLsSurf_multiMat_box.o.meshb
  )


###############################################################################
#####
#####         Check snapping (prevision of non-manifold situations)
#####
###############################################################################
#####
SET(nmRegex "unsnap at least 1 point")

ADD_TEST(NAME mmgs_LSSnapval_manifold1
  COMMAND ${EXECUT_MMGS} -v 5  -ls
  -in ${MMGS_CI_TESTS}/LSSnapval/8elts1.mesh
  -sol ${MMGS_CI_TESTS}/LSSnapval/manifold.sol
  -out ${CTEST_OUTPUT_DIR}/mmgs_LSSnapval_manifold1.o.mesh
  )

ADD_TEST(NAME mmgs_LSSnapval_manifold2
  COMMAND ${EXECUT_MMGS} -v 5  -ls
  -in ${MMGS_CI_TESTS}/LSSnapval/8elts2.mesh
  -sol ${MMGS_CI_TESTS}/LSSnapval/manifold.sol
  -out ${CTEST_OUTPUT_DIR}/mmgs_LSSnapval_manifold2.o.mesh
  )

SET_PROPERTY(TEST mmgs_LSSnapval_manifold1 mmgs_LSSnapval_manifold2
  PROPERTY FAIL_REGULAR_EXPRESSION "${nmRegex}")

ADD_TEST(NAME mmgs_LSSnapval_non-manifold1
  COMMAND ${EXECUT_MMGS} -v 5  -ls
  -in ${MMGS_CI_TESTS}/LSSnapval/8elts1.mesh
  -sol ${MMGS_CI_TESTS}/LSSnapval/8elts1-nm.sol
  -out ${CTEST_OUTPUT_DIR}/mmgs_LSSnapval_non-manifold1.o.mesh
  )

ADD_TEST(NAME mmgs_LSSnapval_non-manifold2
  COMMAND ${EXECUT_MMGS} -v 5  -ls
  -in ${MMGS_CI_TESTS}/LSSnapval/8elts2.mesh
  -sol ${MMGS_CI_TESTS}/LSSnapval/8elts2-nm.sol
  -out ${CTEST_OUTPUT_DIR}/mmgs_LSSnapval_non-manifold2.o.mesh
  )
SET_PROPERTY(TEST mmgs_LSSnapval_non-manifold1 mmgs_LSSnapval_non-manifold2
  PROPERTY PASS_REGULAR_EXPRESSION "${nmRegex}")
