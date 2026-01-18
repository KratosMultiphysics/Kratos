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

# Simple test: must already pass (-d option allows to cover chkmsh function)
ADD_TEST(NAME mmg2d_SimpleCircle
  COMMAND ${EXECUT_MMG2D} -v 5 -d
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_SimpleCircle-cercle.o.meshb)

###############################################################################
#####
#####         Options
#####
###############################################################################
ADD_TEST(NAME mmg2d_help
  COMMAND ${EXECUT_MMG2D} -h
  )
SET_PROPERTY(TEST mmg2d_help
  PROPERTY PASS_REGULAR_EXPRESSION "File specifications")

ADD_TEST(NAME mmg2d_memOption
  COMMAND ${EXECUT_MMG2D} -v 5 -m 100
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_memOption.o.meshb)

ADD_TEST(NAME mmg2d_val
  COMMAND ${EXECUT_MMG2D} -v val
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_val.o.meshb)
SET_PROPERTY(TEST mmg2d_val
  PROPERTY WILL_FAIL TRUE)

ADD_TEST(NAME mmg2d_locParamCrea
  COMMAND ${EXECUT_MMG2D} -v 5 -default
  ${MMG2D_CI_TESTS}/LocParamsCrea/circle2refs.mesh)

SET_TESTS_PROPERTIES ( mmg2d_locParamCrea
  PROPERTIES FIXTURES_SETUP mmg2d_locParamCrea )
ADD_TEST(NAME mmg2d_locParamClean
  COMMAND ${CMAKE_COMMAND} -E remove -f
  ${MMG2D_CI_TESTS}/LocParamsCrea/circle2refs.mmg2d)
SET_TESTS_PROPERTIES ( mmg2d_locParamClean
  PROPERTIES FIXTURES_REQUIRED mmg2d_locParamCrea )

ADD_TEST(NAME mmg2d_hsizOption
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsiz-circle.o.meshb)

ADD_TEST(NAME mmg2d_hsizAni
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -A
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsizAni-circle.o.meshb)

ADD_TEST(NAME mmg2d_hsizAndNosurfOption
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -nosurf
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsizNosurf-circle.o.meshb)

ADD_TEST(NAME mmg2d_hsizAndNosurfAni
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -nosurf -A
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsizNosurfAni-circle.o.meshb)

ADD_TEST(NAME mmg2d_hsizAndNosurfOption2
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -nosurf -3dMedit 2
  ${MMG2D_CI_TESTS}/2squares/2squares
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsizNosurf-2squares.o.meshb)

ADD_TEST(NAME mmg2d_hsizHmax
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -hmax 0.05
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsizHmax-circle.o.meshb)
SET(passRegex "Mismatched options")
SET_PROPERTY(TEST mmg2d_hsizHmax
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

ADD_TEST(NAME mmg2d_hsizHmin
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -sol 2 -hmin 0.2
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_hsizHmin-circle.o.meshb)
SET(passRegex "Mismatched options")
SET_PROPERTY(TEST mmg2d_hsizHmin
  PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

ADD_TEST(NAME mmg2d_reqEntities-ref
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.02
  ${MMG2D_CI_TESTS}/Disk_ReqEntities/disk.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_reqEntities-ref.o.meshb)

ADD_TEST(NAME mmg2d_orphanPoint
  COMMAND ${EXECUT_MMG2D} -v 5 -hausd 10 -hgradreq -1 -nosizreq
  ${MMG2D_CI_TESTS}/Disk_ReqEntities/disk.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_orphan.o.meshb)

ADD_TEST(NAME mmg2d_reqEntitiesAni-ref
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.02 -A
  ${MMG2D_CI_TESTS}/Disk_ReqEntities/disk.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_reqEntitiesAni-ref.o.meshb)

ADD_TEST(NAME mmg2d_reqEntities-unref
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1
  ${MMG2D_CI_TESTS}/Disk_ReqEntities/disk-tiny.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_reqEntities-unref.o.meshb)

ADD_TEST(NAME mmg2d_reqEntitiesAni-unref
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.1 -A
  ${MMG2D_CI_TESTS}/Disk_ReqEntities/disk-tiny.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_reqEntitiesAni-unref.o.meshb)

ADD_TEST(NAME mmg2d_locParam
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/LocParams/circle2refs.mesh
  -out ${CTEST_OUTPUT_DIR}/locParams.o.meshb)

ADD_TEST(NAME mmg2d_locParam_ani
  COMMAND ${EXECUT_MMG2D} -v 5 -A
  ${MMG2D_CI_TESTS}/LocParams/circle2refs.mesh
  -out ${CTEST_OUTPUT_DIR}/locParams-ani.o.meshb)

ADD_TEST(NAME mmg2d_opnbdy_yes
  COMMAND ${EXECUT_MMG2D} -v 5 -opnbdy -hausd 0.001 -d
  ${MMG2D_CI_TESTS}/Opnbdy/opnbdy-mesh.msh
  -out ${CTEST_OUTPUT_DIR}/mmg2d-opnbdy-mesh-yes.o.meshb)

ADD_TEST(NAME mmg2d_opnbdy_no
  COMMAND ${EXECUT_MMG2D} -v 5 -hausd 0.001
  ${MMG2D_CI_TESTS}/Opnbdy/opnbdy-mesh.msh
  -out ${CTEST_OUTPUT_DIR}/mmg2d-opnbdy-mesh-no.o.meshb)

ADD_TEST(NAME mmg2d_opnbdy_ls
  COMMAND ${EXECUT_MMG2D} -v 5 -opnbdy -ls 3.4 -hausd 0.001 -d
  ${MMG2D_CI_TESTS}/Opnbdy/opnbdy.mesh
  -sol  ${MMG2D_CI_TESTS}/Opnbdy/ls.sol
  -out ${CTEST_OUTPUT_DIR}/mmg2d-opnbdy-ls.o.meshb)

ADD_TEST(NAME mmg2d_opnbdy_lssurf
  COMMAND ${EXECUT_MMG2D} -v 5 -opnbdy -lssurf 0.6
  ${MMG2D_CI_TESTS}/Opnbdy/opnbdy.mesh
  -sol  ${MMG2D_CI_TESTS}/Opnbdy/ls.sol
  -out ${CTEST_OUTPUT_DIR}/mmg2d-opnbdy-lssurf.o.meshb)

ADD_TEST(NAME mmg2d_opnbdy_yes_ani
  COMMAND ${EXECUT_MMG2D} -v 5 -hausd 0.001 -A -opnbdy
  ${MMG2D_CI_TESTS}/Opnbdy/opnbdy-mesh.msh
  -out ${CTEST_OUTPUT_DIR}/mmg2d-opnbdy-mesh-yes-ani.o.meshb)

# default hybrid
ADD_TEST(NAME mmg2d_hybrid_2d
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_hybrid_2d-default)

# hybrid opnbdy
ADD_TEST(NAME mmg2d_hybrid_opnbdy_2d
  COMMAND ${EXECUT_MMG2D} -v 5 -opnbdy
  ${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_hybrid_2d-opnbdy)

# hybrid hsiz
ADD_TEST(NAME mmg2d_hybrid_hsiz_2d
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.05 -hgradreq -1
  ${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_hybrid_2d-opnbdy)

ADD_TEST(NAME mmg2d_hybrid_nosizreq_nohgradreq_2d
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh -nosizreq -hgradreq -1
  ${CTEST_OUTPUT_DIR}/mmg2d_hybrid_2d-nosizreq)

# hybrid nsd: remove the triangular domain as it is of ref 1001
ADD_TEST(NAME mmg2d_hybrid-nsd1
  COMMAND ${EXECUT_MMG2D} -v 5 -nsd 1
  ${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_hybrid_2d-nsd1)

###############################################################################
#####
#####         Input/Output
#####
###############################################################################

# Binary gmsh no metric
ADD_TEST(NAME mmg2d_binary_gmsh_2d
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/GmshInout/cercle1.mshb
  ${CTEST_OUTPUT_DIR}/mmg2d_binary_gmsh_2d-cercle.mshb)

# Ascii gmsh no metric
ADD_TEST(NAME mmg2d_ascii_gmsh_2d
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/GmshInout/cercle1.msh
  ${CTEST_OUTPUT_DIR}/mmg2d_ascii_gmsh_2d-cercle)

# Ascii gmsh no metric hybrid
ADD_TEST(NAME mmg2d_gmsh_hybrid_2d
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/Hybrid/hybrid.msh
  ${CTEST_OUTPUT_DIR}/mmg2d_hybrid_gmsh_2d-hybrid)

# Binary gmsh iso metric
ADD_TEST(NAME mmg2d_binary_gmsh_iso
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/GmshInout/iso.mshb
  ${CTEST_OUTPUT_DIR}/mmg2d_binary_gmsh_iso.mshb)

# Ascii gmsh iso metric
ADD_TEST(NAME mmg2d_ascii_gmsh_iso
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/GmshInout/iso.msh
  ${CTEST_OUTPUT_DIR}/mmg2d_ascii_gmsh_iso)

# Binary gmsh iso metric
ADD_TEST(NAME mmg2d_binary_gmsh_ani
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/GmshInout/ani.mshb
  ${CTEST_OUTPUT_DIR}/mmg2d_binary_gmsh_ani.mshb)

# Ascii gmsh iso metric
ADD_TEST(NAME mmg2d_ascii_gmsh_ani
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/GmshInout/ani.msh
  ${CTEST_OUTPUT_DIR}/mmg2d_ascii_gmsh_ani)

# VTK .vtk no metric
ADD_TEST(NAME mmg2d_vtkvtk
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/cercle.vtk
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtk)

# VTK .vtp no metric
ADD_TEST(NAME mmg2d_vtkvtp
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/cercle.vtp
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtp)

# VTK .vtu no metric
ADD_TEST(NAME mmg2d_vtkvtu
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/cercle.vtu
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtu)

# VTK .vtk with iso metric
ADD_TEST(NAME mmg2d_vtkvtk_iso
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/iso.vtk
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtk_iso)

# VTK .vtp with iso metric
ADD_TEST(NAME mmg2d_vtkvtp_iso
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/iso.vtp
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtp_iso)

# VTK .vtu with iso metric
ADD_TEST(NAME mmg2d_vtkvtu_iso
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/iso.vtu
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtu_iso)

# VTK .vtk with aniso metric
ADD_TEST(NAME mmg2d_vtkvtk_ani
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/ani.vtk
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtk_ani)

# VTK .vtp with aniso metric
ADD_TEST(NAME mmg2d_vtkvtp_ani
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/ani.vtp
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtp_ani)

# VTK .vtu with aniso metric
ADD_TEST(NAME mmg2d_vtkvtu_ani
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/VtkInout/ani.vtu
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtu_ani)

# VTK .vtk with ls
ADD_TEST(NAME mmg2d_vtkvtk_ls
COMMAND ${EXECUT_MMG2D} -v 5 -ls 0.8
${MMG2D_CI_TESTS}/VtkInout/cercle_ls.vtk
${CTEST_OUTPUT_DIR}/mmg2d_vtkvtk_ls)

# VTK .vtu with ls
ADD_TEST(NAME mmg2d_vtkvtu_ls
COMMAND ${EXECUT_MMG2D} -v 5 -ls 0.8
${MMG2D_CI_TESTS}/VtkInout/cercle_ls.vtu
${CTEST_OUTPUT_DIR}/mmg2d_vtkvtu_ls)

# VTK .vtp with ls
ADD_TEST(NAME mmg2d_vtkvtp_ls
COMMAND ${EXECUT_MMG2D} -v 5 -ls 0.8
${MMG2D_CI_TESTS}/VtkInout/cercle_ls.vtp
${CTEST_OUTPUT_DIR}/mmg2d_vtkvtp_ls)

# VTK .vtk with ls and metric
ADD_TEST(NAME mmg2d_vtkvtk_ls_metric
  COMMAND ${EXECUT_MMG2D} -v 5 -ls 0.8
  ${MMG2D_CI_TESTS}/VtkInout/cercle_ls_metric.vtk
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtk_ls_metric)

# VTK .vtu with ls and metric
ADD_TEST(NAME mmg2d_vtkvtu_ls_metric
COMMAND ${EXECUT_MMG2D} -v 5 -ls 0.8
${MMG2D_CI_TESTS}/VtkInout/cercle_ls_metric.vtu
${CTEST_OUTPUT_DIR}/mmg2d_vtkvtu_ls_metric)

# VTK .vtp with ls and metric
ADD_TEST(NAME mmg2d_vtkvtp_ls_metric
  COMMAND ${EXECUT_MMG2D} -v 5 -ls 0.8
  ${MMG2D_CI_TESTS}/VtkInout/cercle_ls_metric.vtp
  ${CTEST_OUTPUT_DIR}/mmg2d_vtkvtp_ls_metric)

# VTK .vtk with metric and ls
ADD_TEST(NAME mmg2d_vtkvtk_metric_ls
COMMAND ${EXECUT_MMG2D} -v 5 -ls 0.8
${MMG2D_CI_TESTS}/VtkInout/cercle_metric_ls.vtk
${CTEST_OUTPUT_DIR}/mmg2d_vtkvtk_metric_ls)

IF ( (NOT VTK_FOUND) OR USE_VTK MATCHES OFF )
  SET(expr "VTK library not founded")
  SET_PROPERTY(TEST
    mmg2d_vtkvtk
    mmg2d_vtkvtp
    mmg2d_vtkvtu
    mmg2d_vtkvtk_iso
    mmg2d_vtkvtp_iso
    mmg2d_vtkvtu_iso
    mmg2d_vtkvtk_ani
    mmg2d_vtkvtp_ani
    mmg2d_vtkvtu_ani
    mmg2d_vtkvtk_ls
    mmg2d_vtkvtu_ls
    mmg2d_vtkvtp_ls
    mmg2d_vtkvtk_ls_metric
    mmg2d_vtkvtu_ls_metric
    mmg2d_vtkvtp_ls_metric
    mmg2d_vtkvtk_metric_ls
    PROPERTY PASS_REGULAR_EXPRESSION "${expr}")
ENDIF()

# Triangle output
#
# Respect the default Tetgen behaviour: saves only boundary edges in
# .edge file.
ADD_TEST(NAME mmg2d_Circle-triangle
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_Circle.o.node)

###############################################################################
#####
#####         Isotropic cases
#####
###############################################################################
ADD_TEST(NAME mmg2d_SquareIso
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/SquareIso/carretest
  -out ${CTEST_OUTPUT_DIR}/mmg2d_SquareIso-carretest.o.meshb)

ADD_TEST(NAME mmg2d_SquareIso_nonConstant
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/SquareIso/non-constant
  -out ${CTEST_OUTPUT_DIR}/mmg2d_non-constant.o.meshb)

ADD_TEST(NAME mmg2d_SquareIso_nonConstant2
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/SquareIso/non-constant-2
  -out ${CTEST_OUTPUT_DIR}/mmg2d_non-constant-2.o.meshb)

####### -nosurf option
ADD_TEST(NAME mmg2d_2squares
  COMMAND ${EXECUT_MMG2D} -3dMedit 2 -hmax 1 -nosurf -v 5
  ${MMG2D_CI_TESTS}/2squares/2squares
  -out ${CTEST_OUTPUT_DIR}/mmg2d_2squares.o.meshb)

####### -nsd
ADD_TEST(NAME mmg2d_2squares-nsd16
  COMMAND ${EXECUT_MMG2D} -3dMedit 2 -v 5 -nsd 16
  ${MMG2D_CI_TESTS}/2squares/2squares
  -out ${CTEST_OUTPUT_DIR}/mmg2d_2squares-nsd16.o.meshb)

####### orphan
ADD_TEST(NAME mmg2d_2squares-orphan
  COMMAND ${EXECUT_MMG2D} -3dMedit 2 -v 5 -nsd 10
  ${MMG2D_CI_TESTS}/2squares/2squares
  -out ${CTEST_OUTPUT_DIR}/mmg2d_2squares-nsd10.o.meshb)

####### -met option
ADD_TEST(NAME mmg2d_2squares-withMet
  COMMAND ${EXECUT_MMG2D} -3dMedit 2  -v 5
  ${MMG2D_CI_TESTS}/2squares/2squares -met ${MMG2D_CI_TESTS}/2squares/2s.sol
  -out ${CTEST_OUTPUT_DIR}/mmg2d_2squares-met.o.meshb)

####### -sol option
ADD_TEST(NAME mmg2d_2squares-withSol
  COMMAND ${EXECUT_MMG2D} -3dMedit 2  -v 5
  ${MMG2D_CI_TESTS}/2squares/2squares -sol ${MMG2D_CI_TESTS}/2squares/2s.sol
  -out ${CTEST_OUTPUT_DIR}/mmg2d_2squares-sol.o.meshb)

# -nreg
ADD_TEST(NAME mmg2d_nreg
  COMMAND ${EXECUT_MMG2D} -v 5 -nreg
  ${MMG2D_CI_TESTS}/SquareIso/carretest
  -out ${CTEST_OUTPUT_DIR}/mmg2d_nreg.o.meshb)

###############################################################################
#####
#####         Anisotropic cases
#####
###############################################################################
ADD_TEST(NAME mmg2d_SquareAniso
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/SquareAniso/adap1
  ${CTEST_OUTPUT_DIR}/mmg2d_SquareAniso-mmg2d_SquareAniso-adap1.o.meshb)

# optim
ADD_TEST(NAME mmg2d_Circle-optimAni
  COMMAND ${EXECUT_MMG2D} -v 5 -optim -A -sol 2
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_Circle-optimAni.o.mesh)

ADD_TEST(NAME mmg2d_Circle-hsizAni
  COMMAND ${EXECUT_MMG2D} -v 5 -hsiz 0.01 -A -sol 2
  ${MMG2D_CI_TESTS}/Circle/cercle
  -out ${CTEST_OUTPUT_DIR}/mmg2d_Circle-hsizAni.o.mesh)

# optim + ani + oprhan + unused point
ADD_TEST(NAME mmg2d_Disk-optimAni
  COMMAND ${EXECUT_MMG2D} -v 5 -optim -A -sol 2
  ${MMG2D_CI_TESTS}/Disk/disk-orphan
  -out ${CTEST_OUTPUT_DIR}/mmg2d_disk-optimAni.o.mesh)

# optim + iso + oprhan + unused point
ADD_TEST(NAME mmg2d_Disk-optim
  COMMAND ${EXECUT_MMG2D} -v 5 -optim -sol 2
  ${MMG2D_CI_TESTS}/Disk/disk-orphan
  -out ${CTEST_OUTPUT_DIR}/mmg2d_disk-optim.o.mesh)

###############################################################################
#####
#####         Mesh generation
#####
###############################################################################
ADD_TEST(NAME mmg2d_SquareGeneration
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/SquareGeneration/carretest
  ${CTEST_OUTPUT_DIR}/mmg2d_SquareGeneration-carretest.o.meshb)

ADD_TEST(NAME mmg2d_NacaGeneration
  COMMAND ${EXECUT_MMG2D} -v 5 -hausd 0.001
  ${MMG2D_CI_TESTS}/NacaGeneration/naca
  -out ${CTEST_OUTPUT_DIR}/mmg2d_NacaGeneration-naca.o.meshb)

ADD_TEST(NAME mmg2d_NacaGenerationAni
  COMMAND ${EXECUT_MMG2D} -v 5 -hausd 0.001 -A
  ${MMG2D_CI_TESTS}/NacaGeneration/naca
  -out ${CTEST_OUTPUT_DIR}/mmg2d_NacaGeneration-naca.o.meshb)

# optim
ADD_TEST(NAME mmg2d_NacaGeneration-optim
  COMMAND ${EXECUT_MMG2D} -v 5 -hausd 0.001 -optim
  ${MMG2D_CI_TESTS}/NacaGeneration/naca
  -out ${CTEST_OUTPUT_DIR}/mmg2d_NacaGeneration-optim.o.meshb)

# hsiz
ADD_TEST(NAME mmg2d_NacaGeneration-hsiz
  COMMAND ${EXECUT_MMG2D} -v 5 -hausd 0.001 -hsiz 0.01
  ${MMG2D_CI_TESTS}/NacaGeneration/naca
  -out ${CTEST_OUTPUT_DIR}/mmg2d_NacaGeneration-hsiz.o.meshb)

# hsiz + ani
ADD_TEST(NAME mmg2d_NacaGeneration-hsizAni
  COMMAND ${EXECUT_MMG2D} -v 5 -hausd 0.001 -hsiz 0.01 -A
  ${MMG2D_CI_TESTS}/NacaGeneration/naca
  -out ${CTEST_OUTPUT_DIR}/mmg2d_NacaGeneration-hsizAni.o.meshb)

# non convex test cases
ADD_TEST(NAME mmg2d_ACDCGeneration
  COMMAND ${EXECUT_MMG2D} -v 5 -d
  ${MMG2D_CI_TESTS}/ACDCGeneration/acdcBdy.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_ACDCGeneration.o.meshb)

# nsd option: keep only domain of ref 2
ADD_TEST(NAME mmg2d_ACDCGeneration-nsd2
  COMMAND ${EXECUT_MMG2D} -v 5 -nsd 2 -d
  ${MMG2D_CI_TESTS}/ACDCGeneration/acdcBdy.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_ACDCGeneration-nds2.o.meshb)

ADD_TEST(NAME mmg2d_GaronneGeneration
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/GaronneGeneration/garonneEdges.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_GaronneGeneration.o.meshb)

ADD_TEST(NAME mmg2d_GaronneGeneration2
  COMMAND ${EXECUT_MMG2D} -v 5
  ${MMG2D_CI_TESTS}/GaronneGeneration/garonne.mesh
  -out ${CTEST_OUTPUT_DIR}/mmg2d_GaronneGeneration2.o.meshb)

###############################################################################
#####
#####         Implicit domain discretization
#####
###############################################################################
ADD_TEST(NAME mmg2d_LSMultiMat_val
  COMMAND ${EXECUT_MMG2D} -val -v 5 -ls -hausd 0.001
  -met ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_multi-mat-val.o.meshb
  )

ADD_TEST(NAME mmg2d_OptLs_Bridge
  COMMAND ${EXECUT_MMG2D} -v 5 -ls
  -sol ${MMG2D_CI_TESTS}/OptLs_bridge/bridge.sol
  ${MMG2D_CI_TESTS}/OptLs_bridge/bridge
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLs_bridge.o.meshb
  )

# lssurf: discretization of boundaries only
ADD_TEST(NAME mmg2d_OptLsSurf_box
  COMMAND ${EXECUT_MMG2D} -v 5 -lssurf
  -sol ${MMG2D_CI_TESTS}/OptLsSurf_box/box.sol
  ${MMG2D_CI_TESTS}/OptLsSurf_box/box.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLsSurf_box.o.meshb
  )

# lssurf + multimat: discretization of boundaries only
ADD_TEST(NAME mmg2d_OptLsSurf_multiMat_box
  COMMAND ${EXECUT_MMG2D} -v 5 -lssurf
  -sol ${MMG2D_CI_TESTS}/OptLsSurf_box/box.sol
  ${MMG2D_CI_TESTS}/OptLsSurf_box/box_multiMat.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLsSurf_multiMat_box.o.meshb
  )

#multi-mat + opnbdy + non-manifold check
ADD_TEST(NAME mmg2d_LSMultiMat_nm
  COMMAND ${EXECUT_MMG2D} -v 5 -ls 3 -opnbdy -nr
  ${MMG2D_CI_TESTS}/LSMultiMat/2d-opn.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_2d-opn.o.meshb
  )

####### -nsd
ADD_TEST(NAME mmg2d_LSMultiMat-nsd22
  COMMAND ${EXECUT_MMG2D} -nsd 22 -v 5 -ls
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_multi-mat-nsd22.o.mesh
  )

#ADD_TEST(NAME mmg2d_LSMultiMat_default
#  COMMAND ${EXECUT_MMG2D} -val -v 5 -ls -hausd 0.001
#  -met ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-met.sol
#  -sol ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
#  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
#  ${CTEST_OUTPUT_DIR}/mmg2d_multi-mat-default.o.meshb
#  )
SET_PROPERTY(TEST mmg2d_LSMultiMat_val #mmg2d_LSMultiMat_default
  PROPERTY WILL_FAIL TRUE)

ADD_TEST(NAME mmg2d_LSDiscretization
  COMMAND ${EXECUT_MMG2D} -v 5 -ls
  ${MMG2D_CI_TESTS}/LSDiscretization/dom
  -out ${CTEST_OUTPUT_DIR}/mmg2d_LSDiscretization-dom.o.meshb)

ADD_TEST(NAME mmg2d_LSDiscretization2
  COMMAND ${EXECUT_MMG2D} -v 5 -ls
  ${MMG2D_CI_TESTS}/LSDiscretization/nacai
  -out ${CTEST_OUTPUT_DIR}/mmg2d_LSDiscretization2-nacai.o.meshb)

ADD_TEST(NAME mmg2d_LSMultiMat
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -hmin 0.005 -hmax 0.1 -hausd 0.001 -hgrad 1.3
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_LSMultiMat.o.meshb)

# non 0 ls
ADD_TEST(NAME mmg2d_LSMultiMat_nonzero
  COMMAND ${EXECUT_MMG2D} -v 5 -ls 0.01 -hausd 0.001
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_LSMultiMat-nonzero.o.meshb)

# ls + rmc
ADD_TEST(NAME mmg2d_OptLs_dom_withbub
  COMMAND ${EXECUT_MMG2D} -v 5 -ls
  ${MMG2D_CI_TESTS}/LSDiscretization/dom
  -sol ${MMG2D_CI_TESTS}/LSDiscretization/bub.sol
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLs_dom-withbub.o.meshb)

# ls + rmc + LSBaseReference
ADD_TEST(NAME mmg2d_OptLs_LSBaseReferences-rmc
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -rmc
  ${MMG2D_CI_TESTS}/LSBaseReferences/box
  -sol ${MMG2D_CI_TESTS}/LSBaseReferences/box.sol
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLs_LSBaseReferences-rmc.o.meshb)

ADD_TEST(NAME mmg2d_OptLs_LSBaseReferences-normc
  COMMAND ${EXECUT_MMG2D} -v 5 -ls
  ${MMG2D_CI_TESTS}/LSBaseReferences/box
  -sol ${MMG2D_CI_TESTS}/LSBaseReferences/box.sol
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLs_LSBaseReferences-normc.o.meshb)

# ls + rmc: max pile size bug
ADD_TEST(NAME mmg2d_OptLs_dom_rmcmaxpile
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -rmc
  ${MMG2D_CI_TESTS}/LSDiscretization/dom
  -sol ${MMG2D_CI_TESTS}/LSDiscretization/whole.sol
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLs_dom-rmcmaxpile.o.meshb)

ADD_TEST(NAME mmg2d_OptLs_dom_rembub
  COMMAND ${EXECUT_MMG2D} -v 5 -ls
  ${MMG2D_CI_TESTS}/LSDiscretization/dom
  -sol ${MMG2D_CI_TESTS}/LSDiscretization/bub.sol
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLs_dom-rembub.o.meshb -rmc)

ADD_TEST(NAME mmg2d_OptLs_dom_rembub2
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -rmc 0.1
  ${MMG2D_CI_TESTS}/LSDiscretization/dom
  -sol ${MMG2D_CI_TESTS}/LSDiscretization/bub.sol
  ${CTEST_OUTPUT_DIR}/mmg2d_OptLs_dom-rembub2.o.meshb)

add_test(
  NAME mmg2d_OptLs_isoref_defaut
  COMMAND ${EXECUT_MMG2D} -v 5 -ls ${MMG2D_CI_TESTS}/OptLs_isoref/2d-mesh.mesh
  -sol ${MMG2D_CI_TESTS}/OptLs_isoref/2d-mesh.sol
  ${CTEST_OUTPUT_DIR}/mmg2d_isoref.o.mesh
  )
add_test(
  NAME mmg2d_OptLs_isoref_5
  COMMAND ${EXECUT_MMG2D} -v 5 -isoref 5 -ls
  ${MMG2D_CI_TESTS}/OptLs_isoref/2d-mesh-isoref5.mesh
  -sol ${MMG2D_CI_TESTS}/OptLs_isoref/2d-mesh.sol
  ${CTEST_OUTPUT_DIR}/mmg2d_isoref5.o.mesh
  )

if (BASH)
  add_test(
    NAME mmg2d_optLs_isoref
    COMMAND ${BASH} -c "diff <(wc -wl ${CTEST_OUTPUT_DIR}/mmg2d_isoref.o.mesh  | awk '{print $1 $2}') <(wc -wl ${CTEST_OUTPUT_DIR}/mmg2d_isoref5.o.mesh | awk '{print $1 $2}')"
    )
endif()

# ls discretisation + parameter file
ADD_TEST(NAME mmg2d_ParsOpName
  COMMAND ${EXECUT_MMG2D} -v 5 -ls
  -f ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-refs.mmg2d
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_ParsOpName.o.meshb)

SET(parsopName "multi-mat-refs.mmg2d OPENED")
SET_PROPERTY(TEST mmg2d_ParsOpName
  PROPERTY PASS_REGULAR_EXPRESSION "${parsopName}")

# ls discretisation + wrong name of parameter file
ADD_TEST(NAME mmg2d_ParsOpName_wrongFile
  COMMAND ${EXECUT_MMG2D} -v 5 -ls
  -f ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-false.mmg2d
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_ParsOpName_wrongFile.o.meshb)

SET(parsopNameWrong "multi-mat-false.mmg2d file NOT FOUND.")
SET_PROPERTY(TEST mmg2d_ParsOpName_wrongFile
  PROPERTY PASS_REGULAR_EXPRESSION "${parsopNameWrong}")

# ls discretisation + no name of parameter file
ADD_TEST(NAME mmg2d_ParsOpName_NoFileName
  COMMAND ${EXECUT_MMG2D} -v 5 -f -ls
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat.mesh
  ${CTEST_OUTPUT_DIR}/mmg2d_ParsOpName_NoFileName.o.meshb)

SET(parsopNameNo "Missing filename for -f")
SET_PROPERTY(TEST mmg2d_ParsOpName_NoFileName
  PROPERTY PASS_REGULAR_EXPRESSION "${parsopNameNo}")

  # ls discretisation + optim option
ADD_TEST(NAME mmg2d_LSMultiMat_optim
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -optim -hausd 0.001
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_LSMultiMat-optim.o.meshb)

# ls discretisation + optim + aniso option
ADD_TEST(NAME mmg2d_LSMultiMat_optimAni
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -optim -A -hausd 0.001
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_LSMultiMat-optimAni.o.meshb)

# ls discretisation + hsiz option
ADD_TEST(NAME mmg2d_LSMultiMat_hsiz
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -hsiz 0.05 -hausd 0.001
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_LSMultiMat-hsiz.o.meshb)

# ls discretisation + hsiz Ani option
ADD_TEST(NAME mmg2d_LSMultiMat_hsizAni
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -hsiz 0.05 -A -hausd 0.001
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_LSMultiMat-hsizAni.o.meshb)

# ls discretisation + metric
ADD_TEST(NAME mmg2d_LSMultiMat_withMet
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -hausd 0.001
  -met ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_LSMultiMat-withMet.o.meshb)

# ls discretisation + metric + ls
ADD_TEST(NAME mmg2d_LSMultiMat_withMetAndLs
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -hausd 0.001
  -met ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-met.sol
  -sol ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat-sol.sol
  ${MMG2D_CI_TESTS}/LSMultiMat/multi-mat
  ${CTEST_OUTPUT_DIR}/mmg2d_LSMultiMat-withMetAndLs.o.meshb)

# ls discretisation + xreg
ADD_TEST(NAME mmg2d_CoorRegularization_apple
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -xreg
  ${MMG2D_CI_TESTS}/CoorRegularization_apple/apple
  -out ${CTEST_OUTPUT_DIR}/CoorRegularization_apple.o.meshb)

# ls discretisation + xreg + nr
ADD_TEST(NAME mmg2d_CoorRegularization_appleNR
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -xreg -nr
  ${MMG2D_CI_TESTS}/CoorRegularization_apple/apple
  -out ${CTEST_OUTPUT_DIR}/CoorRegularization_appleNR.o.meshb)

# ls discretisation + xreg + nr + check of negative areas
ADD_TEST(NAME mmg2d_CoorRegularizationNegativeArea
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -xreg -hmax 0.1
  ${MMG2D_CI_TESTS}/CoorRegularizationNegativeArea/CoorRegularizationNegativeArea
  -out ${CTEST_OUTPUT_DIR}/CoorRegularizationNegativeArea.o.meshb)

# ls discretisation + xreg + choice of value for xreg
ADD_TEST(NAME mmg2d_CoorRegularization_apple_value
  COMMAND ${EXECUT_MMG2D} -v 5 -ls -xreg 0.9
  ${MMG2D_CI_TESTS}/CoorRegularization_apple/apple
  -out ${CTEST_OUTPUT_DIR}/CoorRegularization_apple_value.o.meshb)

###############################################################################
#####
#####         Check Lagrangian motion option
#####
###############################################################################
#####
IF ( ELAS_FOUND AND NOT USE_ELAS MATCHES OFF )
  ADD_TEST(NAME mmg2d_LagMotion0_circle
    COMMAND ${EXECUT_MMG2D} -v 5  -lag 0
    -in ${MMG2D_CI_TESTS}/LagMotion_circle/circle
    -out ${CTEST_OUTPUT_DIR}/mmg2d_LagMotion0_circle-circle.o.meshb
    )
  ADD_TEST(NAME mmg2d_LagMotion1_circle
    COMMAND ${EXECUT_MMG2D} -v 5  -lag 1
    -in ${MMG2D_CI_TESTS}/LagMotion_circle/circle
    -out ${CTEST_OUTPUT_DIR}/mmg2d_LagMotion1_circle-circle.o.meshb
    )
  ADD_TEST(NAME mmg2d_LagMotion2_circle
    COMMAND ${EXECUT_MMG2D} -v 5  -lag 2 -d
    -in ${MMG2D_CI_TESTS}/LagMotion_circle/circle
    -out ${CTEST_OUTPUT_DIR}/mmg2d_LagMotion2_circle-circle.o.meshb
    )

  # nsd
  ADD_TEST(NAME mmg2d_LagMotion2_circle-nsd3
    COMMAND ${EXECUT_MMG2D} -v 5  -lag 2 -nsd 3
    -in ${MMG2D_CI_TESTS}/LagMotion_circle/circle
    -out ${CTEST_OUTPUT_DIR}/mmg2d_LagMotion2_circle-nsd3.o.mesh
    )

  IF (${MMG5_INT} MATCHES int64_t )
    SET(passElasRegex "## Error: MMG2D_velextLS: impossible to call elasticity library with int64 integers")
    SET_PROPERTY(TEST mmg2d_LagMotion0_circle mmg2d_LagMotion1_circle mmg2d_LagMotion2_circle mmg2d_LagMotion2_circle-nsd3
      PROPERTY PASS_REGULAR_EXPRESSION "${passElasRegex}")
  ENDIF()

ENDIF()

###############################################################################
#####
#####         Check snapping (prevision of non-manifold situations)
#####
###############################################################################
#####
SET(nmRegex "unsnap at least 1 point")

ADD_TEST(NAME mmg2d_LSSnapval_manifold1
  COMMAND ${EXECUT_MMG2D} -v 5  -ls
  -in ${MMG2D_CI_TESTS}/LSSnapval/8elts1.mesh
  -sol ${MMG2D_CI_TESTS}/LSSnapval/manifold.sol
  -out ${CTEST_OUTPUT_DIR}/mmg2d_LSSnapval_manifold1.o.mesh
  )

ADD_TEST(NAME mmg2d_LSSnapval_manifold2
  COMMAND ${EXECUT_MMG2D} -v 5  -ls
  -in ${MMG2D_CI_TESTS}/LSSnapval/8elts2.mesh
  -sol ${MMG2D_CI_TESTS}/LSSnapval/manifold.sol
  -out ${CTEST_OUTPUT_DIR}/mmg2d_LSSnapval_manifold2.o.mesh
  )

SET_PROPERTY(TEST mmg2d_LSSnapval_manifold1 mmg2d_LSSnapval_manifold2
  PROPERTY FAIL_REGULAR_EXPRESSION "${nmRegex}")

ADD_TEST(NAME mmg2d_LSSnapval_non-manifold1
  COMMAND ${EXECUT_MMG2D} -v 5  -ls
  -in ${MMG2D_CI_TESTS}/LSSnapval/8elts1.mesh
  -sol ${MMG2D_CI_TESTS}/LSSnapval/8elts1-nm.sol
  -out ${CTEST_OUTPUT_DIR}/mmg2d_LSSnapval_non-manifold1.o.mesh
  )

ADD_TEST(NAME mmg2d_LSSnapval_non-manifold2
  COMMAND ${EXECUT_MMG2D} -v 5  -ls
  -in ${MMG2D_CI_TESTS}/LSSnapval/8elts2.mesh
  -sol ${MMG2D_CI_TESTS}/LSSnapval/8elts2-nm.sol
  -out ${CTEST_OUTPUT_DIR}/mmg2d_LSSnapval_non-manifold2.o.mesh
  )
SET_PROPERTY(TEST mmg2d_LSSnapval_non-manifold1 mmg2d_LSSnapval_non-manifold2
  PROPERTY PASS_REGULAR_EXPRESSION "${nmRegex}")
