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
#####         Aniso test case
#####
##############################################################################
#####
ADD_TEST(NAME TorusholesAni_${EXECUT_MMGS}_2
 COMMAND ${EXECUT_MMGS} -v 5 -hgrad 1.15
 ${MMG_CI_TESTS}/TorusholesAni/torusholes.o.meshb
 -out ${MMG_CI_TESTS}/TorusholesAni/torusholes.o.o.meshb)

SET_TESTS_PROPERTIES(TorusholesAni_${EXECUT_MMGS}_2 PROPERTIES DEPENDS TorusholesAni_${EXECUT_MMGS})

ADD_TEST(NAME TorusholesAni_chocCyl_${EXECUT_MMGS}_2
 COMMAND ${EXECUT_MMGS} -v 5 -hgrad 1.15
 ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny.o.meshb
 -out ${MMG_CI_TESTS}/TorusholesAni_chocCyl/torusholesTiny.o.o.meshb)

SET_TESTS_PROPERTIES(TorusholesAni_chocCyl_${EXECUT_MMGS}_2 PROPERTIES DEPENDS TorusholesAni_chocCyl_${EXECUT_MMGS})
