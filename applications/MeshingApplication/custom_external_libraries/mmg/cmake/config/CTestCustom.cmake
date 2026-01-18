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

SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 5000)
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 2000)

SET ( CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
  ".*: warning: array subscript has type 'char'.*")
SET ( CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
  ".*: warning: array subscript has type ‘char’.*")
SET (CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 0 )
SET (CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 0 )
