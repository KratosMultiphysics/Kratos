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

SET(CTEST_PROJECT_NAME "Mmg3d5")
SET(UPDATE_COMMAND "${GITCOMMAND}")

SET(CTEST_CONTINUOUS_DURATION 600) # duration of Continuous dashboard
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE TRUE)
SET(CTEST_CONTINUOUS_MINIMUM_INTERVAL 60) # Frequency between two check
SET(CTEST_NIGHTLY_START_TIME "01:00AM")

SET(CTEST_DROP_METHOD "http")
SET(CTEST_DROP_SITE "localhost:2000")
#SET(CTEST_DROP_SITE "vulcain.bordeaux.inria.fr")
SET(CTEST_DROP_LOCATION "/CDash/submit.php?project=Mmg3d5")
SET(CTEST_DROP_SITE_CDASH TRUE)
