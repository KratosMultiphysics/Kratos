/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

To download the library MMG go to:

http://www.mmgtools.org/mmg-remesher-downloads

or:

git clone https://github.com/MmgTools/mmg.git

Add the following to the main Kratos configure.sh: 

-DINCLUDE_MMG=ON                                                                  \
-DMMG_INCLUDE_DIR="kratos_dir/applications/MeshingApplication/custom_external_libraries/mmg/include/"\
-DMMG2D_INCLUDE_DIR="kratos_dir/applications/MeshingApplication/custom_external_libraries/mmg/include/mmg/mmg2d/"\
-DMMG3D_INCLUDE_DIR="kratos_dir/applications/MeshingApplication/custom_external_libraries/mmg/include/mmg/mmg3d/"\
-DMMGS_INCLUDE_DIR="kratos_dir/applications/MeshingApplication/custom_external_libraries/mmg/include/mmg/mmgs/"\
-DMMG_LIBRARY="kratos_dir/applications/MeshingApplication/custom_external_libraries/mmg/lib/libmmg.a"     \
-DMMG2D_LIBRARY="kratos_dir/applications/MeshingApplication/custom_external_libraries/mmg/lib/libmmg2d.a" \
-DMMG3D_LIBRARY="kratos_dir/applications/MeshingApplication/custom_external_libraries/mmg/lib/libmmg3d.a" \
-DMMGS_LIBRARY="kratos_dir/applications/MeshingApplication/custom_external_libraries/mmg/lib/libmmgs.a"   \

