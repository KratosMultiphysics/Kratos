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

############################################################################################################
# AddCompileLinkFlags.cmake
# (copy from http://binglongx.com/2013/06/30/add-arbitrary-compilelink-flags-in-cmake/)
############################################################################################################

############################################################################################################
# Append str to a string property of a target.
# target:      string: target name.
# property:            name of target’s property. e.g: COMPILE_FLAGS, or LINK_FLAGS
# str:         string: string to be appended to the property
macro(my_append_target_property target property str)
  get_target_property(current_property ${target} ${property})
  if(NOT current_property) # property non-existent or empty
      set_target_properties(${target} PROPERTIES ${property} ${str})
  else()
      set_target_properties(${target} PROPERTIES ${property} "${current_property} ${str}")
  endif()
endmacro(my_append_target_property)

############################################################################################################
# Add/append compile flags to a target.
# target: string: target name.
# flags : string: compile flags to be appended
macro(my_add_compile_flags target flags)
  my_append_target_property(${target} COMPILE_FLAGS ${flags})
endmacro(my_add_compile_flags)

############################################################################################################
# Add/append link flags to a target.
# target: string: target name.
# flags : string: link flags to be appended
macro(my_add_link_flags target flags)
  my_append_target_property(${target} LINK_FLAGS ${flags})
endmacro(my_add_link_flags)
