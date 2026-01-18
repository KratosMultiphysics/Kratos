/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
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

#ifndef MMGEIGENV_H
#define MMGEIGENV_H

#ifdef __cplusplus
extern "C" {
#endif

#define MMG5_EPSD      1.e-30
#define MMG5_EPS       1.e-06

int MMG5_eigenv3d(int symmat,double *mat,double lambda[3],double v[3][3]);
int MMG5_eigenv2d(int symmat,double *mat,double lambda[2],double v[2][2]);
int MMG5_eigen2(double *mm,double *lambda,double vp[2][2]);
extern int MMG5_eigensym(double m[3],double lambda[2],double vp[2][2]);

#ifdef __cplusplus
}
#endif

#endif
