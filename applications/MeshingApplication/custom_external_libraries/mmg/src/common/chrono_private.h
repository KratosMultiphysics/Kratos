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

#ifndef CHRONO_H
#define CHRONO_H

#include <time.h>
#include "mmgcommon_private.h"

#ifndef MMG_POSIX
#include <windows.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include "mmg_core_export_private.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef  ON
#define  RESET  0
#define  ON     1
#define  OFF    2
#endif

#define  BIG      1e6
#define  BIG1     1e-6
#define  TIMEMAX  12

/**
 * \struct mytime
 * \brief Chrono object.
 *
 * mytime allow storage of chronograph informations.
 *
 */
typedef struct mytime {
  double  gini,gend,gdif,uini,uend,udif,sini,send,sdif;
#ifdef MMG_POSIX
  struct  timeval rutim;
  struct  rusage  ru;
#else
  HANDLE        thisProcess;
  FILETIME      ftIni, ftEnd, ftSys, ftUser;
  SYSTEMTIME    stSys, stUser;
  LARGE_INTEGER frequency;
  LARGE_INTEGER rutim;
#endif
  int     call;
} mytime;


/* prototypes */
LIBMMG_CORE_EXPORT void   chrono(int cmode,mytime *ptt);
LIBMMG_CORE_EXPORT void   tminit(mytime *t,int maxtim);
LIBMMG_CORE_EXPORT void   printim(double ,char *);

#ifdef __cplusplus
}
#endif

#endif
