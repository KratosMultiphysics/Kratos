/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
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

/**
 * \file mmgs/roots_s.c
 * \brief Functions to compute roots of degree 2 and 3 polynomial.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgs.h"
#define EPSRO     1.e-3

#ifdef GNU
/* Returns the 2 complex roots of a degree 2 polynomial with real coefficients a[2]T^2 + ... + a[0]
   By convention, the real roots are stored first (same thing for multiple roots) :
   return value = number of roots, counted with multiplicity */
int rootDeg2(double complex a[3], double complex r[2]){
  double complex Delta,delta,r1,r2;

  if( cabs(a[2])<_MMG5_EPSD ) {
    if( cabs(a[1])<_MMG5_EPSD ) {
      r[0] = r[1] = 0.0;
      return(0);
    }

    else{
      r[0] = r[1] = -a[0]/a[1];
      return(1);
    }
  }

  Delta = a[1]*a[1] - 4.0*a[2]*a[0];
  /* ONE square root of Delta */
  delta = cpow(Delta,0.5);
  r1 = 0.5/a[2]*(-a[1] - delta);
  r2 = 0.5/a[2]*(-a[1] + delta);

  if( fabs(cimag(r1)) < EPSRO ){
    r[0] = creal(r1);
    r[1] = ( fabs(cimag(r2)) < EPSRO ) ? creal(r2) : r2;
  }
  else {
    r[0] = ( fabs(cimag(r2)) < EPSRO ) ? creal(r2) : r2;
    r[1] = r1;
  }

  return(2);
}

#else
/* Returns the 2 complex roots of a degree 2 polynomial with real coefficients a[2]T^2 + ... + a[0]
   By convention, the real roots are stored first (same thing for multiple roots) :
   return value = number of roots, counted with multiplicity */
int rootDeg2(DOUBLE_COMPLEX a[3], DOUBLE_COMPLEX r[2]){
  DOUBLE_COMPLEX Delta,delta,r1,r2;
  double         real_tmp, imag_tmp, denom,real_r,imag_r;

  if( cabs(a[2])<_MMG5_EPSD ) {
    if( cabs(a[1])<_MMG5_EPSD ) {
      r[0] = r[1] = _DCOMPLEX_(0.0,0.0);
      return(0);
    }

    else{
	  real_tmp = creal(a[0])*creal(a[1]) + cimag(a[0])*cimag(a[1]);
	  imag_tmp = cimag(a[0])*creal(a[1]) - creal(a[0])*cimag(a[1]);
	  denom    = creal(a[1])*creal(a[1]) + cimag(a[1])*cimag(a[1]);

      r[0] = r[1] = _DCOMPLEX_(-real_tmp/denom,-imag_tmp/denom);
      return(1);
    }
  }

  real_tmp = creal(a[1])*creal(a[1]) - cimag(a[1])*cimag(a[1])
	  -4.0*(creal(a[2])*creal(a[0]) - cimag(a[2])*cimag(a[0]));
  imag_tmp = creal(a[1])*cimag(a[1]) + cimag(a[1])*creal(a[1])
	  - 4.0*(creal(a[2])*cimag(a[0]) + cimag(a[2])*creal(a[0]));
  Delta = _DCOMPLEX_(real_tmp, imag_tmp);

  /* ONE square root of Delta */
  delta = cpow(Delta,_DCOMPLEX_(0.5,0.0));
  a[2] = _DCOMPLEX_(2. * creal(a[2]), 2. * cimag(a[2]));

  real_tmp = -creal(a[1]) - creal(delta);
  imag_tmp = -cimag(a[1]) - cimag(delta);
  real_r = real_tmp*creal(a[2]) + imag_tmp*cimag(a[2]);
  imag_r = imag_tmp*creal(a[2]) - real_tmp*cimag(a[2]);
  denom = creal(a[2])*creal(a[2]) + cimag(a[2])*cimag(a[2]);
  r1 = _DCOMPLEX_(real_r / denom, imag_r / denom);
  //r1 = 0.5/a[2]*(-a[1] - delta);

  real_tmp = -creal(a[1]) + creal(delta);
  imag_tmp = -cimag(a[1]) + cimag(delta);
  real_r = real_tmp*creal(a[2]) + imag_tmp*cimag(a[2]);
  imag_r = imag_tmp*creal(a[2]) - real_tmp*cimag(a[2]);
  denom = creal(a[2])*creal(a[2]) + cimag(a[2])*cimag(a[2]);
  r2 = _DCOMPLEX_(real_r / denom, imag_r / denom);
 // r2 = 0.5/a[2]*(-a[1] + delta);

  if( fabs(cimag(r1)) < EPSRO ){
    r[0] = _DCOMPLEX_(creal(r1),0.0);
    r[1] = ( fabs(cimag(r2)) < EPSRO ) ? _DCOMPLEX_(creal(r2),0.0) : r2;
  }
  else {
    r[0] = ( fabs(cimag(r2)) < EPSRO ) ? _DCOMPLEX_(creal(r2),0.0) : r2;
    r[1] = r1;
  }

  return(2);
}
#endif
