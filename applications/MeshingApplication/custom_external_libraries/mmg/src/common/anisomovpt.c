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

/**
 * \file common/anisomovpt.c
 * \brief Functions to move a point in the mesh (with anisotropic metric).
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "mmgcommon_private.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param pt pointer toward the tria on which we integrate.
 * \param p0 pointer toward the point that we want to move.
 * \param pb bezier patch of the triangle.
 * \param r rotation matrix that sends the normal at point \a p0 to e_z.
 * \param gv centre of mass that we want to update using the computed element
 * weight.
 *
 * \return 0 if fail, 1 otherwise.
 *
 * Compute integral of sqrt(T^J(xi)  M(P(xi)) J(xi)) * P(xi) over the triangle.
 *
 */
int MMG5_elementWeight(MMG5_pMesh mesh,MMG5_pSol met, MMG5_pTria pt,
                        MMG5_pPoint p0, MMG5_Bezier *pb,double r[3][3],
                        double gv[2])
{
  MMG5_pPoint    p1,p2;
  double         Jacsigma[3][2],Jactmp[3][2],m[6],mo[6],density,to[3],no[3],ll;
  double         dens[3],*n1,*n2,ps1,ps2,intpt[2],ux,uy,uz;
  int8_t         i0,i1,i2,j,nullDens;
  static int8_t  mmgErr=0;

  i0 = 0;
  i1 = 1;
  i2 = 2;

  nullDens = 0;
  for (j=0; j<3; j++) {
    /* Set jacobian matrix of parametric Bezier patch, for each quadrature point */
    if ( j == 0 ) {  //w,u = 1/2, v = 0
      Jacsigma[0][0] = 0.75*(pb->b[7][0] - pb->b[0][0]) + 0.75*(pb->b[1][0] - pb->b[8][0]) + 1.5*(pb->b[8][0] - pb->b[7][0]);
      Jacsigma[1][0] = 0.75*(pb->b[7][1] - pb->b[0][1]) + 0.75*(pb->b[1][1] - pb->b[8][1]) + 1.5*(pb->b[8][1] - pb->b[7][1]);
      Jacsigma[2][0] = 0.75*(pb->b[7][2] - pb->b[0][2]) + 0.75*(pb->b[1][2] - pb->b[8][2]) + 1.5*(pb->b[8][2] - pb->b[7][2]);

      Jacsigma[0][1] = 0.75*(pb->b[6][0] - pb->b[0][0]) + 0.75*(pb->b[3][0] - pb->b[8][0]) + 1.5*(pb->b[9][0] - pb->b[7][0]);
      Jacsigma[1][1] = 0.75*(pb->b[6][1] - pb->b[0][1]) + 0.75*(pb->b[3][1] - pb->b[8][1]) + 1.5*(pb->b[9][1] - pb->b[7][1]);
      Jacsigma[2][1] = 0.75*(pb->b[6][2] - pb->b[0][2]) + 0.75*(pb->b[3][2] - pb->b[8][2]) + 1.5*(pb->b[9][2] - pb->b[7][2]);
    }
    else if( j == 1 ) { //u,v = 1/2, w = 0
      Jacsigma[0][0] = 0.75*(pb->b[1][0] - pb->b[8][0]) + 0.75*(pb->b[4][0] - pb->b[5][0]) + 1.5*(pb->b[3][0] - pb->b[9][0]);
      Jacsigma[1][0] = 0.75*(pb->b[1][1] - pb->b[8][1]) + 0.75*(pb->b[4][1] - pb->b[5][1]) + 1.5*(pb->b[3][1] - pb->b[9][1]);
      Jacsigma[2][0] = 0.75*(pb->b[1][2] - pb->b[8][2]) + 0.75*(pb->b[4][2] - pb->b[5][2]) + 1.5*(pb->b[3][2] - pb->b[9][2]);

      Jacsigma[0][1] = 0.75*(pb->b[3][0] - pb->b[8][0]) + 0.75*(pb->b[2][0] - pb->b[5][0]) + 1.5*(pb->b[4][0] - pb->b[9][0]);
      Jacsigma[1][1] = 0.75*(pb->b[3][1] - pb->b[8][1]) + 0.75*(pb->b[2][1] - pb->b[5][1]) + 1.5*(pb->b[4][1] - pb->b[9][1]);
      Jacsigma[2][1] = 0.75*(pb->b[3][2] - pb->b[8][2]) + 0.75*(pb->b[2][2] - pb->b[5][2]) + 1.5*(pb->b[4][2] - pb->b[9][2]);
    }
    else { //w,v = 1/2, u= 0
      Jacsigma[0][0] = 0.75*(pb->b[7][0] - pb->b[0][0]) + 0.75*(pb->b[4][0] - pb->b[5][0]) + 1.5*(pb->b[9][0] - pb->b[6][0]);
      Jacsigma[1][0] = 0.75*(pb->b[7][1] - pb->b[0][1]) + 0.75*(pb->b[4][1] - pb->b[5][1]) + 1.5*(pb->b[9][1] - pb->b[6][1]);
      Jacsigma[2][0] = 0.75*(pb->b[7][2] - pb->b[0][2]) + 0.75*(pb->b[4][2] - pb->b[5][2]) + 1.5*(pb->b[9][2] - pb->b[6][2]);

      Jacsigma[0][1] = 0.75*(pb->b[6][0] - pb->b[0][0]) + 0.75*(pb->b[2][0] - pb->b[5][0]) + 1.5*(pb->b[5][0] - pb->b[6][0]);
      Jacsigma[1][1] = 0.75*(pb->b[6][1] - pb->b[0][1]) + 0.75*(pb->b[2][1] - pb->b[5][1]) + 1.5*(pb->b[5][1] - pb->b[6][1]);
      Jacsigma[2][1] = 0.75*(pb->b[6][2] - pb->b[0][2]) + 0.75*(pb->b[2][2] - pb->b[5][2]) + 1.5*(pb->b[5][2] - pb->b[6][2]);
    }

    /* Take metric at control point */
    if ( !MG_GEO_OR_NOM(pt->tag[i2]) ) {
      int ier = MMG5_interpreg_ani(mesh,met,pt,i2,0.5,m);
      if ( !ier ) {
        if ( mesh->info.ddebug ) {
          fprintf(stdout,"  ## Warning:%s:%d: unable to move point (interpreg_ani failure).\n",
                  __func__,__LINE__);
        }
        return 0;
      }
    }
    else {
      if ( !MMG5_nortri(mesh,pt,no) )  return 0;
      if ( !MMG5_intridmet(mesh,met,pt->v[i0],pt->v[i1],0.5,no,mo) )  return 0;

      p1 = &mesh->point[pt->v[i0]];
      p2 = &mesh->point[pt->v[i1]];

      to[0] = p2->c[0] - p1->c[0];
      to[1] = p2->c[1] - p1->c[1];
      to[2] = p2->c[2] - p1->c[2];

      ll = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
      if ( ll < MMG5_EPSD )  return 0;
      ll = 1.0 / sqrt(ll);
      to[0] *= ll;
      to[1] *= ll;
      to[2] *= ll;

      if ( ( MG_SIN(p1->tag) || (p1->tag & MG_NOM) )
           && ( MG_SIN(p2->tag) || (p2->tag & MG_NOM) ) ) {
             if ( !MMG5_buildridmetfic(mesh,to,no,mo[0],mo[0],mo[0],m) )
               return 0;
      }
      else if ( !(MG_SIN(p1->tag) || (p1->tag & MG_NOM)) ) {
        n1 = &mesh->xpoint[p1->xp].n1[0];
        n2 = &mesh->xpoint[p1->xp].n2[0];
        ps1 = n1[0]*no[0] + n1[1]*no[1] + n1[2]*no[2];
        ps2 = n2[0]*no[0] + n2[1]*no[1] + n2[2]*no[2];
        if ( fabs(ps1) > fabs(ps2) ) {
          if ( !MMG5_buildridmetfic(mesh,to,no,mo[0],mo[1],mo[3],m) )  return 0;
        }
        else {
          if ( !MMG5_buildridmetfic(mesh,to,no,mo[0],mo[2],mo[4],m) )  return 0;
        }
      }
      else {
        assert(!(MG_SIN(p2->tag)||(p2->tag & MG_NOM)));
        n1 = &mesh->xpoint[p2->xp].n1[0];
        n2 = &mesh->xpoint[p2->xp].n2[0];
        ps1 = n1[0]*no[0] + n1[1]*no[1] + n1[2]*no[2];
        ps2 = n2[0]*no[0] + n2[1]*no[1] + n2[2]*no[2];
        if ( fabs(ps1) > fabs(ps2) ) {
          if ( !MMG5_buildridmetfic(mesh,to,no,mo[0],mo[1],mo[3],m) )  return 0;
        }
        else {
          if ( !MMG5_buildridmetfic(mesh,to,no,mo[0],mo[2],mo[4],m) )  return 0;
        }
      }
    }

    /* Compute density matrix {^t}Jacsigma * M * Jacsigma */
    Jactmp[0][0] = m[0]*Jacsigma[0][0] + m[1]*Jacsigma[1][0] + m[2]*Jacsigma[2][0];
    Jactmp[1][0] = m[1]*Jacsigma[0][0] + m[3]*Jacsigma[1][0] + m[4]*Jacsigma[2][0];
    Jactmp[2][0] = m[2]*Jacsigma[0][0] + m[4]*Jacsigma[1][0] + m[5]*Jacsigma[2][0];

    Jactmp[0][1] = m[0]*Jacsigma[0][1] + m[1]*Jacsigma[1][1] + m[2]*Jacsigma[2][1];
    Jactmp[1][1] = m[1]*Jacsigma[0][1] + m[3]*Jacsigma[1][1] + m[4]*Jacsigma[2][1];
    Jactmp[2][1] = m[2]*Jacsigma[0][1] + m[4]*Jacsigma[1][1] + m[5]*Jacsigma[2][1];

    dens[0] = Jacsigma[0][0]*Jactmp[0][0] + Jacsigma[1][0]*Jactmp[1][0] + Jacsigma[2][0]*Jactmp[2][0];
    dens[1] = Jacsigma[0][0]*Jactmp[0][1] + Jacsigma[1][0]*Jactmp[1][1] + Jacsigma[2][0]*Jactmp[2][1];
    dens[2] = Jacsigma[0][1]*Jactmp[0][1] + Jacsigma[1][1]*Jactmp[1][1] + Jacsigma[2][1]*Jactmp[2][1];

    density = dens[0]*dens[2] - dens[1]*dens[1];
    if ( density <= MMG5_EPSD2 ) {
#ifndef NDEBUG
      if ( !mmgErr ) {
        fprintf(stderr,"\n  ## Warning: %s: at least 1 negative or null density.\n",
                __func__);
        mmgErr = 1;
      }
#endif
      ++nullDens;
      continue;
    }

    density = sqrt(density);
    /* Coordinates of the integration point */
    p1 = &mesh->point[pt->v[i0]];
    p2 = &mesh->point[pt->v[i1]];

    ux = 0.5*(p1->c[0] + p2->c[0]) - p0->c[0];
    uy = 0.5*(p1->c[1] + p2->c[1]) - p0->c[1];
    uz = 0.5*(p1->c[2] + p2->c[2]) - p0->c[2];

    intpt[0] =  r[0][0]*ux + r[0][1]*uy + r[0][2]*uz;
    intpt[1] =  r[1][0]*ux + r[1][1]*uy + r[1][2]*uz;

    gv[0] += density*MMG5_ATHIRD*intpt[0];
    gv[1] += density*MMG5_ATHIRD*intpt[1];

    i0 = MMG5_inxt2[i0];
    i1 = MMG5_inxt2[i1];
    i2 = MMG5_inxt2[i2];
  }

  if ( nullDens==3 ) return 0;

  return 1;
}
