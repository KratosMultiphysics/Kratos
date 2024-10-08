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
 * \file mmg3d/tools_3d.c
 * \brief Various algorithmic and algebraic tools.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"

extern int8_t ddb;

/** Return 1 if reference ref is in the br table, 0 otherwise */
int MMG5_isbr(MMG5_pMesh mesh,MMG5_int ref) {
  int k;

  for(k=0; k<mesh->info.nbr; k++)
    if ( ref == mesh->info.br[k] ) return(1);

  return(0);
}

/** Compute normal to face iface of tetra k, exterior to tetra k */
int MMG5_norface(MMG5_pMesh mesh,MMG5_int k,int iface,double n[3]) {
  MMG5_pTetra     pt;

  pt = &mesh->tetra[k];

  return MMG5_norpts(mesh,
                      pt->v[MMG5_idir[iface][0]],
                      pt->v[MMG5_idir[iface][1]],
                      pt->v[MMG5_idir[iface][2]],n);
}

/** If need be, invert the travelling sense of surfacic ball so that it is travelled in
    the direct sense with respect to direction n anchored at point ip (ip = global num.):
    return 2 = orientation reversed, 1 otherwise */
inline int MMG5_directsurfball(MMG5_pMesh mesh, MMG5_int ip, MMG5_int *list, int ilist, double n[3]){
    int      j;
    MMG5_int iel,aux;
    double   nt[3],ps;
    uint8_t  iface;

    iel   = list[0] / 4;
    iface = list[0] % 4;

    if ( !MMG5_norface(mesh,iel,iface,nt) ) return 0;
    ps = nt[0]*n[0] +  nt[1]*n[1] +  nt[2]*n[2];
    if ( ps > 0.0 )  return 1;

    for (j=1; j<=(ilist-1)/2; j++) {
        aux = list[j];
        list[j] = list[ilist -j];
        list[ilist -j] = aux;
    }

    return 2;
}

/** If need be, reorder the surfacic ball of point ip, so that its first element has
    edge (p,q) (nump,q = global num) as edge MMG5_iprv2[ip] of face iface.
    return 2 = orientation reversed, 1 otherwise */
int MMG5_startedgsurfball(MMG5_pMesh mesh,MMG5_int nump,MMG5_int numq,MMG5_int *list,int ilist) {
    MMG5_pTetra pt;
    int         l;
    MMG5_int    tmp,iel;
    uint8_t     iface,ip,ipt;

    iel = list[0]/4;
    iface = list[0]%4;
    pt = &mesh->tetra[iel];

    for(ip=0;ip<4;ip++){
        if(pt->v[ip] == nump) break;
    }
    assert(ip<4);

    ipt = MMG5_idirinv[iface][ip]; // index of ip in face iface
    ipt = MMG5_inxt2[ipt];         // next index in this face
    ipt = MMG5_idir[iface][ipt];  // index of this point in local num of tetra

    if(pt->v[ipt] == numq) return 1;

    else{
        ipt = MMG5_idir[iface][MMG5_iprv2[MMG5_idirinv[iface][ip]]];
        assert(pt->v[ipt] == numq);

        tmp = list[0];
        for(l=0;l<ilist-1;l++){
            list[l] = list[l+1];
        }
        list[ilist-1] = tmp;
    }

    return 2;
}

/**
 * \param mesh mesh
 * \param ip0 first edge extremity
 * \param ip1 second edge extremity
 * \param s parameter value
 * \param o point coordinates
 * \param no1 normal at point \a o (to fill)
 * \param no2 normal at point \a o (to fill)
 * \param to tangent at point \a o (to fill)
 *
 * \return 0 if fail,  1 otherwise
 *
 * \warning return 1 without filling \a no1 and \a no2 if the edge has 2
 * singular extremities.
 *
 * Compute point located at parameter value step from point ip0, as
 * well as interpolate of normals, tangent for a RIDGE edge
 *
 */
inline int MMG5_BezierRidge ( MMG5_pMesh mesh,MMG5_int ip0,MMG5_int ip1,double s,double *o,
                              double *no1,double *no2,double *to ) {
    MMG5_pPoint    p0,p1;
    double         ux,uy,uz,n01[3],n02[3],n11[3],n12[3],t0[3],t1[3];
    double         ps,ps2,b0[3],b1[3],bn[3],ll,il,ntemp[3],dd,alpha;

    p0 = &mesh->point[ip0];  /* Ref point, from which step is counted */
    p1 = &mesh->point[ip1];
    if ( !(MG_GEO & p0->tag) || !(MG_GEO & p1->tag) ) {
      return 0;
    }

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];
    ll = ux*ux + uy*uy + uz*uz;
    if ( ll < MMG5_EPSD2 )  {
      return 0;
    }
    il = 1.0 / sqrt(ll);

    if ( MG_SIN(p0->tag) ) {
        t0[0] = ux * il;
        t0[1] = uy * il;
        t0[2] = uz * il;
    }
    else {
      /* It is ok to pass here for nm points because they have always a tangent */
      memcpy(t0,&(p0->n[0]),3*sizeof(double));
      ps = t0[0]*ux + t0[1]*uy + t0[2]*uz;
      if ( ps < 0.0 ) {
        t0[0] *= -1.0;
        t0[1] *= -1.0;
        t0[2] *= -1.0;
      }
    }

    if ( MG_SIN(p1->tag) ) {
        t1[0] = -ux * il;
        t1[1] = -uy * il;
        t1[2] = -uz * il;
    }
    else {
        /* It is ok to pass here for nm points because they have always a tangent */
        memcpy(t1,&(p1->n[0]),3*sizeof(double));
        ps = - ( t1[0]*ux + t1[1]*uy + t1[2]*uz );
        if ( ps < 0.0 ) {
            t1[0] *= -1.0;
            t1[1] *= -1.0;
            t1[2] *= -1.0;
        }
    }
    alpha = MMG5_BezierGeod(p0->c,p1->c,t0,t1);

    b0[0] = p0->c[0] + alpha * t0[0];
    b0[1] = p0->c[1] + alpha * t0[1];
    b0[2] = p0->c[2] + alpha * t0[2];

    b1[0] = p1->c[0] + alpha * t1[0];
    b1[1] = p1->c[1] + alpha * t1[1];
    b1[2] = p1->c[2] + alpha * t1[2];

    o[0] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[0] + 3.0*s*(1.0-s)*(1.0-s)*b0[0] + \
        3.0*s*s*(1.0-s)*b1[0] + s*s*s*p1->c[0];

    o[1] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[1] + 3.0*s*(1.0-s)*(1.0-s)*b0[1] + \
        3.0*s*s*(1.0-s)*b1[1] + s*s*s*p1->c[1];

    o[2] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[2] + 3.0*s*(1.0-s)*(1.0-s)*b0[2] + \
        3.0*s*s*(1.0-s)*b1[2] + s*s*s*p1->c[2];

   if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
      /* In this case, the tangent orientation depends on the triangle from
       * which we see the ridge */
      memcpy(to,t0,3*sizeof(double));
      return 1;
    }
    else if ( MG_SIN(p0->tag) ) {
        memcpy(n11,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
        memcpy(n12,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));
        memcpy(n01,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
        memcpy(n02,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));
    }
    else if ( MG_SIN(p1->tag) ) {
        memcpy(n01,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
        memcpy(n02,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
        memcpy(n11,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
        memcpy(n12,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
    }
    else {
        memcpy(n01,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
        memcpy(n02,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
        memcpy(n11,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
        memcpy(n12,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));

        /* Switch normals of p1 for pairing */
        ps  = n01[0] * n11[0] + n01[1] * n11[1] + n01[2] * n11[2];
        ps2 = n01[0] * n12[0] + n01[1] * n12[1] + n01[2] * n12[2];
        if ( ps2 > ps ) {
            memcpy(ntemp,n11,3*sizeof(double));
            memcpy(n11,n12,3*sizeof(double));
            memcpy(n12,ntemp,3*sizeof(double));
        }
    }

    /* Normal n1 interpolation */
    ps = ux*(n01[0] + n11[0]) + uy*(n01[1] + n11[1]) + uz*(n01[2] + n11[2]);
    ps = 2.0*ps / ll;

    bn[0] = n01[0] + n11[0] -ps*ux;
    bn[1] = n01[1] + n11[1] -ps*uy;
    bn[2] = n01[2] + n11[2] -ps*uz;

    dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
    if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        bn[0] *= dd;
        bn[1] *= dd;
        bn[2] *= dd;
    }
    no1[0] = (1.0-s)*(1.0-s)*n01[0] + 2.0*s*(1.0-s)*bn[0] + s*s*n11[0];
    no1[1] = (1.0-s)*(1.0-s)*n01[1] + 2.0*s*(1.0-s)*bn[1] + s*s*n11[1];
    no1[2] = (1.0-s)*(1.0-s)*n01[2] + 2.0*s*(1.0-s)*bn[2] + s*s*n11[2];
    dd = no1[0]*no1[0] + no1[1]*no1[1] + no1[2]*no1[2];
    if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        no1[0] *= dd;
        no1[1] *= dd;
        no1[2] *= dd;
    }

    /* Normal n2 interpolation */
    ps = ux*(n02[0] + n12[0]) + uy*(n02[1] + n12[1]) + uz*(n02[2] + n12[2]);
    ps = 2.0*ps/ll;

    bn[0] = n02[0] + n12[0] -ps*ux;
    bn[1] = n02[1] + n12[1] -ps*uy;
    bn[2] = n02[2] + n12[2] -ps*uz;

    dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
    if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        bn[0] *= dd;
        bn[1] *= dd;
        bn[2] *= dd;
    }
    no2[0] = (1.0-s)*(1.0-s)*n02[0] + 2.0*s*(1.0-s)*bn[0] + s*s*n12[0];
    no2[1] = (1.0-s)*(1.0-s)*n02[1] + 2.0*s*(1.0-s)*bn[1] + s*s*n12[1];
    no2[2] = (1.0-s)*(1.0-s)*n02[2] + 2.0*s*(1.0-s)*bn[2] + s*s*n12[2];
    dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
    if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        no2[0] *= dd;
        no2[1] *= dd;
        no2[2] *= dd;

        to[0] = no1[1]*no2[2] - no1[2]*no2[1];
        to[1] = no1[2]*no2[0] - no1[0]*no2[2];
        to[2] = no1[0]*no2[1] - no1[1]*no2[0];
    }
    else {
      /* Open boundary: tangent interpolation (possibly flip (back) t1) */
      ps = t0[0]*t1[0] + t0[1]*t1[1] + t0[2]*t1[2];
      if ( ps < 0.0 ) {
        t1[0] *= -1.0;
        t1[1] *= -1.0;
        t1[2] *= -1.0;
      }
      to[0] = (1.0-s)*t0[0] + s*t1[0];
      to[1] = (1.0-s)*t0[1] + s*t1[1];
      to[2] = (1.0-s)*t0[2] + s*t1[2];

      /* Projection of the tangent in the tangent plane defined by no */
      ps = to[0]*no1[0] + to[1]*no1[1] + to[2]*no1[2];
      to[0] = to[0] -ps*no1[0];
      to[1] = to[1] -ps*no1[1];
      to[2] = to[2] -ps*no1[2];
    }

    dd = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
    if ( dd > MMG5_EPSD2 ) {
      dd = 1.0/sqrt(dd);
      to[0] *= dd;
      to[1] *= dd;
      to[2] *= dd;
    }

    return 1;
}

/**
 * \param mesh mesh
 * \param ip0 first edge extremity
 * \param ip1 second edge extremity
 * \param s parameter value
 * \param o point coordinates
 * \param no normal at point \a o (to fill)
 * \param to tangent at point \a o (to fill)
 *
 * \return 0 if fail,  1 otherwise
 *
 * \warning return 1 without filling \a no if the edge has 2 singular extremities.
 *
 * Compute point located at parameter value step from point ip0, as
 * well as interpolate of normals, tangent for a REF edge.
 *
 */
inline int MMG5_BezierRef ( MMG5_pMesh mesh,MMG5_int ip0,MMG5_int ip1,double s,double *o,
                            double *no,double *to) {
    MMG5_pPoint     p0,p1;
    double          ux,uy,uz,n01[3],n02[3],n11[3],n12[3],ntemp[3],t0[3],t1[3];
    double          ps,b0[3],b1[3],bn[3],ll,il,dd,alpha,ps2;

    p0 = &mesh->point[ip0];  /* Ref point, from which step is counted */
    p1 = &mesh->point[ip1];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];
    ll = ux*ux + uy*uy + uz*uz;
    if ( ll < MMG5_EPSD2 ) {
      return 0;
    }
    il = 1.0 / sqrt(ll);
    assert( (MG_REF & p0->tag) && (MG_REF & p1->tag) );

    /* Coordinates of the new point */
    if ( MG_SIN(p0->tag) ) {
        t0[0] = ux * il;
        t0[1] = uy * il;
        t0[2] = uz * il;
    }
    else {
        /* nm points have a tangent so its ok to get it here */
        memcpy(t0,&(p0->n[0]),3*sizeof(double));
        ps = t0[0]*ux + t0[1]*uy + t0[2]*uz;
        if ( ps < 0.0 ) {
            t0[0] *= -1.0;
            t0[1] *= -1.0;
            t0[2] *= -1.0;
        }
    }
    if ( MG_SIN(p1->tag) ) {
        t1[0] = -ux * il;
        t1[1] = -uy * il;
        t1[2] = -uz * il;
    }
    else {
        /* nm points have a tangent so its ok to get it here */
        memcpy(t1,&(p1->n[0]),3*sizeof(double));
        ps = - ( t1[0]*ux + t1[1]*uy + t1[2]*uz );
        if ( ps < 0.0 ) {
            t1[0] *= -1.0;
            t1[1] *= -1.0;
            t1[2] *= -1.0;
        }
    }

    alpha = MMG5_BezierGeod(p0->c,p1->c,t0,t1);

    b0[0] = p0->c[0] + alpha * t0[0];
    b0[1] = p0->c[1] + alpha * t0[1];
    b0[2] = p0->c[2] + alpha * t0[2];

    b1[0] = p1->c[0] + alpha * t1[0];
    b1[1] = p1->c[1] + alpha * t1[1];
    b1[2] = p1->c[2] + alpha * t1[2];

    o[0] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[0] + 3.0*s*(1.0-s)*(1.0-s)*b0[0] + \
        3.0*s*s*(1.0-s)*b1[0] + s*s*s*p1->c[0];

    o[1] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[1] + 3.0*s*(1.0-s)*(1.0-s)*b0[1] + \
        3.0*s*s*(1.0-s)*b1[1] + s*s*s*p1->c[1];

    o[2] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[2] + 3.0*s*(1.0-s)*(1.0-s)*b0[2] + \
        3.0*s*s*(1.0-s)*b1[2] + s*s*s*p1->c[2];

    /* Coordinates of the new tangent and normal */
    if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
        memcpy(to,t0,3*sizeof(double));
        return 1;
    }
    else if ( MG_SIN(p0->tag) ) {
        memcpy(n11,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
        memcpy(n01,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
        memcpy(n12,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));
        memcpy(n02,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));
    }
    else if ( MG_SIN(p1->tag) ) {
        memcpy(n01,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
        memcpy(n11,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
        memcpy(n02,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
        memcpy(n12,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
    }
    else {
        memcpy(n01,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
        memcpy(n11,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
        memcpy(n02,&(mesh->xpoint[p0->xp].n2[0]),3*sizeof(double));
        memcpy(n12,&(mesh->xpoint[p1->xp].n2[0]),3*sizeof(double));
    }

    /* Switch normals of p1 for pairing */
    ps  = n01[0] * n11[0] + n01[1] * n11[1] + n01[2] * n11[2];
    ps2 = n01[0] * n12[0] + n01[1] * n12[1] + n01[2] * n12[2];
    if ( ps2 > ps ) {
        memcpy(ntemp,n11,3*sizeof(double));
        memcpy(n11,n12,3*sizeof(double));
        memcpy(n12,ntemp,3*sizeof(double));
    }

    /* Normal interpolation */
    ps = ux*(n01[0] + n11[0]) + uy*(n01[1] + n11[1]) + uz*(n01[2] + n11[2]);
    ps = 2.0*ps/ll;

    bn[0] = n01[0] + n11[0] -ps*ux;
    bn[1] = n01[1] + n11[1] -ps*uy;
    bn[2] = n01[2] + n11[2] -ps*uz;

    dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
    if ( dd > MMG5_EPSD ) {
        dd = 1.0 / sqrt(dd);
        bn[0] *= dd;
        bn[1] *= dd;
        bn[2] *= dd;
    }
    no[0] = (1.0-s)*(1.0-s)*n01[0] + 2.0*s*(1.0-s)*bn[0] + s*s*n11[0];
    no[1] = (1.0-s)*(1.0-s)*n01[1] + 2.0*s*(1.0-s)*bn[1] + s*s*n11[1];
    no[2] = (1.0-s)*(1.0-s)*n01[2] + 2.0*s*(1.0-s)*bn[2] + s*s*n11[2];

    dd = no[0]*no[0] + no[1]*no[1] + no[2]*no[2];
    if ( dd > MMG5_EPSD2 ) {
        dd = 1.0/sqrt(dd);
        no[0] *= dd;
        no[1] *= dd;
        no[2] *= dd;
    }

    /* Tangent interpolation : possibly flip (back) t1 */
    ps = t0[0]*t1[0] + t0[1]*t1[1] + t0[2]*t1[2];
    if ( ps < 0.0 ) {
        t1[0] *= -1.0;
        t1[1] *= -1.0;
        t1[2] *= -1.0;
    }
    to[0] = (1.0-s)*t0[0] + s*t1[0];
    to[1] = (1.0-s)*t0[1] + s*t1[1];
    to[2] = (1.0-s)*t0[2] + s*t1[2];

    /* Projection of the tangent in the tangent plane defined by no */
    ps = to[0]*no[0] + to[1]*no[1] + to[2]*no[2];
    to[0] = to[0] -ps*no[0];
    to[1] = to[1] -ps*no[1];
    to[2] = to[2] -ps*no[2];

    dd = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
    if ( dd > MMG5_EPSD2) {
        dd = 1.0 / sqrt(dd);
        to[0] *= dd;
        to[1] *= dd;
        to[2] *= dd;
    }

    return 1;
}

/**
 * \param mesh mesh
 * \param ip0 first edge extremity
 * \param ip1 second edge extremity
 * \param s parameter value
 * \param o point coordinates
 * \param no normal at point \a o (to fill)
 * \param to tangent at point \a o along edge ip0 ip1 (to fill)
 *
 * \return 0 if fail,  1 otherwise
 *
 * \warning return 1 without filling \a no if the edge has 2 singular extremities.
 *
 * Compute point located at parameter value \a s from point ip0, as
 * well as interpolate of normals, tangent for a NOM edge
 *
 */
inline int MMG5_BezierNom(MMG5_pMesh mesh,MMG5_int ip0,MMG5_int ip1,double s,double *o,double *no,double *to) {
    MMG5_pPoint p0,p1;
    double      ux,uy,uz,il,ll,ps,alpha,dd;
    double      t0[3],t1[3],b0[3],b1[3],n0[3],n1[3],bn[3];
    int8_t      intnom;

    p0 = &mesh->point[ip0];
    p1 = &mesh->point[ip1];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];
    ll = ux*ux + uy*uy + uz*uz;

    if(ll < MMG5_EPSD2) return 0;
    il = 1.0 / sqrt(ll);

    assert(( p0->tag & MG_NOM ) && ( p1->tag & MG_NOM ));

    /* Coordinates of the new point */
    if ( MG_SIN(p0->tag) ) {
        t0[0] = ux * il;
        t0[1] = uy * il;
        t0[2] = uz * il;
    }
    else {
        memcpy(t0,&(p0->n[0]),3*sizeof(double));
        ps = t0[0]*ux + t0[1]*uy + t0[2]*uz;
        if ( ps < 0.0 ) {
            t0[0] *= -1.0;
            t0[1] *= -1.0;
            t0[2] *= -1.0;
        }
    }
    if ( MG_SIN(p1->tag) ) {
        t1[0] = -ux * il;
        t1[1] = -uy * il;
        t1[2] = -uz * il;
    }
    else {
        memcpy(t1,&(p1->n[0]),3*sizeof(double));
        ps = - ( t1[0]*ux + t1[1]*uy + t1[2]*uz );
        if ( ps < 0.0 ) {
            t1[0] *= -1.0;
            t1[1] *= -1.0;
            t1[2] *= -1.0;
        }
    }

    alpha = MMG5_BezierGeod(p0->c,p1->c,t0,t1);

    b0[0] = p0->c[0] + alpha * t0[0];
    b0[1] = p0->c[1] + alpha * t0[1];
    b0[2] = p0->c[2] + alpha * t0[2];

    b1[0] = p1->c[0] + alpha * t1[0];
    b1[1] = p1->c[1] + alpha * t1[1];
    b1[2] = p1->c[2] + alpha * t1[2];

    o[0] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[0] + 3.0*s*(1.0-s)*(1.0-s)*b0[0] + \
        3.0*s*s*(1.0-s)*b1[0] + s*s*s*p1->c[0];

    o[1] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[1] + 3.0*s*(1.0-s)*(1.0-s)*b0[1] + \
        3.0*s*s*(1.0-s)*b1[1] + s*s*s*p1->c[1];

    o[2] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[2] + 3.0*s*(1.0-s)*(1.0-s)*b0[2] + \
        3.0*s*s*(1.0-s)*b1[2] + s*s*s*p1->c[2];

    /* Coordinates of the new tangent and normal */
    if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {  // function should not be used in that case
        memcpy(to,t0,3*sizeof(double));
        return 1;
    }
    else if ( MG_SIN(p0->tag) ) {
        memcpy(n1,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
        memcpy(n0,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
    }
    else if ( MG_SIN(p1->tag) ) {
        memcpy(n0,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
        memcpy(n1,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
    }
    else {
        memcpy(n0,&(mesh->xpoint[p0->xp].n1[0]),3*sizeof(double));
        memcpy(n1,&(mesh->xpoint[p1->xp].n1[0]),3*sizeof(double));
    }

    /* Check for internal non manifold edge */
    intnom = ( mesh->xpoint[p0->xp].nnor ) || ( mesh->xpoint[p1->xp].nnor );

    /* Normal interpolation */
    if ( !intnom ) {
      ps = ux*(n0[0] + n1[0]) + uy*(n0[1] + n1[1]) + uz*(n0[2] + n1[2]);
      ps = 2.0*ps/ll;

      bn[0] = n0[0] + n1[0] -ps*ux;
      bn[1] = n0[1] + n1[1] -ps*uy;
      bn[2] = n0[2] + n1[2] -ps*uz;

      dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
      if ( dd > MMG5_EPSD ) {
          dd = 1.0 / sqrt(dd);
          bn[0] *= dd;
          bn[1] *= dd;
          bn[2] *= dd;
      }
      no[0] = (1.0-s)*(1.0-s)*n0[0] + 2.0*s*(1.0-s)*bn[0] + s*s*n1[0];
      no[1] = (1.0-s)*(1.0-s)*n0[1] + 2.0*s*(1.0-s)*bn[1] + s*s*n1[1];
      no[2] = (1.0-s)*(1.0-s)*n0[2] + 2.0*s*(1.0-s)*bn[2] + s*s*n1[2];

      dd = no[0]*no[0] + no[1]*no[1] + no[2]*no[2];
      if ( dd > MMG5_EPSD2 ) {
          dd = 1.0/sqrt(dd);
          no[0] *= dd;
          no[1] *= dd;
          no[2] *= dd;
      }
    }

    /* Tangent interpolation : possibly flip (back) t1 */
    ps = t0[0]*t1[0] + t0[1]*t1[1] + t0[2]*t1[2];
    if ( ps < 0.0 ) {
        t1[0] *= -1.0;
        t1[1] *= -1.0;
        t1[2] *= -1.0;
    }
    to[0] = (1.0-s)*t0[0] + s*t1[0];
    to[1] = (1.0-s)*t0[1] + s*t1[1];
    to[2] = (1.0-s)*t0[2] + s*t1[2];

    /* Projection of the tangent in the tangent plane defined by no */
    if ( !intnom ) {
      ps = to[0]*no[0] + to[1]*no[1] + to[2]*no[2];
      to[0] = to[0] -ps*no[0];
      to[1] = to[1] -ps*no[1];
      to[2] = to[2] -ps*no[2];
    }

    dd = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];
    if ( dd > MMG5_EPSD2) {
        dd = 1.0 / sqrt(dd);
        to[0] *= dd;
        to[1] *= dd;
        to[2] *= dd;
    }

    return 1;
}


/**
 * \param mesh mesh
 * \param ip0 first edge extremity
 * \param ip1 second edge extremity
 * \param s parameter value
 * \param v reference normal
 * \param o point coordinates
 * \param no normal at point \a o (to fill)
 *
 * \return 0 if fail,  1 otherwise
 *
 * \warning return 1 without filling \a no if the edge has 2 singular extremities.
 *
 * Compute point located at parameter value step from point ip0, as
 * well as interpolate of normals, tangent for a regular edge ; v =
 * ref vector (normal) for choice of normals if need be
 *
 */
inline int MMG5_BezierReg(MMG5_pMesh mesh,MMG5_int ip0, MMG5_int ip1, double s, double v[3], double *o, double *no){
    MMG5_pPoint p0,p1;
    double      b0[3],b1[3],bn[3],t0[3],t1[3],np0[3],np1[3],alpha,ux,uy,uz,ps1,ps2,ll;
    double      il,dd,*n1,*n2;

    p0 = &mesh->point[ip0];
    p1 = &mesh->point[ip1];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    uz = p1->c[2] - p0->c[2];

    ll = ux*ux + uy*uy + uz*uz;

    np0[0] = np0[1] = np0[2] = 0;
    np1[0] = np1[1] = np1[2] = 0;

    /* Pathological case : don't move in that case ! */
    if ( (MG_SIN_OR_NOM(p0->tag) && MG_SIN_OR_NOM(p1->tag)) || (ll<MMG5_EPSD) ){
        o[0] = 0.5*( p0->c[0] + p1->c[0] );
        o[1] = 0.5*( p0->c[1] + p1->c[1] );
        o[2] = 0.5*( p0->c[2] + p1->c[2] );

        memcpy(no,v,3*sizeof(double));
        return 1;
    }

    il = 1.0 /sqrt(ll);

    /* Coordinates of the new tangent and normal */
    if ( MG_SIN_OR_NOM(p0->tag) ) {
        if(p1->tag & MG_GEO){
            n1 = &mesh->xpoint[p1->xp].n1[0];
            n2 = &mesh->xpoint[p1->xp].n2[0];
            ps1 = n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2];
            ps2 = n2[0]*v[0] + n2[1]*v[1] + n2[2]*v[2];

            if(fabs(ps1) < fabs(ps2)){
                memcpy(np1,&mesh->xpoint[p1->xp].n2[0],3*sizeof(double));
            }
            else{
                memcpy(np1,&mesh->xpoint[p1->xp].n1[0],3*sizeof(double));
            }
        }
        else{
            memcpy(np1,&mesh->xpoint[p1->xp].n1[0],3*sizeof(double));
        }
        memcpy(np0,np1,3*sizeof(double));
    }
    else if ( MG_SIN_OR_NOM(p1->tag) ) {
        if(p0->tag & MG_GEO){
            n1 = &mesh->xpoint[p0->xp].n1[0];
            n2 = &mesh->xpoint[p0->xp].n2[0];
            ps1 = n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2];
            ps2 = n2[0]*v[0] + n2[1]*v[1] + n2[2]*v[2];

            if(fabs(ps1) < fabs(ps2)){
                memcpy(np0,&mesh->xpoint[p0->xp].n2[0],3*sizeof(double));
            }
            else{
                memcpy(np0,&mesh->xpoint[p0->xp].n1[0],3*sizeof(double));
            }
        }
        else{
            memcpy(np0,&mesh->xpoint[p0->xp].n1[0],3*sizeof(double));
        }
        memcpy(np1,np0,3*sizeof(double));
    }
    else{
        if(p0->tag & MG_GEO){
            n1 = &mesh->xpoint[p0->xp].n1[0];
            n2 = &mesh->xpoint[p0->xp].n2[0];
            ps1 = n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2];
            ps2 = n2[0]*v[0] + n2[1]*v[1] + n2[2]*v[2];

            if(fabs(ps1) < fabs(ps2)){
                memcpy(np0,&mesh->xpoint[p0->xp].n2[0],3*sizeof(double));
            }
            else{
                memcpy(np0,&mesh->xpoint[p0->xp].n1[0],3*sizeof(double));
            }
        }
        else{
            memcpy(np0,&mesh->xpoint[p0->xp].n1[0],3*sizeof(double));
        }

        if(p1->tag & MG_GEO){
            n1 = &mesh->xpoint[p1->xp].n1[0];
            n2 = &mesh->xpoint[p1->xp].n2[0];
            ps1 = n1[0]*v[0] + n1[1]*v[1] + n1[2]*v[2];
            ps2 = n2[0]*v[0] + n2[1]*v[1] + n2[2]*v[2];

            if(fabs(ps1) < fabs(ps2)){
                memcpy(np1,&mesh->xpoint[p1->xp].n2[0],3*sizeof(double));
            }
            else{
                memcpy(np1,&mesh->xpoint[p1->xp].n1[0],3*sizeof(double));
            }
        }
        else{
            memcpy(np1,&mesh->xpoint[p1->xp].n1[0],3*sizeof(double));
        }
    }

    /* Normal interpolation */
    ps1 = ux*(np0[0] + np1[0]) + uy*(np0[1] + np1[1]) + uz*(np0[2] + np1[2]);
    ps1 = 2.0*ps1/ll;

    bn[0] = np0[0] + np1[0] -ps1*ux;
    bn[1] = np0[1] + np1[1] -ps1*uy;
    bn[2] = np0[2] + np1[2] -ps1*uz;

    dd = bn[0]*bn[0] + bn[1]*bn[1] + bn[2]*bn[2];
    if(dd > MMG5_EPSD){
        dd = 1.0/sqrt(dd);
        bn[0] *= dd;
        bn[1] *= dd;
        bn[2] *= dd;
    }

    no[0] = (1.0-s)*(1.0-s)*np0[0] + 2.0*s*(1.0-s)*bn[0] + s*s*np1[0];
    no[1] = (1.0-s)*(1.0-s)*np0[1] + 2.0*s*(1.0-s)*bn[1] + s*s*np1[1];
    no[2] = (1.0-s)*(1.0-s)*np0[2] + 2.0*s*(1.0-s)*bn[2] + s*s*np1[2];

    dd = no[0]*no[0] + no[1]*no[1] + no[2]*no[2];
    if(dd > MMG5_EPSD){
        dd = 1.0/sqrt(dd);
        no[0] *= dd;
        no[1] *= dd;
        no[2] *= dd;
    }

    /* vertex position interpolation */
    if(!MMG5_BezierTgt(p0->c,p1->c,np0,np1,t0,t1)){
        t0[0] = ux * il;
        t0[1] = uy * il;
        t0[2] = uz * il;

        t1[0] = - ux * il;
        t1[1] = - uy * il;
        t1[2] = - uz * il;
    }

    alpha = MMG5_BezierGeod(p0->c,p1->c,t0,t1);

    b0[0] = p0->c[0] + alpha * t0[0];
    b0[1] = p0->c[1] + alpha * t0[1];
    b0[2] = p0->c[2] + alpha * t0[2];

    b1[0] = p1->c[0] + alpha * t1[0];
    b1[1] = p1->c[1] + alpha * t1[1];
    b1[2] = p1->c[2] + alpha * t1[2];

    o[0] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[0] + 3.0*s*(1.0-s)*(1.0-s)*b0[0] + \
        3.0*s*s*(1.0-s)*b1[0] + s*s*s*p1->c[0];

    o[1] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[1] + 3.0*s*(1.0-s)*(1.0-s)*b0[1] + \
        3.0*s*s*(1.0-s)*b1[1] + s*s*s*p1->c[1];

    o[2] = (1.0-s)*(1.0-s)*(1.0-s)*p0->c[2] + 3.0*s*(1.0-s)*(1.0-s)*b0[2] + \
        3.0*s*s*(1.0-s)*b1[2] + s*s*s*p1->c[2];

    return 1;
}

/** find the element number in packed numerotation */
MMG5_int MMG3D_indElt(MMG5_pMesh mesh, MMG5_int kel) {
    MMG5_pTetra pt;
    MMG5_int    ne, k;

#ifndef USE_SCOTCH
    ne = mesh->ne+1;
    k  = 1;

    assert ( MG_EOK(&mesh->tetra[kel]) );

    do {
      pt = &mesh->tetra[k];

      if ( !MG_EOK(pt) ) {
        --ne;
        pt = &mesh->tetra[ne];
        assert( pt );

        /* Search last used tetra */
        while ( !MG_EOK(pt) && k < ne ) {
          --ne;
          pt = &mesh->tetra[ne];
        }

        /* Found elt */
        if ( ne == kel ) {
          return k;
        }
      }
      else {
        /* Found elt */
        if ( k==kel ) {
          return k;
        }
      }
      /* All elts have been treated end of loop */
      if ( k==ne ) {
        break;
      }
    }
    while ( ++k < ne );

#else
    ne = 0;
    for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( MG_EOK(pt) ) {
            ne++;
            if ( k == kel )  return ne;
        }
    }
#endif
    return 0;
}

/** find the point number in packed numerotation */
MMG5_int MMG3D_indPt(MMG5_pMesh mesh, MMG5_int kp) {
    MMG5_pPoint ppt;
    MMG5_int    np, k;

#ifndef USE_SCOTCH
    np = mesh->np+1;
    k  = 1;

    assert ( MG_VOK(&mesh->point[kp]) );

    do {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) ) {
        --np;
        ppt = &mesh->point[np];
        assert ( ppt );

        /* Search the last used point */
        while ( !MG_VOK(ppt) && k < np ) {
          --np;
          ppt = &mesh->point[np];
        }

        /* Found pt */
        if ( np == kp ) {
          return k;
        }
      }
      else {
        /* Found pt */
        if ( k==kp ) {
          return k;
        }
      }

      /* All points have been treated end of loop */
      if ( k==np ) {
        break;
      }
    }
    while ( ++k < np );

#else
    np = 0;
    for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( MG_VOK(ppt) ) {
            np++;
            if ( k == kp )  return np;
        }
    }
#endif

    return 0;
}

/** Debug function (not use in clean code): print mesh->tetra structure */
void MMG5_printTetra(MMG5_pMesh mesh,char* fileName) {
    MMG5_pTetra  pt;
    MMG5_pxTetra pxt;
    MMG5_int     k;
    FILE         *inm;

    inm = fopen(fileName,"w");

    fprintf(inm,"----------> %" MMG5_PRId " MMG5_TETRAHEDRAS <----------\n",mesh->ne);
    for(k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        fprintf(inm,"num %" MMG5_PRId " -> %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",k,pt->v[0],pt->v[1],
                pt->v[2],pt->v[3]);
        fprintf(inm,"ref,tag,xt  -> %" MMG5_PRId " %d %" MMG5_PRId "\n",pt->ref,pt->tag,pt->xt);
        if ( pt->xt ) {
            pxt = &mesh->xtetra[pt->xt];
            fprintf(inm,"tag   -> %d %d %d %d %d %d\n",pxt->tag[0],pxt->tag[1],
                    pxt->tag[2],pxt->tag[3],pxt->tag[4],pxt->tag[5]);
            fprintf(inm,"edg   -> %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",pxt->edg[0],pxt->edg[1],
                    pxt->edg[2],pxt->edg[3],pxt->edg[4],pxt->edg[5]);
            fprintf(inm,"ftag  -> %d %d %d %d\n",pxt->ftag[0],pxt->ftag[1],
                    pxt->ftag[2],pxt->ftag[3]);
            fprintf(inm,"ref   -> %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",pxt->ref[0],pxt->ref[1],
                    pxt->ref[2],pxt->ref[3]);
            fprintf(inm,"ori   -> %d \n",pxt->ori);
        }
        fprintf(inm,"\n");
    }
    fprintf(inm,"---------> END MMG5_TETRAHEDRAS <--------\n");
    fclose(inm);
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param ip global index of point in which we want to compute the local parameters
 * \param listv pointer toward the ball of \a ip
 * \param ilistv number of tetra in the ball of \a ip
 * \param lists pointer toward the surface ball of \a ip
 * \param ilists number of tetra in the surface ball of \a ip
 * \param hausd_ip pointer toward the local hausdorff parameter to compute
 * \param hmin_ip pointer toward the local minimal edge size to compute
 * \param hmax_ip pointer toward the local maximal edge size to compute
 *
 * \return 1 if success, 0 if fail
 *
 * Compute the local parameters at point \a ip (the volume and surface
 * ball of point must be provided).
 *
 */
int MMG3D_localParamReg(MMG5_pMesh mesh,MMG5_int ip,int64_t *listv,int ilistv,
                         MMG5_int *lists,int ilists,
                         double* hausd_ip,double *hmin_ip,double *hmax_ip) {

  MMG5_pTetra pt;
  MMG5_pPar   par;
  double      hausd, hmin, hmax;
  int         l,k,isloc,ifac1;

  hausd = mesh->info.hausd;
  hmin  = mesh->info.hmin;
  hmax  = mesh->info.hmax;
  isloc = 0;

  /* travel across the ball of ip to find the minimal local params imposed on
   * tetras */
  if ( mesh->info.parTyp & MG_Tetra ) {
    l = 0;
    do
    {
      if ( isloc )  break;

      par = &mesh->info.par[l];
      if ( par->elt != MMG5_Tetrahedron )  continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/4];
        if ( par->ref == pt->ref ) {
          hausd = par->hausd;
          hmin  = par->hmin;
          hmax  = par->hmax;
          isloc = 1;
          break;
        }
      }
    } while ( ++l<mesh->info.npar );

    for ( ; l<mesh->info.npar; ++l ) {
      par = &mesh->info.par[l];
      if ( par->elt != MMG5_Tetrahedron ) continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/4];
        if ( par->ref == pt->ref ) {
          hausd = MG_MIN(hausd,par->hausd);
          hmin  = MG_MAX(hmin,par->hmin);
          hmax  = MG_MIN(hmax,par->hmax);
          break;
        }
      }
    }
  }
  /* travel across the surface ball of ip to find the minimal local params
   * imposed on trias */
  if ( mesh->info.parTyp & MG_Tria ) {
    l = 0;
    do
    {
      if ( isloc )  break;

      par = &mesh->info.par[l];
      if ( par->elt != MMG5_Triangle )  continue;

      for ( k=0; k<ilists; ++k ) {
        pt = &mesh->tetra[lists[k]/4];
        ifac1 =  lists[k] % 4;
        assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac1] & MG_BDY) );

        if ( par->ref == mesh->xtetra[pt->xt].ref[ifac1] ) {
          hausd = par->hausd;
          hmin  = par->hmin;
          hmax  = par->hmax;
          isloc = 1;
          break;
        }
      }
    } while ( ++l<mesh->info.npar );

    for ( ; l<mesh->info.npar; ++l ) {
      par = &mesh->info.par[l];
      if ( par->elt != MMG5_Triangle ) continue;

      for ( k=0; k<ilists; ++k ) {
        pt = &mesh->tetra[lists[k]/4];
        ifac1 =  lists[k] % 4;
        assert(pt->xt && (mesh->xtetra[pt->xt].ftag[ifac1] & MG_BDY) );

        if ( par->ref == mesh->xtetra[pt->xt].ref[ifac1] ) {
          hausd = MG_MIN(hausd,par->hausd);
          hmin  = MG_MAX(hmin,par->hmin);
          hmax  = MG_MIN(hmax,par->hmax);
          break;
        }
      }
    }
  }

  /* Return the wanted values */
  if ( hausd_ip ) *hausd_ip = hausd;
  if ( hmin_ip ) *hmin_ip = hmin;
  if ( hmax_ip ) *hmax_ip = hmax;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param iel index of tetra in which we work
 * \param iface index of face in \a iel
 * \param ia index of edge in \a iel along which we want to compute the local
 * parameters
 * \param hausd_ip pointer toward the local hausdorff parameter to compute
 * \param hmin_ip pointer toward the local minimal edge size to compute
 * \param hmax_ip pointer toward the local maximal edge size to compute
 *
 * \return 1 if success, 0 if fail
 *
 * Compute the local parameters at non manifold point \a ip.
 *
 */
int MMG3D_localParamNm(MMG5_pMesh mesh,MMG5_int iel,int iface,int ia,
                         double* hausd_ip,double *hmin_ip,double *hmax_ip) {

  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPar     par;
  double        hausd, hmin, hmax;
  int           l,k,isloc;
  int           ilistv;
  int64_t       listv[MMG3D_LMAX+2];
  MMG5_int      ifac1,ifac2;
  static int8_t mmgWarn0;


  hausd = mesh->info.hausd;
  hmin  = mesh->info.hmin;
  hmax  = mesh->info.hmax;
  isloc = 0;

  pt = &mesh->tetra[iel];
  pxt = &mesh->xtetra[pt->xt];

  /* Warning : rough eval of the local param at triangles if coquilface
   * fails because we have more than 2 boundaries in the edge shell
   * (non-manifold domain). In this case, we just take into account 2
   * boundaries of the shell */

  if ( pxt->tag[ia] & MG_OPNBDY ) {
    ilistv = 1;
    ifac1  = ifac2 = 4*iel + iface;
  }
  else {
    ilistv = MMG5_coquilface( mesh,iel,iface, ia,listv,&ifac1,&ifac2,1);
  }
  if ( ilistv < 0 )
  {
    if ( mesh->info.ddebug || mesh->info.imprim>5 ) {
      if ( !mmgWarn0 ) {
        mmgWarn0 = 1;
        fprintf(stderr, "  ## Warning: %s: unable to take into account local"
                " parameters at at least 1 vertex.\n",__func__ );
      }
    }

    if ( mesh->info.parTyp & MG_Tria ) {
      for ( l=0; l<mesh->info.npar; ++l) {
        par = &mesh->info.par[l];
        if ( par->elt != MMG5_Triangle ) continue;

        if ( pxt->ref[iface]!=par->ref ) continue;

        if ( isloc ) {
          hausd   = MG_MIN(hausd,par->hausd);
          hmin    = MG_MAX(hmin,par->hmin);
          hmax    = MG_MIN(hmax,par->hmax);
        }
        else {
          hausd   = par->hausd;
          hmin    = par->hmin;
          hmax    = par->hmax;
          isloc   = 1;
        }
      }
    }

  }
  else {

    /* Local params at triangles containing the edge (not optimal) */
    if ( mesh->info.parTyp & MG_Tria ) {
      for ( l=0; l<mesh->info.npar; ++l) {
        par = &mesh->info.par[l];
        if ( par->elt != MMG5_Triangle ) continue;

        if ( mesh->xtetra[mesh->tetra[ifac1/4].xt].ref[ifac1%4]!=par->ref &&
             mesh->xtetra[mesh->tetra[ifac2/4].xt].ref[ifac2%4]!=par->ref )
          continue;

        if ( isloc ) {
          hausd   = MG_MIN(hausd,par->hausd);
          hmin    = MG_MAX(hmin,par->hmin);
          hmax    = MG_MIN(hmax,par->hmax);
        }
        else {
          hausd   = par->hausd;
          hmin    = par->hmin;
          hmax    = par->hmax;
          isloc   = 1;
        }
      }
    }
  }

  /* Local params at tetra of the edge shell */
  if ( mesh->info.parTyp & MG_Tetra ) {
    ilistv/=2;
    l = 0;
    do
    {
      if ( isloc )  break;

      par = &mesh->info.par[l];
      if ( par->elt != MMG5_Tetrahedron ) continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/6];
        if ( par->ref != pt->ref ) continue;

        hausd   = par->hausd;
        hmin    = par->hmin;
        hmax    = par->hmax;
        isloc   = 1;
      }
    } while ( ++l<mesh->info.npar );

    for ( ; l<mesh->info.npar; ++l ) {
      par = &mesh->info.par[l];
      if ( par->elt != MMG5_Tetrahedron ) continue;

      for ( k=0; k<ilistv; ++k ) {
        pt = &mesh->tetra[listv[k]/6];
        if ( par->ref != pt->ref ) continue;

        hausd = MG_MIN(hausd,par->hausd);
        hmin  = MG_MAX(hmin,par->hmin);
        hmax  = MG_MIN(hmax,par->hmax);
        break;
      }
    }
  }
  /* Return the wanted values */
  if ( hausd_ip ) *hausd_ip = hausd;
  if ( hmin_ip ) *hmin_ip = hmin;
  if ( hmax_ip ) *hmax_ip = hmax;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * Mark the mesh vertices that belong to triangles or quadrangles as used (for
 * Mmgs or Mmg2d).
 *
 */
static inline
void MMG3D_mark_usedVertices ( MMG5_pMesh mesh ) {
  MMG5_pTetra pt;
  MMG5_pPrism pq;
  MMG5_pPoint ppt;
  int         i;
  MMG5_int    k;

  /* Preserve isolated required points */
  for ( k=1; k<=mesh->np; k++ ) {
    ppt = &mesh->point[k];
    if ( ppt->flag || !(ppt->tag & MG_REQ) ) {
      continue;
    }
    ppt->tag &= ~MG_NUL;
  }

  /* Mark points used by the connectivity */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<4; i++) {
      ppt = &mesh->point[ pt->v[i] ];
      ppt->tag &= ~MG_NUL;
    }
  }

  for ( k=1; k<=mesh->nprism; k++ ) {
    pq = &mesh->prism[k];
    if ( !MG_EOK(pq) )  continue;

    for (i=0; i<6; i++) {
      ppt = &mesh->point[ pq->v[i] ];
      ppt->tag &= ~MG_NUL;
    }
  }

  /* Finally, clean point array */
  while ( (!MG_VOK(&mesh->point[mesh->np])) && mesh->np ) {
    MMG3D_delPt(mesh,mesh->np);
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param nsd subdomain index.
 *
 * Remove tetra that do not belong to subdomain of index \a nsd
 *
 */
static
void MMG3D_keep_subdomainElts ( MMG5_pMesh mesh, int nsd ) {
  MMG5_pTetra pt;
  int         i,iv;
  MMG5_int    k,*adja,iadr,iadrv;
  int         nfac = 4; // number of faces per elt

  for ( k=1 ; k <= mesh->ne ; k++) {
    pt = &mesh->tetra[k];

    if ( !MG_EOK(pt) ) continue;

    /* Mark tetra vertices as seen to be able to detect isolated points */
    mesh->point[pt->v[0]].flag = 1;
    mesh->point[pt->v[1]].flag = 1;
    mesh->point[pt->v[2]].flag = 1;
    mesh->point[pt->v[3]].flag = 1;

    if ( pt->ref == nsd ) continue;

    /* Update adjacency relationship: we will delete elt k so k adjacent will
     * not be adjacent to k anymore */
    if ( mesh->adja ) {
      iadr = nfac*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for ( i=0; i<nfac; ++i ) {
        iadrv = adja[i];
        if ( !iadrv ) {
          continue;
        }
        iv = iadrv%nfac;
        iadrv /= nfac;
        mesh->adja[nfac*(iadrv-1)+1+iv] = 0;
      }
    }

    /* Delete element */
    MMG3D_delElt(mesh,k);
  }

  return;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param nsd index of subdomain to keep.
 *
 * Keep only subdomain of index \a nsd and remove other subdomains.
 *
 */
void MMG3D_keep_only1Subdomain ( MMG5_pMesh mesh,int nsd ) {

  if ( !nsd ) {
    return;
  }

  if ( mesh->info.imprim > 4 || mesh->info.ddebug ) {
    fprintf(stdout,"\n  -- ONLY KEEP DOMAIN OF REF %d\n",nsd );
  }

  MMG5_mark_verticesAsUnused ( mesh );

  MMG3D_keep_subdomainElts ( mesh, nsd );

  MMG3D_mark_usedVertices ( mesh );

  return;
}
