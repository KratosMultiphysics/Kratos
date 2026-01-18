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
 * \file mmg3d/analys_3d.c
 * \brief Mesh analysis.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"

/**
 * \param mesh pointer towarad the mesh structure.
 *
 * Set all boundary triangles to required and add a tag to detect that they are
 * not realy required.
 *
 */
void MMG3D_set_reqBoundaries(MMG5_pMesh mesh) {
  MMG5_pTria     ptt;
  MMG5_int       k;

  /* The MG_REQ+MG_NOSURF tag mark the boundary edges that we dont want to touch
   * but that are not really required (-nosurf option) */
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];

    if ( mesh->info.nosurf  && (!(ptt->tag[0] & MG_REQ)) ) {
      ptt->tag[0] |= MG_REQ;
      ptt->tag[0] |= MG_NOSURF;
    }

    if ( ptt->tag[0] & MG_PARBDY ) {
      ptt->tag[0] |= MG_NOSURF;
      ptt->tag[0] |= MG_REQ;
    }

    if ( mesh->info.nosurf && (!(ptt->tag[1] & MG_REQ)) ) {
      ptt->tag[1] |= MG_REQ;
      ptt->tag[1] |= MG_NOSURF;
    }

    if ( ptt->tag[1] & MG_PARBDY ) {
      ptt->tag[1] |= MG_NOSURF;
      ptt->tag[1] |= MG_REQ;
    }

    if ( mesh->info.nosurf && (!(ptt->tag[2] & MG_REQ)) ) {
      ptt->tag[2] |= MG_REQ;
      ptt->tag[2] |= MG_NOSURF;
    }

    if ( ptt->tag[2] & MG_PARBDY ) {
      ptt->tag[2] |= MG_NOSURF;
      ptt->tag[2] |= MG_REQ;
    }
  }

  return;
}


/**
 * \param mesh pointer towarad the mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * topology: set tria adjacency, detect Moebius, flip faces, count connected
 * comp.
 *
 * \remark: as all triangles are mesh boundaries, we do not need to mark their
 * adges as MG_BDY so the MG_BDY tag may be used inside geometrical triangles
 * (external non-parallel, or internal parallel) to tag edges on the
 * intersection with purely parallel (non-geometrical) triangles.
 * The MG_PARBDYBDY tag is also added, as it does not have a supporting triangle
 * to inherit this tag from.
 *
 * \remark REQ, NOSURF, etc... tags are added only inside xtetra.
 * \remark In openbdy mode, all non-manifold edges are marked as opnbdy.
 *
 */
int MMG5_setadj(MMG5_pMesh mesh){
  MMG5_pTria   pt,pt1;
  MMG5_int     *adja,*adjb,adji1,adji2,*pile,iad,ipil,ip1,ip2,gen;
  MMG5_int     k,kk,iel,jel,nvf,nf,nr,nm,nt,nre,nreq,ncc,ned,ref;
  uint16_t     tag;
  int8_t       i,ii,i1,i2,ii1,ii2,voy;

  nvf = nf = ncc = ned = 0;

  MMG5_SAFE_MALLOC(pile,mesh->nt+1,MMG5_int,return 0);

  pile[1] = 1;
  ipil    = 1;

  while ( ipil > 0 ) {
    ncc++;

    do {
      k  = pile[ipil--];
      pt = &mesh->tria[k];
      pt->flag = ncc;
      if ( !MG_EOK(pt) )  continue;

      adja = &mesh->adjt[3*(k-1)+1];
      for (i=0; i<3; i++) {
        if( ((pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY)) ||
            (pt->tag[i] & MG_BDY) ) continue;
        i1  = MMG5_inxt2[i];
        i2  = MMG5_iprv2[i];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];

        if ( !mesh->point[ip1].tmp )  mesh->point[ip1].tmp = ++nvf;
        if ( !mesh->point[ip2].tmp )  mesh->point[ip2].tmp = ++nvf;

        if ( MG_EDG(pt->tag[i]) || pt->tag[i] & MG_REQ ) {
          tag = mesh->point[ip1].tag;
          mesh->point[ip1].tag |= pt->tag[i];
          // Remove the MG_NOSURF tag if the vertex is really required.
          if ( (tag & MG_REQ) && !(tag & MG_NOSURF) ) {
            mesh->point[ip1].tag &= ~MG_NOSURF;
          }
          tag = mesh->point[ip2].tag;
          mesh->point[ip2].tag |= pt->tag[i];
          // Remove the MG_NOSURF tag if the vertex is really required.
          if ( (tag & MG_REQ) && !(tag & MG_NOSURF) ) {
            mesh->point[ip2].tag &= ~MG_NOSURF;
          }
        }

        /* open boundary */
        tag = 0;
        if ( mesh->info.opnbdy ) tag += MG_OPNBDY;
        if ( !adja[i] ) {
          /* Mark non-manifold edges and open-boundary ones: note that in open
           * boundary mode, all non-manifold edges are marked as open (because
           * at least one of the triangles that shares the edge has no
           * adjacent). */
          tag += MG_NOM;
          pt->tag[i] |= tag;
          mesh->point[ip1].tag |= tag;
          mesh->point[ip2].tag |= tag;
          ned++;
          continue;
        }

        kk = adja[i] / 3;
        ii = adja[i] % 3;
        if ( kk > k )  ned++;

        /* store adjacent */
        pt1 = &mesh->tria[kk];

        /* non manifold edge */
        if ( pt->tag[i] & MG_NOM ) {
          mesh->point[ip1].tag |= MG_NOM;
          mesh->point[ip2].tag |= MG_NOM;
          continue;
        }

        if ( MMG5_abs(pt1->ref) != MMG5_abs(pt->ref) ) {
          pt->tag[i]   |= MG_REF;
          pt1->tag[ii] |= MG_REF;
          mesh->point[ip1].tag |= MG_REF;
          mesh->point[ip2].tag |= MG_REF;
        }

        /* store adjacent */
        if ( !pt1->flag ) {
          pt1->flag    = ncc;
          pile[++ipil] = kk;
        }

        /* check orientation */
        ii1 = MMG5_inxt2[ii];
        ii2 = MMG5_iprv2[ii];
        if ( pt1->v[ii1] == ip1 ) {
          assert ( pt1->base );
          /* Moebius strip */
          if ( pt1->base < 0 ) {
            fprintf(stderr,"\n  ## Error: %s: Triangle orientation problem (1):"
                    " Moebius strip?\n",__func__);
            MMG5_SAFE_FREE(pile);
            return 0;
          }
          /* flip orientation */
          else {
            pt1->base   = -pt1->base;
            pt1->v[ii1] = ip2;
            pt1->v[ii2] = ip1;

            /* update adj */
            iad   = 3*(kk-1)+1;
            adjb  = &mesh->adjt[iad];
            adji1 = mesh->adjt[iad+ii1];
            adji2 = mesh->adjt[iad+ii2];
            adjb[ii1] = adji2;
            adjb[ii2] = adji1;

            /* modif tag + ref */
            tag = pt1->tag[ii1];
            pt1->tag[ii1] = pt1->tag[ii2];
            pt1->tag[ii2] = tag;
            ref = pt1->edg[ii1];
            pt1->edg[ii1] = pt1->edg[ii2];
            pt1->edg[ii2] = ref;

            /* modif voyeurs */
            if ( adjb[ii1] ) {
              iel = adjb[ii1] / 3;
              voy = adjb[ii1] % 3;
              mesh->adjt[3*(iel-1)+1+voy] = 3*kk + ii1;
            }
            if ( adjb[ii2] ) {
              iel = adjb[ii2] / 3;
              voy = adjb[ii2] % 3;
              mesh->adjt[3*(iel-1)+1+voy] = 3*kk + ii2;
            }
            nf++;
          }
        }
        else {
          /* Mark triangles that have a consistent orientation with their
           * neighbours */
          pt1->base =  -pt1->base;
        }
      }
    }
    while ( ipil > 0 );

    /* find next unmarked triangle */
    ipil = 0;
    for (kk=1; kk<=mesh->nt; kk++) {
      pt = &mesh->tria[kk];
      if ( MG_EOK(pt) && (pt->flag == 0) ) {
        ipil = 1;
        pile[ipil] = kk;
        pt->flag   = ncc+1;
        break;
      }
    }
  }

  /* bilan */
  nr = nm = nre = nreq = nt = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    nt++;
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( ( !MG_EDG(pt->tag[i]) ) && ( !(pt->tag[i] & MG_REQ) ) )  continue;

      jel  = adja[i] / 3;
      if ( !jel || jel > k ) {
        if ( pt->tag[i] & MG_GEO )  nr++;
        if ( pt->tag[i] & MG_NOM )  nm++;
        if ( pt->tag[i] & MG_REF )  nre++;
        if ( pt->tag[i] & MG_REQ )  nreq++;
      }
    }
  }

  if ( mesh->info.ddebug ) {
    fprintf(stdout,"  a- ridges: %" MMG5_PRId " found.\n",nr);
    fprintf(stdout,"  a- nm    : %" MMG5_PRId " found.\n",nm);
    fprintf(stdout,"  a- requir: %" MMG5_PRId " found.\n",nreq);
    fprintf(stdout,"  a- connex: %" MMG5_PRId " connected component(s)\n",ncc);
    fprintf(stdout,"  a- orient: %" MMG5_PRId " flipped\n",nf);
  }
  else if ( abs(mesh->info.imprim) > 3 ) {
    gen = (2 - nvf + ned - nt) / 2;
    fprintf(stdout,"     Connected component: %" MMG5_PRId ",  genus: %" MMG5_PRId ",   reoriented: %" MMG5_PRId "\n",ncc,gen,nf);
    fprintf(stdout,"     Edges: %" MMG5_PRId ",  tagged: %" MMG5_PRId ",  ridges: %" MMG5_PRId ", required: %" MMG5_PRId ", refs: %" MMG5_PRId "\n",
            ned,nr+nre+nreq,nr,nreq,nre);
  }

  MMG5_SAFE_FREE(pile);
  return 1;
}

/** check for ridges: dihedral angle */
int MMG5_setdhd(MMG5_pMesh mesh) {
  MMG5_pTria    pt,pt1;
  double        n1[3],n2[3],dhd;
  MMG5_int      *adja,k,kk,ne,nr,nrrm;
  int8_t        i,ii,i1,i2;
  static int8_t warn=0;

  /** Step 1: check input ridges provided by the user to remove those ones
   * between triangles belonging to the same plane. This step has to be done
   * prior the next one because we want to remove the MG_GEO tag transfered from
   * ridges that we delete toward vertices by \a MMG5_setadj but we want to
   * preserve the MG_GEO tags at vertices of ridges that will be added by the
   * next step. */
  nrrm = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    /* triangle normal */
    MMG5_nortri(mesh,pt,n1);
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( ((pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY)) ||
           (pt->tag[i] & MG_BDY) ) continue;

      if ( pt->tag[i] & MG_NOM ) {
        /* 1. We don't compute ridges along non-manifold edges because: if we
         * choose to analyze their angle, we have to check the normal deviation
         * during mesh adaptation (which is not done for now).
         *
         * 2. Do not analyze if nm edges are MG_REF ones because we can't
         * analyze if adjacent tria have same refs due to non-consistency of
         * adjacency building.
         */
        continue;
      }

      kk  = adja[i] / 3;
      ii  = adja[i] % 3;

      if ( kk && k < kk ) {
        pt1 = &mesh->tria[kk];
        /* check angle w. neighbor. */
        MMG5_nortri(mesh,pt1,n2);
        dhd = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];

        if ( (pt->tag[i] & MG_GEO) || (pt1->tag[ii] & MG_GEO) ) {
          /* Edge is provided as ridge by the user: check that it isn't at the
           * interface of 2 triangles belonging to the same plane */
          if ( fabs(dhd-1.) < MMG5_EPSOK ) {
            pt->tag[i]   &= ~MG_GEO;
            pt1->tag[ii] &= ~MG_GEO;
            i1 = MMG5_inxt2[i];
            i2 = MMG5_inxt2[i1];
            mesh->point[pt->v[i1]].tag &= ~MG_GEO;
            mesh->point[pt->v[i2]].tag &= ~MG_GEO;
            nrrm++;
            if ( !warn ) {
              fprintf(stdout,"\n  ## Warning: %s: at least one ridge along flat angle.\n"
                      "              Ridge tag will be removed from edges but vertices can"
                      " still have tags that interfer with remeshing.\n\n",
                      __func__);
              warn = 1;
            }
          }
        }
      }
    }
  }

  /** Step 2: check ref and angle with neighbour to update ref tags and
   * ridge ones */
  ne = nr = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    /* triangle normal */
    MMG5_nortri(mesh,pt,n1);
    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      if ( ((pt->tag[i] & MG_PARBDY) && !(pt->tag[i] & MG_PARBDYBDY)) ||
           (pt->tag[i] & MG_BDY) ) continue;

      if ( pt->tag[i] & MG_NOM ) {
        /* 1. We don't compute ridges along non-manifold edges because: if we
         * choose to analyze their angle, we have to check the normal deviation
         * during mesh adaptation (which is not done for now).
         *
         * 2. Do not analyze if nm edges are MG_REF ones because we can't
         * analyze if adjacent tria have same refs due to non-consistency of
         * adjacency building.
         */
        continue;
      }

      kk  = adja[i] / 3;
      ii  = adja[i] % 3;

      if ( kk && k < kk ) {
        pt1 = &mesh->tria[kk];
        /* reference curve */
        /* Remark: along non-manifold edges we store only adjacency relationship
         * between 2 surface parts (other are considered without adja). As we
         * don't ensure consistency in the choic of the surface we cannot rely
         * on the current test to detect ref edges. For this reason, all
         * non-manifold edges have to be marked as reference. */
        if ( pt1->ref != pt->ref ) {
          pt->tag[i]   |= MG_REF;
          pt1->tag[ii] |= MG_REF;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].tag |= MG_REF;
          mesh->point[pt->v[i2]].tag |= MG_REF;
          ne++;
        }

        /* check angle w. neighbor. */
        MMG5_nortri(mesh,pt1,n2);
        dhd = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
        if ( dhd <= mesh->info.dhd ) {
          pt->tag[i]   |= MG_GEO;
          pt1->tag[ii] |= MG_GEO;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].tag |= MG_GEO;
          mesh->point[pt->v[i2]].tag |= MG_GEO;
          nr++;
        }
        else if ( (pt->tag[i] & MG_GEO) || (pt1->tag[ii] & MG_GEO) ) {
          /* The MG_GEO tag of a vertex at interface of a "true" ridge and a
           * "spurious" ridge deleted at step 1 may have been removed by
           * error */
          pt->tag[i]   |= MG_GEO;
          pt1->tag[ii] |= MG_GEO;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].tag |= MG_GEO;
          mesh->point[pt->v[i2]].tag |= MG_GEO;
        }
      }
    }
  }
  if ( abs(mesh->info.imprim) > 3 && nr > 0 ) {
    fprintf(stdout,"     %" MMG5_PRId " ridges, %" MMG5_PRId " edges added\n",nr,ne);
  }
  if ( abs(mesh->info.imprim) > 3 && nrrm > 0 ) {
    fprintf(stdout,"     %" MMG5_PRId " ridges removed\n",nrrm);
  }


  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \return 1.
 *
 * check subdomains connected by a vertex and mark these vertex as CRN and REQ.
 *
 */
int MMG5_chkVertexConnectedDomains(MMG5_pMesh mesh){
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   ppt;
  MMG5_int      k,lists[MMG3D_LMAX+2];
  int64_t       listv[MMG3D_LMAX+2];
  int           ilists,ilistv;
  int           i0,ier;
  int8_t        i,j;
  static int8_t mmgWarn = 0;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->s = 0;
    ppt->flag = mesh->mark;
  }
  ++mesh->mark;

  /*count the number of tet around a point*/
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )   continue;

    for (i=0; i<4; i++) {
      mesh->point[pt->v[i]].s++;
    }
  }

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )   continue;

    /* point j on face i */
    for (i=0; i<4; i++) {
      for (j=0; j<3; j++) {
        if ( pt->xt ) {
          pxt = &mesh->xtetra[pt->xt];
        }
        else  pxt = 0;

        i0  = MMG5_idir[i][j];
        ppt = &mesh->point[pt->v[i0]];
        if ( !(ppt->tag & MG_BDY) ) continue;
        if ( ppt->flag == mesh->mark ) continue;

        /* Catch a boundary point by a boundary face */
        if ( (!pt->xt) || !(MG_BDY & pxt->ftag[i]) )  continue;
        if( ppt->tag & MG_NOM ){
          if ( mesh->adja[4*(k-1)+1+i] ) continue;
          ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,1);
        } else {
          ier=MMG5_boulesurfvolp(mesh,k,i0,i,listv,&ilistv,lists,&ilists,0);
        }
        if ( ier != 1 && !mmgWarn ) {
          mmgWarn = 1;
          printf("  ## Warning: %s: unable to check that we don't have"
                 " non-connected domains.\n",__func__);
        }

        if(ilistv != ppt->s) {
          if(!(ppt->tag & MG_REQ) ) {
            ppt->tag |= MG_REQ;
            ppt->tag |= MG_CRN;
          }
        }
        ppt->flag = mesh->mark;
      }
    }
  }
  return 1;
}

/** check for singularities */
int MMG5_singul(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p1,p2;
  double         ux,uy,uz,vx,vy,vz,dd;
  MMG5_int       list[MMG3D_LMAX+2],listref[MMG3D_LMAX+2],k,nc,nre;
  int            xp,nr,ns;
  int8_t         i;

  nre = nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) || ( ppt->tag & MG_NOM ) ||
          ( ppt->tag & MG_PARBDY ) ) continue;
      else if ( MG_EDG(ppt->tag) ) {
        /* Store the number of ridges passing through the point (xp) and the
         * number of ref edges (nr) */
        ns = MMG5_bouler(mesh,mesh->adjt,k,i,list,listref,&xp,&nr,MMG3D_LMAX);

        if ( !ns )  continue;
        if ( (xp+nr) > 2 ) {
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        else if ( (xp == 1) && (nr == 1) ) {
          ppt->tag |= MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
        }
        else if ( xp == 1 && !nr ){
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        else if ( nr == 1 && !xp ){
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        /* check ridge angle */
        else {
          p1 = &mesh->point[list[1]];
          p2 = &mesh->point[list[2]];
          ux = p1->c[0] - ppt->c[0];
          uy = p1->c[1] - ppt->c[1];
          uz = p1->c[2] - ppt->c[2];
          vx = p2->c[0] - ppt->c[0];
          vy = p2->c[1] - ppt->c[1];
          vz = p2->c[2] - ppt->c[2];
          dd = (ux*ux + uy*uy + uz*uz) * (vx*vx + vy*vy + vz*vz);
          if ( fabs(dd) > MMG5_EPSD ) {
            dd = (ux*vx + uy*vy + uz*vz) / sqrt(dd);
            if ( dd > -mesh->info.dhd ) {
              ppt->tag |= MG_CRN;
              nc++;
            }
          }
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 3 && nre > 0 )
    fprintf(stdout,"     %" MMG5_PRId " corners, %" MMG5_PRId " singular points detected\n",nc,nre);
  return 1;
}

/**
 * \param mesh pointer to mesh
 * \return 1 if successful, 0 if failed
 *
 * Compute normals at C1 vertices, for C0: tangents.
 *
 * Following list summerizes computed (and not computed) data depending on
 * point type:
 *  - Corner point (MG_CRN): nothing (and no xpoint);
 *  - Reference point (MG_REF): xp, tangent (ppt->n), 1 normal (pxp->n1);
 *  - Ridge point (MG_GEO): xp, tangent,2 normals (pxp->n1 and n2);
 *  - Non-manifold point (MG_NOM) are not filled here but in \a nmgeom function
 *  - Open boundary points (MG_OPNBDY) are non-manifold so they are filled in nmgeom too.
 *  - Required points are analyzed using their other flags (so they may have an
 * xpoint and a normal) but their normal is not used during remeshing.
 *
 * Normals at regular boundary points can be provided by users but are ignored
 * along featured edges.
 *
 */
int MMG5_norver(MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt;
  MMG5_xPoint    *pxp;
  double         n[3],dd;
  MMG5_int       *adja,k,kk,ng,nn,nt,nf,nnr;
  int8_t         i,ii,i1;

  /** recomputation of normals only if mesh->xpoint has been freed */
  if ( mesh->xpoint ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: %s: no research of boundary points"
              " and normals of mesh. mesh->xpoint must be freed to enforce"
              " analysis.\n",__func__);
    }
    return 1;
  }

  /** Step 1: identify boundary points */
  ++mesh->base;
  mesh->xp = 0;
  nnr      = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( ppt->flag == mesh->base )  continue;
      else {
        ++mesh->xp;
        ppt->flag = mesh->base;
        if ( mesh->nc1 ) {
          if ( ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2] > 0 ) {
            /** input normals are ignored along all type of featured edges (ref,
             * geo, nom) but it is possible to implement their taking into
             * account along non-ridges reference edges and external
             * non-manifold ones. */
            if ( ppt->tag & MG_PARBDY || ppt->tag & MG_CRN || MG_EDG_OR_NOM(ppt->tag) ) {
              ++nnr;
              continue;
            }
            ppt->xp = -1;
          }
        }
      }
    }
  }

  /** Step 2: Allocate memory to store normals for boundary points */
  mesh->xpmax  = MG_MAX( (MMG5_int)(1.5*mesh->xp),mesh->npmax);

  MMG5_ADD_MEM(mesh,(mesh->xpmax+1)*sizeof(MMG5_xPoint),"boundary points",return 0);
  MMG5_SAFE_CALLOC(mesh->xpoint,mesh->xpmax+1,MMG5_xPoint,return 0);

  /** Step 3: compute normals + tangents */
  nn = ng = nt = nf = 0;
  mesh->xp = 0;
  ++mesh->base;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];

      if ( ppt->tag & MG_PARBDY || ppt->tag & MG_CRN || ppt->tag & MG_NOM || ppt->flag == mesh->base )
      {
        continue;
      }

      /** At C1 point */
      if ( !MG_EDG(ppt->tag) ) {

        if ( (!mesh->nc1) ||
             ppt->n[0]*ppt->n[0]+ppt->n[1]*ppt->n[1]+ppt->n[2]*ppt->n[2]<=MMG5_EPSD2 ) {
          if ( !MMG5_boulen(mesh,mesh->adjt,k,i,ppt->n) ) {
            ++nf;
            continue;
          }
          else ++nn;
        }

        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                             "larger xpoint table",
                             mesh->xp--;return 0;);
        }
        ppt->xp = mesh->xp;
        pxp = &mesh->xpoint[ppt->xp];
        memcpy(pxp->n1,ppt->n,3*sizeof(double));
        ppt->n[0] = ppt->n[1] = ppt->n[2] = 0.;
        ppt->flag = mesh->base;

      }

      /** along ridge-curve */
      i1  = MMG5_inxt2[i];
      if ( !MG_EDG(pt->tag[i1]) ) {
        continue;
      }

      /* As we skip non-manifold point, the edge should be manifold */
      assert ( (!(MG_NOM & pt->tag[i1])) && "Unexpected non-manifold edge" );
      if ( !MMG5_boulen(mesh,mesh->adjt,k,i,n) ) {
        ++nf;
        continue;
      }

      ++mesh->xp;
      if(mesh->xp > mesh->xpmax){
        MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                           "larger xpoint table",
                           mesh->xp--;return 0;);
      }
      ppt->xp = mesh->xp;
      pxp = &mesh->xpoint[ppt->xp];
      memcpy(pxp->n1,n,3*sizeof(double));

      if ( (pt->tag[i1] & MG_GEO) && adja[i1] > 0 ) {
        kk = adja[i1] / 3;
        ii = adja[i1] % 3;
        ii = MMG5_inxt2[ii];
        if ( !MMG5_boulen(mesh,mesh->adjt,kk,ii,n) ) {
          ++nf;
          continue;
        }
        memcpy(pxp->n2,n,3*sizeof(double));

        /** Along ridge: compute tangent as intersection of n1 + n2 */
        ppt->n[0] = pxp->n1[1]*pxp->n2[2] - pxp->n1[2]*pxp->n2[1];
        ppt->n[1] = pxp->n1[2]*pxp->n2[0] - pxp->n1[0]*pxp->n2[2];
        ppt->n[2] = pxp->n1[0]*pxp->n2[1] - pxp->n1[1]*pxp->n2[0];
        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
        if ( dd > MMG5_EPSD2 ) {
          dd = 1.0 / sqrt(dd);
          ppt->n[0] *= dd;
          ppt->n[1] *= dd;
          ppt->n[2] *= dd;
        }
        ppt->flag = mesh->base;
        ++nt;
        continue;
      }

      /* compute tgte */
      ppt->flag = mesh->base;
      ++nt;
      if ( !MMG5_boulec(mesh,mesh->adjt,k,i,ppt->n) ) {
        ++nf;
        continue;
      }
      dd = pxp->n1[0]*ppt->n[0] + pxp->n1[1]*ppt->n[1] + pxp->n1[2]*ppt->n[2];
      ppt->n[0] -= dd*pxp->n1[0];
      ppt->n[1] -= dd*pxp->n1[1];
      ppt->n[2] -= dd*pxp->n1[2];
      dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        ppt->n[0] *= dd;
        ppt->n[1] *= dd;
        ppt->n[2] *= dd;
      }
    }
  }
  mesh->nc1 = 0;

  if ( abs(mesh->info.imprim) > 3 && nn+nt > 0 ) {
    if ( nnr )
      fprintf(stdout,"     %" MMG5_PRId " input normals ignored\n",nnr);
    fprintf(stdout,"     %" MMG5_PRId " normals,  %" MMG5_PRId " tangents updated  (%" MMG5_PRId " failed)\n",nn,nt,nf);
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param pt pointer to current triangle
 * \param k number of current point
 * \param c newly computed coordinates (giving negative area)
* \param n normal of triangle before regularization
 *
 * \return 0 if fail, 1 if success
 *
 * In coordinate regularization, performs a dichotomy between previous point /
 * and newly computed point in the case of negative area of a triangle
 *
 */
static inline int MMG3D_dichotomytria(MMG5_pMesh mesh, MMG5_pTria pt, MMG5_int k, double *c, double *n) {

  MMG5_pPoint  ppt;
  double       to,tp,t,nnew[3],p[3],o[3],result;
  int          it,maxit,pos,ier;

  it = 0;
  maxit = 5;
  to = 0.0;
  tp = 1.0;
  t = 0.5;
  pos = 0;

  ppt = &mesh->point[k];

  /* initial coordinates of point before regularization */
  o[0] = ppt->c[0];
  o[1] = ppt->c[1];
  o[2] = ppt->c[2];

  /* initial coordinates of new point */
  p[0] = c[0];
  p[1] = c[1];
  p[2] = c[2];

  do {
    mesh->point[0].c[0] = o[0] + t*(p[0] - o[0]);
    mesh->point[0].c[1] = o[1] + t*(p[1] - o[1]);
    mesh->point[0].c[2] = o[2] + t*(p[2] - o[2]);

    MMG5_nortri(mesh, pt, nnew);
    MMG5_dotprod(3,n,nnew,&result);

    if ( result <= 0.0 ) {
      tp = t;
    }
    else {
      c[0] = mesh->point[0].c[0];
      c[1] = mesh->point[0].c[1];
      c[2] = mesh->point[0].c[2];
      to = t;
      pos = 1;
    }

    t = 0.5*(to + tp);
  }
  while ( ++it < maxit );

  if ( pos ) {
    return 1;
  }
  else
    return 0;
}

/**
 * \param mesh pointer to the mesh
 * \param v list of vertices of current tetrahedron
 * \param k number of current point
 * \param c input : newly computed coordinates (giving negative area), output : coordinates after dichotomy
 *
 * \return 0 if fail, 1 if success
 *
 * In coordinate regularization, performs a dichotomy between previous point /
 * and newly computed point in the case of negative volume
 *
 */
static inline int MMG3D_dichotomytetra(MMG5_pMesh mesh, MMG5_int *v, MMG5_int k, double *c) {

  MMG5_pPoint  ppt;
  double       p[3],o[3],vol,to,tp,t;
  int          it,maxit,pos;

  it = 0;
  maxit = 5;
  to = 0.0;
  tp = 1.0;
  t = 0.5;
  pos = 0;

  ppt = &mesh->point[k];

  /* initial coordinates of point before regularization */
  o[0] = ppt->c[0];
  o[1] = ppt->c[1];
  o[2] = ppt->c[2];

  /* initial coordinates of new point */
  p[0] = c[0];
  p[1] = c[1];
  p[2] = c[2];

  do {
    mesh->point[0].c[0] = o[0] + t*(p[0] - o[0]);
    mesh->point[0].c[1] = o[1] + t*(p[1] - o[1]);
    mesh->point[0].c[2] = o[2] + t*(p[2] - o[2]);

    vol = MMG5_orvol(mesh->point,v);

    if ( vol <= 0.0 ) {
      tp = t;
    }
    else {
      c[0] = mesh->point[0].c[0];
      c[1] = mesh->point[0].c[1];
      c[2] = mesh->point[0].c[2];
      to = t;
      pos = 1;
    }

    t = 0.5*(to + tp);
  }
  while ( ++it < maxit );

  if ( pos ) {
    return 1;
  }
  else
    return 0;
}

/**
 * \param mesh pointer to a MMG5 mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Regularization procedure for vertices coordinates, dual Laplacian in 3D
 *
 */
int MMG3D_regver(MMG5_pMesh mesh) {
  MMG5_pTria    pt;
  MMG5_pTetra   ptet;
  MMG5_pPoint   ppt,p0;
  MMG5_Tria     tnew;
  double        *tabl,c[3],n[3],nnew[3],*cptr,lm1,lm2,cx,cy,cz,res0,res,result;
  int           i,ii,it,nit,ilist,noupdate;
  MMG5_int      k,kt,nn,iel,list[MMG5_LMAX],tlist[MMG5_LMAX],*adja,iad,v[4];
  int64_t       tetlist[MMG5_LMAX];

  /* assign seed to vertex */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !ppt->s )  ppt->s = k;
    }
  }

   for (k=1; k<=mesh->np; k++) {
     ppt = &mesh->point[k];
     ppt->flag = 0;
  }

  for (k=1; k<=mesh->ne; k++) {
    ptet = &mesh->tetra[k];
    if ( !MG_EOK(ptet) )  continue;
    for (i=0; i<4; i++) {
      ppt = &mesh->point[ptet->v[i]];
      if ( !ppt->flag )  ppt->flag = k;
    }
  }

  /* allocate memory for coordinates */
  MMG5_SAFE_CALLOC(tabl,3*mesh->np+1,double,return 0);

  /* Pointer toward the suitable adjacency array for Mmgs and Mmg3d */
  adja = mesh->adjt;

  it   = 0;
  nit  = 10;
  res0 = 0.0;
  lm1  = mesh->info.lxreg;
  lm2  = 0.99*lm1;
  while ( it++ < nit ) {
    /* step 1: laplacian */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];

      iad = 3*(k-1)+1;
      tabl[iad+0] = ppt->c[0];
      tabl[iad+1] = ppt->c[1];
      tabl[iad+2] = ppt->c[2];
      if ( !MG_VOK(ppt) )  continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM || MG_EDG(ppt->tag) ) continue;

      iel = ppt->s;
      if ( !iel ) continue; // Mmg3d

      pt = &mesh->tria[iel];
      i  = 0;
      if ( pt->v[1] == k )  i = 1;
      else if ( pt->v[2] == k ) i = 2;

      ilist = MMG5_boulep(mesh,iel,i,adja,list,tlist);

      /* average coordinates */
      cx = cy = cz = 0.0;
      for (i=1; i<=ilist; i++) {
        p0  = &mesh->point[list[i]];

        cptr = p0->c;
        cx += cptr[0];
        cy += cptr[1];
        cz += cptr[2];
      }
      cx /= ilist;
      cy /= ilist;
      cz /= ilist;

      /* Laplacian */
      cptr = ppt->c;
      tabl[iad+0] = cptr[0] + lm1 * (cx - cptr[0]);
      tabl[iad+1] = cptr[1] + lm1 * (cy - cptr[1]);
      tabl[iad+2] = cptr[2] + lm1 * (cz - cptr[2]);
    }

    /* step 2: anti-laplacian */
    res = 0;
    nn  = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];

      if ( !MG_VOK(ppt) )  continue;
      if ( MG_SIN(ppt->tag) || ppt->tag & MG_NOM || MG_EDG(ppt->tag) ) continue;

      iel = ppt->s;
      if ( !iel ) continue; // Mmg3d

      pt = &mesh->tria[iel];
      i = 0;
      if ( pt->v[1] == k )  i = 1;
      else if ( pt->v[2] == k ) i = 2;

      ilist = MMG5_boulep(mesh,iel,i,adja,list,tlist);

      /* average normal */
      cx = cy = cz = 0.0;
      for (i=1; i<=ilist; i++) {
        iad = 3*(list[i]-1) + 1;
        cx += tabl[iad+0];
        cy += tabl[iad+1];
        cz += tabl[iad+2];
      }
      cx /= ilist;
      cy /= ilist;
      cz /= ilist;

      /* antiLaplacian */
      iad = 3*(k-1)+1;
      c[0] = tabl[iad+0] - lm2 * (cx - tabl[iad+0]);
      c[1] = tabl[iad+1] - lm2 * (cy - tabl[iad+1]);
      c[2] = tabl[iad+2] - lm2 * (cz - tabl[iad+2]);

      cptr = ppt->c;

      mesh->point[0].c[0] = c[0];
      mesh->point[0].c[1] = c[1];
      mesh->point[0].c[2] = c[2];

      /* check for negative areas */
      noupdate = 0;
      for (kt = 0 ; kt<ilist ; kt++) {
        pt = &mesh->tria[tlist[kt]];

        if ( !MG_EOK(pt) ) continue;

        MMG5_nortri(mesh, pt, n);

        for (i=0;i<3;i++) {
          tnew.v[i] = pt->v[i];
        }

        i = 0;
        if ( pt->v[1] == k ) i = 1;
        if ( pt->v[2] == k ) i = 2;

        tnew.v[i] = 0;

        MMG5_nortri(mesh, &tnew, nnew);
        MMG5_dotprod(3,n,nnew,&result);
        if ( result < 0.0 ) {
          if (!MMG3D_dichotomytria(mesh,&tnew,k,c,n))
            noupdate = 1;
          continue;
        }
      }

      /* check for negative volumes */
      if ( !noupdate ) {
        iel = ppt->flag;
        ptet = &mesh->tetra[iel];

        if ( !MG_EOK(ptet) ) continue;

        i = 0;
        if ( ptet->v[1] == k )  i = 1;
        else if ( ptet->v[2] == k ) i = 2;
        else if ( ptet->v[3] == k ) i = 3;
        ilist = MMG5_boulevolp(mesh, iel, i, tetlist);

        for( kt=0 ; kt<ilist ; kt++ ) {
          iel = tetlist[kt] / 4;
          i   = tetlist[kt] % 4;
          ptet = &mesh->tetra[iel];

          for ( ii=0 ; ii<4 ; ii++ ) {
            v[ii] = ptet->v[ii];
          }
          v[i] = 0;
          result = MMG5_orvol(mesh->point,v);

          if ( result <= 0.0 ) {
            if (!MMG3D_dichotomytetra(mesh,v,k,c))
              noupdate = 1;
            continue;
          }
        }
      }
      if ( !noupdate ) {
        res += (cptr[0]-c[0])*(cptr[0]-c[0]) + (cptr[1]-c[1])*(cptr[1]-c[1]) + (cptr[2]-c[2])*(cptr[2]-c[2]);
        cptr[0] = c[0];
        cptr[1] = c[1];
        cptr[2] = c[2];
        nn++;
      }
    }
      if ( it == 1 )  res0 = res;
      if ( res0 > MMG5_EPSD )  res  = res / res0;
      if ( mesh->info.imprim < -1 || mesh->info.ddebug ) {
        fprintf(stdout,"     iter %5d  res %.3E\r",it,res);
        fflush(stdout);
      }
      if ( it > 1 && res < MMG5_EPS )  break;
    }
  /* reset the ppt->s and ppt->flag tags */
  for (k=1; k<=mesh->np; ++k) {
    mesh->point[k].s    = 0;
    mesh->point[k].flag = 0;
  }

  if ( mesh->info.imprim < -1 || mesh->info.ddebug )  fprintf(stdout,"\n");

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"     %" MMG5_PRId " coordinates regularized: %.3e\n",nn,res);

  MMG5_SAFE_FREE(tabl);
  return 1;
}


/**
 * \param mesh pointer to the mesh
 *
 * \return 0 if fail, 1 otherwise
 *
 * Define continuous geometric support at non manifold vertices, using volume
 * information.
 *
 * Following list summerizes computed data depending on point type:
 *  - Non-manifold external point along connex mesh part (MG_NOM): xp, tangent, normal n1, nnor=0
 *  - Non-manifold external point along edge or point connected mesh part (MG_NOM): no xp,no tangent,no normal,nnor=0
 *  - Non-manifold not required internal point (MG_NOM+MG_REQ): xp, tangent, no normal, nnor=1.
 *  - Non-manifold required internal point (MG_NOM+MG_REQ): no xp, no tangent, no normal, nnor=1.
 *  - Non-manifold open-boundary (MG_OPNBDY) edges are treated like others:
 *    - if we are along a "truly" open boundary, it is an internal nm edge so
 * we don't have any normal;
 *    - if the edge is marked as open while being only non-manifold, we may
 * compute a normal (along an external boundary) or not (along an internal one).
 *  - We don't compute anything along corner points.
 */
int MMG3D_nmgeom(MMG5_pMesh mesh){
  MMG5_pTetra     pt;
  MMG5_pPoint     p0;
  MMG5_pxPoint    pxp;
  MMG5_int        k;
  MMG5_int        *adja,base;
  double          n[3],t[3];
  int8_t          i,j,ip,ier;

  base = ++mesh->base;
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if( !MG_EOK(pt) ) continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      if ( adja[i] ) continue;
      for (j=0; j<3; j++) {
        ip = MMG5_idir[i][j];
        p0 = &mesh->point[pt->v[ip]];
        if ( p0->flag == base )  continue;
        else if ( (!(p0->tag & MG_NOM)) || (p0->tag & MG_PARBDY) ) {
          continue;
        }

        p0->flag = base;
        ier = MMG5_boulenm(mesh,k,ip,i,n,t);

        if ( ier < 0 )
          return 0;
        else if ( !ier ) {
          p0->tag |= MG_REQ;
          p0->tag &= ~MG_NOSURF;
        }
        else {
          if ( !p0->xp ) {
            ++mesh->xp;
            if(mesh->xp > mesh->xpmax){
              MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                                 "larger xpoint table",
                                 mesh->xp--;
                                 fprintf(stderr,"  Exit program.\n");return 0;);
            }
            p0->xp = mesh->xp;
          }
          pxp = &mesh->xpoint[p0->xp];
          memcpy(pxp->n1,n,3*sizeof(double));
          memcpy(p0->n,t,3*sizeof(double));
        }
      }
    }
  }
  /* Deal with the non-manifold points that do not belong to a surface
   * tetra (a tetra that has a face without adjacent)*/
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if( !MG_EOK(pt) ) continue;

    for (i=0; i<4; i++) {
      p0 = &mesh->point[pt->v[i]];
      if ( p0->tag & MG_REQ || !(p0->tag & MG_NOM) || p0->xp || (p0->tag & MG_PARBDY) ) {
        /* CRN points are skipped because we skip REQ points.
         * For connex meshes it is possible
         * to compute the tangent at required points but it is not if the mesh
         * contains some edge or point-connections: all points along such a feature line are marked
         * as required by the previous loop because \ref MMB5_boulnm fails, thus we pass here
         * even if the line belongs to a surface. Then \ref MMG5_boulenmInt fails too and
         * we end up without normals and tangents. For sake of simplicity and because normals
         * and tangents at required points are not used for remahsing, we choose to skip all
         * internal REQ points. */
        continue;
      }
      ier = MMG5_boulenmInt(mesh,k,i,t);
      if ( ier ) {
        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                            "larger xpoint table",
                            mesh->xp--;
                            fprintf(stderr,"  Exit program.\n");return 0;);
        }
        p0->xp = mesh->xp;
        pxp = &mesh->xpoint[p0->xp];
        memcpy(p0->n,t,3*sizeof(double));
        pxp->nnor = 1;
      }
      else {
        p0->tag |= MG_REQ;
        p0->tag &= ~MG_NOSURF;
      }
    }
  }

  /*for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !(p0->tag & MG_NOM) || p0->xp ) continue;
    p0->tag |= MG_REQ;
    p0->tag &= ~MG_NOSURF;
  }*/
  
  return 1;
}

/**
 * \param mesh pointer to mesh
 * \return 1 if successful, 0 if fail
 *
 * preprocessing stage: mesh analysis.
 *
 * At the end of this function:
 *   1. triangles have been deleted and stored in xtetra;
 *   2. Boundary point have been analyzed and associated to xpoints (if needed).
 *      If computed tangent is stored in ppt->n and normals in pxp->n1 and pxp->n2.
 *      Following list summerizes computed (and not computed) data depending on
 *      point type:
 *         - Corner point (MG_CRN): no computation (and no xpoint);
 *         - Required points have different treatment depending if they are along
 * a 'regular' surf point, along an internal or external non-manifold edge, along
 * a non connex mesh part(so they may have a xpoint, a tangent and a normal or not)
 * but their normal is not used during remeshing. See comments in \ref MMG3D_nmgeom
 *         - Reference point (MG_REF): xpoint, tangent, 1 normal;
 *         - Ridge point (MG_GEO): xp, tangent,2 normals;
 *         - Non-manifold point (MG_NOM): xp, tangent, no normal if internal,
 *           1 normal if external
 *         - Open boundary points (MG_OPNBDY) are treated as nm points.
 *
 */
int MMG3D_analys(MMG5_pMesh mesh) {
  MMG5_Hash hash;
  int       ier;
  
  /**--- stage 1: data structures for surface */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** SURFACE ANALYSIS\n");

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
    return 0;
  }

  if ( mesh->info.iso && mesh->info.opnbdy ) {
    ier = MMG3D_update_xtetra ( mesh );
    if ( !ier ) {
      fprintf(stderr,"\n  ## Problem when updating the xtetra data after ls discretization."
              " Exit program.\n");
      return 0;
    }
  }

  /* create prism adjacency */
  if ( !MMG3D_hashPrism(mesh) ) {
    fprintf(stderr,"\n  ## Prism hashing problem. Exit program.\n");
    return 0;
  }

  /* compatibility triangle orientation w/r tetras */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }
  
  /* identify surface mesh */
  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }

  MMG5_freeXTets(mesh);
  MMG5_freeXPrisms(mesh);

  /* Set surface triangles to required in nosurf mode or for parallel boundaries */
  MMG3D_set_reqBoundaries(mesh);

  /* create surface adjacency */
  memset ( &hash, 0x0, sizeof(MMG5_Hash));
  if ( !MMG3D_hashTria(mesh,&hash) ) {
    MMG5_DEL_MEM(mesh,hash.item);
    fprintf(stderr,"\n  ## Hashing problem (2). Exit program.\n");
    return 0;
  }

  /* build hash table for geometric edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /**--- stage 2: surface analysis */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  /* identify connexity, flip orientation of faces if needed and transfer
   * triangle edge tags toward vertices. */
  if ( !MMG5_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* check for ridges */
  if ( mesh->info.dhd > MMG5_ANGLIM && !MMG5_setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* identify singularities */
  if ( !MMG5_singul(mesh) ) {
    fprintf(stderr,"\n  ## MMG5_Singularity problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

  /* regularize vertices coordinates*/
  if ( mesh->info.xreg && !MMG3D_regver(mesh) ) {
    fprintf(stderr,"\n  ## Coordinates regularization problem. Exit program.\n");
    return 0;
  }

  /* define (and regularize) normals */
  if ( !MMG5_norver(mesh) ) {
    fprintf(stderr,"\n  ## Normal problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }
  if ( mesh->info.nreg && !MMG5_regnor(mesh) ) {
    fprintf(stderr,"\n  ## Normal regularization problem. Exit program.\n");
    return 0;
  }

  /* set bdry entities to tetra and fill the orientation field */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* set non-manifold edges sharing non-intersecting multidomains as required */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** UPDATING TOPOLOGY AT NON-MANIFOLD POINTS\n");

  if ( !MMG5_setNmTag(mesh,&hash) ) {
    fprintf(stderr,"\n  ## Non-manifold topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* check subdomains connected by a vertex and mark these vertex as corner and
     required */
  MMG5_chkVertexConnectedDomains(mesh);

  /* build hash table for geometric edges */
  if ( !mesh->na && !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /* Update edges tags and references for xtetras */
  if ( !MMG5_bdryUpdate(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* define geometry for non manifold points */
  if ( !MMG3D_nmgeom(mesh) ) return 0;

#ifndef NDEBUG
  MMG3D_chkfacetags(mesh);
#endif

#ifdef USE_POINTMAP
  /* Initialize source point with input index */
  MMG5_int ip;
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].src = ip;
#endif

  /* release memory */
  MMG5_DEL_MEM(mesh,mesh->htab.geom);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  if ( mesh->nprism ) MMG5_DEL_MEM(mesh,mesh->adjapr);

  return 1;
}
