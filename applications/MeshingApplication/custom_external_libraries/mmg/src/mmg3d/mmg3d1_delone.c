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
 * \file mmg3d/mmg3d1_delone.c
 * \brief Perform volume and surface mesh adaptation in delaunay mode.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo Doxygen documentation
 *
 * Perform volume and surface mesh adaptation in delaunay mode (\a
 * PATTERN preprocessor flag set to OFF).
 *
 * \todo Clean the boucle for (code copy...)
 */
#include "mmg3d.h"

#ifndef PATTERN

char  ddb;

#define _MMG5_LOPTL_MMG5_DEL     1.41
#define _MMG5_LOPTS_MMG5_DEL     0.6

/* Decomment this part to debug */
//int MMG_npuiss,MMG_nvol,MMG_npres,MMG_npd;


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param ne number of elements.
 * \param ifilt pointer to store the number of vertices filtered by the octree.
 * \param ns pointer to store the number of vertices insertions.
 * \param nc pointer to store the number of collapse.
 * \param warn pointer to store a flag that warn the user in case of
 * reallocation difficulty.
 * \param it iteration index.
 * \return -1 if fail and we don't save the mesh, 0 if fail but we try to save
 * the mesh, 1 otherwise.
 *
 * \a adpsplcol loop: split edges longer than \ref _MMG5_LOPTL_MMG5_DEL and
 * collapse edges shorter than \ref _MMG5_LOPTS_MMG5_DEL.
 *
 */
static inline int
_MMG5_boucle_for(MMG5_pMesh mesh, MMG5_pSol met,_MMG3D_pOctree octree,int ne,
                 int* ifilt,int* ns,int* nc,int* warn,int it) {
  MMG5_pTetra     pt;
  MMG5_pxTetra    pxt;
  MMG5_Tria       ptt;
  MMG5_pPoint     p0,p1,ppt;
  MMG5_pxPoint    pxp;
  double     dd,len,lmax,o[3],to[3],no1[3],no2[3],v[3];
  int        k,ip,ip1,ip2,list[MMG3D_LMAX+2],ilist,lists[MMG3D_LMAX+2],ilists,ref;
  int16_t    tag;
  char       imax,j,i,i1,i2,ifa0,ifa1;
  int        lon,ret,ier;
  double     lmin,lfilt;
  int        imin,iq;
  int        ii;
  double     lmaxtet,lmintet;
  int        imaxtet,imintet,base;

  base  = ++mesh->mark;
  for (k=1; k<=ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt)  || (pt->tag & MG_REQ) )   continue;
    else if ( pt->mark < base-1 )  continue;
    pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

    /* 1) find longest and shortest edge  and try to manage it*/
    imax = -1; lmax = 0.0;
    imin = -1; lmin = DBL_MAX;
    for (ii=0; ii<6; ii++) {
      if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
      len = _MMG5_lenedg(mesh,met,ii,pt);

      if ( len > lmax ) {
        lmax = len;
        imax = ii;
      }
      if ( len < lmin ) {
        lmin = len;
        imin = ii;
      }
    }
    if ( imax==-1 && (mesh->info.ddebug || mesh->info.imprim > 5 ) )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
              __FILE__,__LINE__,k);
    if ( imin==-1  && (mesh->info.ddebug || mesh->info.imprim > 5 ) )
      fprintf(stdout,"%s:%d: Warning: all edges of tetra %d are boundary and required\n",
              __FILE__,__LINE__,k);

    if ( lmax >= _MMG5_LOPTL_MMG5_DEL )  {
      /* proceed edges according to lengths */
      ifa0 = _MMG5_ifar[imax][0];
      ifa1 = _MMG5_ifar[imax][1];
      i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
      j  = _MMG5_iarfinv[i][imax];
      i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
      i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
      ip1 = pt->v[i1];
      ip2 = pt->v[i2];
      p0  = &mesh->point[ip1];
      p1  = &mesh->point[ip2];

      /* Case of a boundary face */
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
        if ( !(MG_GET(pxt->ori,i)) ) continue;
        ref = pxt->edg[_MMG5_iarf[i][j]];
        tag = pxt->tag[_MMG5_iarf[i][j]];
        if ( tag & MG_REQ )  continue;
        tag |= MG_BDY;
        ilist = _MMG5_coquil(mesh,k,imax,list);
        if ( !ilist )  continue;
        else if ( ilist<0 ) return(-1);
        if ( tag & MG_NOM ){
          if( !_MMG5_BezierNom(mesh,ip1,ip2,0.5,o,no1,to) )
            continue;
          else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
            _MMG5_tet2tri(mesh,k,i,&ptt);
            _MMG5_nortri(mesh,&ptt,no1);
            if ( !MG_GET(pxt->ori,i) ) {
              no1[0] *= -1.0;
              no1[1] *= -1.0;
              no1[2] *= -1.0;
            }
          }
        }
        else if ( tag & MG_GEO ) {
          if ( !_MMG5_BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) )
            continue;
          if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
            _MMG5_tet2tri(mesh,k,i,&ptt);
            _MMG5_nortri(mesh,&ptt,no1);
            no2[0] = to[1]*no1[2] - to[2]*no1[1];
            no2[1] = to[2]*no1[0] - to[0]*no1[2];
            no2[2] = to[0]*no1[1] - to[1]*no1[0];
            dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
            if ( dd > _MMG5_EPSD2 ) {
              dd = 1.0 / sqrt(dd);
              no2[0] *= dd;
              no2[1] *= dd;
              no2[2] *= dd;
            }
          }
        }
        else if ( tag & MG_REF ) {
          if ( !_MMG5_BezierRef(mesh,ip1,ip2,0.5,o,no1,to) )
            goto collapse;
        }
        else {
          if ( !_MMG5_norface(mesh,k,i,v) )  goto collapse;
          if ( !_MMG5_BezierReg(mesh,ip1,ip2,0.5,v,o,no1) ) goto collapse;

        }
        ip = _MMG3D_newPt(mesh,o,tag);
        if ( !ip ) {
          /* reallocation of point table */
          _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                              *warn=1;
                              goto collapse,
                              o,tag);
        }
        if ( met->m ) {
          if ( _MMG5_intmet(mesh,met,k,imax,ip,0.5) <=0 ) {
            _MMG3D_delPt(mesh,ip);
            goto collapse;
          }
        }
        ier = _MMG3D_simbulgept(mesh,met,list,ilist,ip);
        if ( !ier ) {
          ier = _MMG3D_dichoto1b(mesh,met,list,ilist,ip);
        }

        ier = _MMG5_split1b(mesh,met,list,ilist,ip,1,1);

        /* if we realloc memory in _MMG5_split1b pt and pxt pointers are not valid */
        pt = &mesh->tetra[k];
        pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

        if ( ier < 0 ) {
          fprintf(stderr,"  ## Error: unable to split.\n");
          _MMG3D_delPt(mesh,ip);
          return(-1);
        }
        else if ( !ier ) {
          _MMG3D_delPt(mesh,ip);
          goto collapse;
        } else {
          (*ns)++;

          ppt = &mesh->point[ip];
          if ( MG_EDG(tag) || (tag & MG_NOM) )
            ppt->ref = ref;
          else
            ppt->ref = pxt->ref[i];
          ppt->tag = tag;

          pxp = &mesh->xpoint[ppt->xp];
          if ( tag & MG_NOM ){
            memcpy(pxp->n1,no1,3*sizeof(double));
            memcpy(ppt->n,to,3*sizeof(double));
          }
          else if ( tag & MG_GEO ) {
            memcpy(pxp->n1,no1,3*sizeof(double));
            memcpy(pxp->n2,no2,3*sizeof(double));
            memcpy(ppt->n,to,3*sizeof(double));
          }
          else if ( tag & MG_REF ) {
            memcpy(pxp->n1,no1,3*sizeof(double));
            memcpy(ppt->n,to,3*sizeof(double));
          }
          else
            memcpy(pxp->n1,no1,3*sizeof(double));
        }
        continue;
      }
      else if(pt->xt){
        if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) {
          continue;
        }
        ilist = _MMG5_coquil(mesh,k,imax,list);
        if ( !ilist )    continue;
        else if ( ilist<0 ) return(-1);
        o[0] = 0.5*(p0->c[0] + p1->c[0]);
        o[1] = 0.5*(p0->c[1] + p1->c[1]);
        o[2] = 0.5*(p0->c[2] + p1->c[2]);
        ip = _MMG3D_newPt(mesh,o,MG_NOTAG);

        if ( !ip )  {
          /* reallocation of point table */
          _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                              *warn=1;
                              goto collapse,
                              o,MG_NOTAG);
        }
        ppt = &mesh->point[ip];
        if ( met->m ) {
          if ( _MMG5_intmet(mesh,met,k,imax,ip,0.5)<=0 ) {
            _MMG3D_delPt(mesh,ip);
            goto collapse;
          }
        }
        ier = _MMG5_split1b(mesh,met,list,ilist,ip,1,1);
        if ( ier < 0 ) {
          fprintf(stderr,"  ## Error: unable to split.\n");
          _MMG3D_delPt(mesh,ip);
          return(-1);
        }
        else if ( !ier ) {
          _MMG3D_delPt(mesh,ip);
          goto collapse;
        }
        else {
          if ( octree ) {
            _MMG3D_addOctree(mesh,octree,ip);
          }

          (*ns)++;
          continue;
        }

        /* Case of an internal face */
      } else {
        ilist = _MMG5_coquil(mesh,k,imax,list);
        if ( !ilist )    continue;
        else if ( ilist<0 ) return(-1);
        else if(ilist%2) goto collapse; //bdry edge
        o[0] = 0.5*(p0->c[0] + p1->c[0]);
        o[1] = 0.5*(p0->c[1] + p1->c[1]);
        o[2] = 0.5*(p0->c[2] + p1->c[2]);
        ip = _MMG3D_newPt(mesh,o,MG_NOTAG);

        if ( !ip )  {
          /* reallocation of point table */
          _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                              *warn=1;
                              goto collapse,
                              o,MG_NOTAG);
        }
        ppt = &mesh->point[ip];
        if ( met->m ) {
          if ( _MMG5_intmet(mesh,met,k,imax,ip,0.5)<=0 ) {
            _MMG3D_delPt(mesh,ip);
            goto collapse;
          };
        }

        /* Delaunay */
        if ( lmax<1.6 ) {
          lfilt = 0.7;
        }
        else lfilt = 0.2;

        if ( octree && !_MMG3D_octreein(mesh,met,octree,ip,lfilt) ) {
          _MMG3D_delPt(mesh,ip);
          (*ifilt)++;
          goto collapse;
        } else {
          lon = _MMG5_cavity(mesh,met,k,ip,list,ilist/2);
          if ( lon < 1 ) {
            // MMG_npd++; // decomment to debug
            _MMG3D_delPt(mesh,ip);
            goto collapse;
          } else {
            ret = _MMG5_delone(mesh,met,ip,list,lon);
            if ( ret > 0 ) {
              if ( octree ) {
                _MMG3D_addOctree(mesh,octree,ip);
              }
              (*ns)++;
              continue;
            }
            else if ( ret == 0 ) {
              // MMG_npd++; // decomment to debug
              _MMG3D_delPt(mesh,ip);
              goto collapse;//continue;
            }
            else { /*allocation problem ==> saveMesh*/
              _MMG3D_delPt(mesh,ip);
              return(0);
            }
          }
        }
      }
    }
  collapse:
    if(lmin <= _MMG5_LOPTS_MMG5_DEL) {
      // Case of an internal tetra with 4 ridges vertices.
      if ( lmin == 0 ) continue;

      ifa0 = _MMG5_ifar[imin][0];
      ifa1 = _MMG5_ifar[imin][1];
      i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
      j  = _MMG5_iarfinv[i][imin];
      i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
      i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
      ip = pt->v[i1];
      iq = pt->v[i2];
      p0 = &mesh->point[ip];
      p1 = &mesh->point[iq];

      if ( (p0->tag > p1->tag) || (p0->tag & MG_REQ) )  continue;

      /* Case of a boundary face */
      ilist = 0;
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
        tag = pxt->tag[_MMG5_iarf[i][j]];
        if ( tag & MG_REQ )  continue;
        tag |= MG_BDY;
        if ( p0->tag > tag )   continue;
        if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;

        if (_MMG5_boulesurfvolp(mesh,k,i1,i,
                                list,&ilist,lists,&ilists,(p0->tag & MG_NOM)) < 0 )
          return(-1);

        ilist = _MMG5_chkcol_bdy(mesh,met,k,i,j,list,ilist,lists,ilists,2);
        if ( ilist > 0 ) {
          ier = _MMG5_colver(mesh,met,list,ilist,i2,2);

          if ( ier < 0 ) return(-1);
          else if(ier) {
            _MMG3D_delPt(mesh,ier);
            (*nc)++;
            continue;
          }
        }
        else if (ilist < 0 )  return(-1);
      }
      /* Case of an internal face */
      else {
        ilist = _MMG5_boulevolp(mesh,k,i1,list);

        if ( p0->tag & MG_BDY )  continue;
        ilist = _MMG5_chkcol_int(mesh,met,k,i,j,list,ilist,2);
        if ( ilist > 0 ) {
          ier = _MMG5_colver(mesh,met,list,ilist,i2,2);
          if ( ilist < 0 ) continue;
          if ( ier < 0 ) return(-1);
          else if(ier) {
            if ( octree )
              _MMG3D_delOctree(mesh, octree, ier);
            _MMG3D_delPt(mesh,ier);
            (*nc)++;
            continue;
          }
        }
        else if (ilist < 0 )  return(-1);
      }
    } //end if lmin < _MMG5_LOPTS_MMG5_DEL

    /*2) longest and shortest edges are stucked => try another edges*/
    imaxtet = imax;
    imintet = imin;
    lmaxtet = lmax;
    lmintet = lmin;
    assert(lmin);

    for (ii=0; ii<6; ii++) {
      if ( pt->xt && (pxt->tag[ii] & MG_REQ) )  continue;
      if ( (ii==imintet) && (lmintet < _MMG5_LOPTS_MMG5_DEL)) continue;
      if ( (ii==imaxtet) && (lmaxtet > _MMG5_LOPTL_MMG5_DEL) ) continue;

      ip1  = _MMG5_iare[ii][0];
      ip2  = _MMG5_iare[ii][1];
      len = _MMG5_lenedg(mesh,met,ii,pt);

      imax = ii;
      lmax = len;
      imin = ii;
      lmin = len;
      if ( lmax >= _MMG5_LOPTL_MMG5_DEL )  {
        /* proceed edges according to lengths */
        ifa0 = _MMG5_ifar[imax][0];
        ifa1 = _MMG5_ifar[imax][1];
        i  = (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
        j  = _MMG5_iarfinv[i][imax];
        i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
        i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
        ip1 = pt->v[i1];
        ip2 = pt->v[i2];
        p0  = &mesh->point[ip1];
        p1  = &mesh->point[ip2];

        /* Case of a boundary face */
        if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
          if ( !(MG_GET(pxt->ori,i)) ) continue;
          ref = pxt->edg[_MMG5_iarf[i][j]];
          tag = pxt->tag[_MMG5_iarf[i][j]];
          if ( tag & MG_REQ )  continue;
          tag |= MG_BDY;
          ilist = _MMG5_coquil(mesh,k,imax,list);
          if ( !ilist )  continue;
          else if ( ilist<0 ) return(-1);
          if ( tag & MG_NOM ){
            if( !_MMG5_BezierNom(mesh,ip1,ip2,0.5,o,no1,to) )
              continue;
            else if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
              _MMG5_tet2tri(mesh,k,i,&ptt);
              _MMG5_nortri(mesh,&ptt,no1);
              if ( !MG_GET(pxt->ori,i) ) {
                no1[0] *= -1.0;
                no1[1] *= -1.0;
                no1[2] *= -1.0;
              }
            }
          }
          else if ( tag & MG_GEO ) {
            if ( !_MMG5_BezierRidge(mesh,ip1,ip2,0.5,o,no1,no2,to) )
              continue;
            if ( MG_SIN(p0->tag) && MG_SIN(p1->tag) ) {
              _MMG5_tet2tri(mesh,k,i,&ptt);
              _MMG5_nortri(mesh,&ptt,no1);
              no2[0] = to[1]*no1[2] - to[2]*no1[1];
              no2[1] = to[2]*no1[0] - to[0]*no1[2];
              no2[2] = to[0]*no1[1] - to[1]*no1[0];
              dd = no2[0]*no2[0] + no2[1]*no2[1] + no2[2]*no2[2];
              if ( dd > _MMG5_EPSD2 ) {
                dd = 1.0 / sqrt(dd);
                no2[0] *= dd;
                no2[1] *= dd;
                no2[2] *= dd;
              }
            }
          }
          else if ( tag & MG_REF ) {
            if ( !_MMG5_BezierRef(mesh,ip1,ip2,0.5,o,no1,to) )
              goto collapse2;
          }
          else {
            if ( !_MMG5_norface(mesh,k,i,v) )  goto collapse2;
            if ( !_MMG5_BezierReg(mesh,ip1,ip2,0.5,v,o,no1) ) goto collapse2;

          }
          ip = _MMG3D_newPt(mesh,o,tag);
          if ( !ip ){
            /* reallocation of point table */
            _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                                *warn=1;
                                goto collapse2//break
                                ,o,tag);
          }
          if ( met->m ) {
            if ( _MMG5_intmet(mesh,met,k,imax,ip,0.5)<=0 ) {
              _MMG3D_delPt(mesh,ip);
              goto collapse2;
            }
          }
          ier = _MMG3D_simbulgept(mesh,met,list,ilist,ip);
          if ( !ier ) {
            ier = _MMG3D_dichoto1b(mesh,met,list,ilist,ip);
          }

          ier = _MMG5_split1b(mesh,met,list,ilist,ip,1,1);
          /* if we realloc memory in _MMG5_split1b pt and pxt pointers are not valid */
          pt = &mesh->tetra[k];
          pxt = pt->xt ? &mesh->xtetra[pt->xt] : 0;

          if ( ier < 0 ) {
            fprintf(stderr,"  ## Error: unable to split.\n");
            return(-1);
          }
          else if ( !ier ) {
            _MMG3D_delPt(mesh,ip);
            goto collapse2;//continue;
          } else {
            (*ns)++;
            //~ if ( octree )
              //~ _MMG3D_addOctree(mesh,octree,ip);

            ppt = &mesh->point[ip];

            if ( MG_EDG(tag) || (tag & MG_NOM) )
              ppt->ref = ref;
            else
              ppt->ref = pxt->ref[i];
            ppt->tag = tag;

            pxp = &mesh->xpoint[ppt->xp];
            if ( tag & MG_NOM ){
              memcpy(pxp->n1,no1,3*sizeof(double));
              memcpy(ppt->n,to,3*sizeof(double));
            }
            else if ( tag & MG_GEO ) {
              memcpy(pxp->n1,no1,3*sizeof(double));
              memcpy(pxp->n2,no2,3*sizeof(double));
              memcpy(ppt->n,to,3*sizeof(double));
            }
            else if ( tag & MG_REF ) {
              memcpy(pxp->n1,no1,3*sizeof(double));
              memcpy(ppt->n,to,3*sizeof(double));
            }
            else
              memcpy(pxp->n1,no1,3*sizeof(double));
          }
          break;//imax continue;
        }
        else if(pt->xt){
          if ( (p0->tag & MG_BDY) && (p1->tag & MG_BDY) ) {
            continue;
          }
          ilist = _MMG5_coquil(mesh,k,imax,list);
          if ( !ilist )    continue;
          else if ( ilist<0 ) return(-1);
          o[0] = 0.5*(p0->c[0] + p1->c[0]);
          o[1] = 0.5*(p0->c[1] + p1->c[1]);
          o[2] = 0.5*(p0->c[2] + p1->c[2]);
          ip = _MMG3D_newPt(mesh,o,MG_NOTAG);

          if ( !ip )  {
            /* reallocation of point table */
            _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                                *warn=1;
                                goto collapse2
                                ,o,MG_NOTAG);
          }
          ppt = &mesh->point[ip];
          if ( met->m ) {
            if ( _MMG5_intmet(mesh,met,k,imax,ip,0.5)<=0 ) {
              _MMG3D_delPt(mesh,ip);
              goto collapse2;
            }
          }
          ier = _MMG5_split1b(mesh,met,list,ilist,ip,1,1);
          if ( ier < 0 ) {
            fprintf(stderr,"  ## Error: unable to split.\n");
            _MMG3D_delPt(mesh,ip);
            return(-1);
          }
          else if ( !ier ) {
            _MMG3D_delPt(mesh,ip);
            goto collapse2;
          }
          else {
            if ( octree )
              _MMG3D_addOctree(mesh,octree,ip);
            (*ns)++;
            break;//imax continue;
          }

          /* Case of an internal face */
        } else {
          ilist = _MMG5_coquil(mesh,k,imax,list);
          if ( !ilist )    continue;
          else if ( ilist<0 ) return(-1);
          else if(ilist%2) goto collapse2; //bdry edge
          o[0] = 0.5*(p0->c[0] + p1->c[0]);
          o[1] = 0.5*(p0->c[1] + p1->c[1]);
          o[2] = 0.5*(p0->c[2] + p1->c[2]);
          ip = _MMG3D_newPt(mesh,o,MG_NOTAG);

          if ( !ip )  {
            /* reallocation of point table */
            _MMG5_POINT_REALLOC(mesh,met,ip,mesh->gap,
                                *warn=1;
                                goto collapse2,
                                o,MG_NOTAG);
          }
          ppt = &mesh->point[ip];
          if ( met->m ) {
            if ( _MMG5_intmet(mesh,met,k,imax,ip,0.5)<=0 ) {
              _MMG3D_delPt(mesh,ip);
              goto collapse2;
            }
          }

          if ( lmaxtet<1.6 ) {
            lfilt = 0.7;
          }
          else lfilt = 0.2;

          if (  /*it &&*/  octree && !_MMG3D_octreein(mesh,met,octree,ip,lfilt) ) {
            _MMG3D_delPt(mesh,ip);
            (*ifilt)++;
            goto collapse2;
          } else {
            lon = _MMG5_cavity(mesh,met,k,ip,list,ilist/2);
            if ( lon < 1 ) {
              // MMG_npd++; // decomment to debug
              _MMG3D_delPt(mesh,ip);
              goto collapse2;
            } else {
              ret = _MMG5_delone(mesh,met,ip,list,lon);
              if ( ret > 0 ) {
                if ( octree )
                  _MMG3D_addOctree(mesh,octree,ip);
                (*ns)++;
                break;//imax continue;
              }
              else if ( ret == 0 ) {
                // MMG_npd++; // decomment to debug
                _MMG3D_delPt(mesh,ip);
                goto collapse2;//continue;
              }
              else { /*allocation problem ==> savemesh*/
                _MMG3D_delPt(mesh,ip);
                return(0);
              }
            }
          }
        }
      }
    collapse2:
      if(lmin > _MMG5_LOPTS_MMG5_DEL) continue;
      // Case of an internal tetra with 4 ridges vertices.
      if ( lmin == 0 ) continue;

      ifa0 = _MMG5_ifar[imin][0];
      ifa1 = _MMG5_ifar[imin][1];
      i  =  (pt->xt && (pxt->ftag[ifa1] & MG_BDY)) ? ifa1 : ifa0;
      j  = _MMG5_iarfinv[i][imin];
      i1 = _MMG5_idir[i][_MMG5_inxt2[j]];
      i2 = _MMG5_idir[i][_MMG5_iprv2[j]];
      ip = pt->v[i1];
      iq = pt->v[i2];
      p0 = &mesh->point[ip];
      p1 = &mesh->point[iq];

      if ( (p0->tag > p1->tag) || (p0->tag & MG_REQ) )  continue;

      /* Case of a boundary face */
      ilist = 0;
      if ( pt->xt && (pxt->ftag[i] & MG_BDY) ) {
        tag = pxt->tag[_MMG5_iarf[i][j]];
        if ( tag & MG_REQ )  continue;
        tag |= MG_BDY;
        if ( p0->tag > tag )   continue;
        if ( ( tag & MG_NOM ) && (mesh->adja[4*(k-1)+1+i]) ) continue;

        if (_MMG5_boulesurfvolp(mesh,k,i1,i,
                                list,&ilist,lists,&ilists,(p0->tag & MG_NOM)) < 0 )
          return(-1);

        ilist = _MMG5_chkcol_bdy(mesh,met,k,i,j,list,ilist,lists,ilists,2);
        if ( ilist > 0 ) {
          ier = _MMG5_colver(mesh,met,list,ilist,i2,2);
          if ( ier < 0 ) return(-1);
          else if(ier) {
            _MMG3D_delPt(mesh,ier);
            (*nc)++;
            break;
          }
        }
        else if (ilist < 0 )  return(-1);
      }
      /* Case of an internal face */
      else {
        if ( p0->tag & MG_BDY )  continue;
        ilist = _MMG5_boulevolp(mesh,k,i1,list);
        ilist = _MMG5_chkcol_int(mesh,met,k,i,j,list,ilist,2);
        if ( ilist > 0 ) {
          ier = _MMG5_colver(mesh,met,list,ilist,i2,2);
          if ( ilist < 0 ) continue;
          if ( ier < 0 ) return(-1);
          else if(ier) {
            if ( octree )
              _MMG3D_delOctree(mesh,octree,ier);
            _MMG3D_delPt(mesh,ier);
            (*nc)++;
            break;
          }
        }
        else if (ilist < 0 )  return(-1);
      }
    }//end for ii
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \param warn set to 1 if we can't insert point due to lack of memory.
 * \return -1 if fail and we dont try to end the remesh process,
 * 0 if fail but we try to end the remesh process and 1 if success.
 *
 * Split edges longer than \ref _MMG5_LOPTL_MMG5_DEL and collapse edges shorter
 * than \ref _MMG5_LOPTS_MMG5_DEL.
 *
 */
static int
_MMG5_adpsplcol(MMG5_pMesh mesh,MMG5_pSol met,_MMG3D_pOctree octree, int* warn) {
  int        nfilt,ifilt,ne,ier;
  int        ns,nc,it,nnc,nns,nnf,nnm,maxit,nf,nm;
  double     maxgap;

  /* Iterative mesh modifications */
  it = nnc = nns = nnf = nnm = nfilt = 0;
  maxit = 10;
  mesh->gap = maxgap = 0.5;
  // MMG_npuiss = MMG_nvol = MMG_npres = MMG_npd = 0; // decomment to debug
  do {
    if ( !mesh->info.noinsert ) {
      *warn=0;
      ns = nc = 0;
      nf = nm = 0;
      ifilt = 0;
      ne = mesh->ne;
      ier = _MMG5_boucle_for(mesh,met,octree,ne,&ifilt,&ns,&nc,warn,it);
      if(ier<0) exit(EXIT_FAILURE);
      else if(!ier) return(-1);
    } /* End conditional loop on mesh->info.noinsert */
    else  ns = nc = ifilt = 0;

    if ( !mesh->info.noswap ) {
      nf = _MMG5_swpmsh(mesh,met,octree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;
      if(it==2 || it==6/*&& it==1 || it==3 || it==5 || it > 8*/) {
        nf += _MMG5_swptet(mesh,met,1.053,octree,2);
      } else {
        nf += 0;
      }
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      /*perform only boundary moves*/
      nm = _MMG5_movtet(mesh,met,octree,-1);
      if ( nm < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh.\n");
        return(0);
      }
    }
    else  nm = 0;

    nnm += nm;
    nnc += nc;
    nns += ns;
    nnf += nf;
    nfilt += ifilt;

    /* decrease size of gap for reallocation */

    if ( mesh->gap > maxgap/(double)maxit )
      mesh->gap -= maxgap/(double)maxit;
    else
      mesh->gap -= mesh->gap/(double)maxit;


    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && ns+nc+nm+nf > 0)
      fprintf(stdout,"     %8d filtered, %8d splitted, %8d collapsed,"
              " %8d swapped, %8d moved\n",ifilt,ns,nc,nf,nm);

    if ( ns < 10 && abs(nc-ns) < 3 )  break;
    else if ( it > 3 && abs(nc-ns) < 0.3 * MG_MAX(nc,ns) )  break;

  }
  while( ++it < maxit && nc+ns > 0 );

  if ( mesh->info.imprim ) {
    if ( (abs(mesh->info.imprim) < 5) && ( nnc || nns ) ) {
      fprintf(stdout,"     %8d filtered, %8d splitted, %8d collapsed,"
              " %8d swapped, %8d moved, %d iter.\n",nfilt,nns,nnc,nnf,nnm, it);
    }
  }

  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization for LES computation (improve the element skewness).
 *
 */
static int
_MMG5_optetLES(MMG5_pMesh mesh, MMG5_pSol met,_MMG3D_pOctree octree) {
  int it,nnm,nnf,maxit,nm,nf,nw;
  double declic;

  it = nnm = nnf = 0;
  maxit = 10;
  declic = 1.01;
  do {
    /* treatment of bad elements*/
    if(it < 5) {
      nw = MMG3D_opttyp(mesh,met,octree);
    }
    else
      nw = 0;
    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = _MMG5_swptet(mesh,met,declic,octree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      nm = _MMG5_movtet(mesh,met,octree, 3);
      if ( nm < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh.\n");
        return(0);
      }
    }
    else  nm = 0;
    nnm += nm;

//be careful, this procedure can degrade the worst elt
    if ( !mesh->info.nomove && (it==2)) {
      _MMG3D_optlap(mesh,met);
    }

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nw+nf+nm > 0 ){
      fprintf(stdout,"                                          ");
      fprintf(stdout,"  %8d improved, %8d swapped, %8d moved\n",nw,nf,nm);
    }
  }
  while( ++it < maxit && nw+nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = _MMG5_movtet(mesh,met,octree, 3);
    if ( nm < 0 ) {
      fprintf(stderr,"  ## Unable to improve mesh.\n");
      return(0);
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ) {
    fprintf(stdout,"                                            "
            "                                ");
    fprintf(stdout,"     %8d moved\n",nm);
  }


  if ( mesh->info.imprim ) {
    if ( abs(mesh->info.imprim) < 5 && (nnf > 0 || nnm > 0) )
      fprintf(stdout,"                                                 "
              "        "
              "      %8d swapped, %8d moved, %d iter. \n",nnf,nnm,it);
  }
  return(1);
}
/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Mesh optimization using egde swapping and point relocation.
 *
 */
static int
_MMG5_optet(MMG5_pMesh mesh, MMG5_pSol met,_MMG3D_pOctree octree) {
  int it,nnm,nnf,maxit,nm,nf,nw;
  double declic;

  /* shape optim */
  it = nnm = nnf = 0;
  maxit = 4;
  declic = 1.053;
  do {
    /* treatment of bad elements*/
    if(it < 5) {
      nw = MMG3D_opttyp(mesh,met,octree);
    }
    else
      nw = 0;
    /* badly shaped process */
    if ( !mesh->info.noswap ) {
      nf = _MMG5_swpmsh(mesh,met,octree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
      nnf += nf;

      nf = _MMG5_swptet(mesh,met,declic,octree,2);
      if ( nf < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
        return(0);
      }
    }
    else  nf = 0;

    if ( !mesh->info.nomove ) {
      nm = _MMG5_movtet(mesh,met,octree,0);
      if ( nm < 0 ) {
        fprintf(stderr,"  ## Unable to improve mesh.\n");
        return(0);
      }
    }
    else  nm = 0;
    nnm += nm;

    if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nw+nf+nm > 0 ){
      fprintf(stdout,"                                          ");
      fprintf(stdout,"  %8d improved, %8d swapped, %8d moved\n",nw,nf,nm);
    }
  }
  while( ++it < maxit && nw+nm+nf > 0 );

  if ( !mesh->info.nomove ) {
    nm = _MMG5_movtet(mesh,met,octree,3);
    if ( nm < 0 ) {
      fprintf(stderr,"  ## Unable to improve mesh.\n");
      return(0);
    }
  }
  else  nm = 0;
  nnm += nm;
  if ( (abs(mesh->info.imprim) > 4 || mesh->info.ddebug) && nm > 0 ) {
    fprintf(stdout,"                                            "
            "                                ");
    fprintf(stdout,"     %8d moved\n",nm);
  }

  if ( mesh->info.imprim ) {
    if ( abs(mesh->info.imprim) < 5 && (nnf > 0 || nnm > 0) )
      fprintf(stdout,"                                                 "
              "        "
              "      %8d swapped, %8d moved, %d iter. \n",nnf,nnm,it);
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param octree pointer toward the octree structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Analyze tetrahedra and split long / collapse short, according to
 * prescribed metric.
 *
 */
static int
_MMG5_adptet_delone(MMG5_pMesh mesh,MMG5_pSol met,_MMG3D_pOctree octree) {
  int      nnf,ns,nf;
  int      warn;

  /*initial swap*/
  if ( !mesh->info.noswap ) {
    nf = _MMG5_swpmsh(mesh,met,octree,2);
    if ( nf < 0 ) {
      fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf = nf;
    nf = _MMG5_swptet(mesh,met,1.053,octree,2);
    if ( nf < 0 ) {
      fprintf(stderr,"  ## Unable to improve mesh. Exiting.\n");
      return(0);
    }
    nnf+=nf;
  } else  nnf = nf = 0;

#ifdef DEBUG
  fprintf(stdout,"$$$$$$$$$$$$$$$$$$ INITIAL SWAP %7d\n",nnf);
  _MMG3D_outqua(mesh,met);
#endif

  /* Iterative mesh modifications */
  warn = 0;

  ns = _MMG5_adpsplcol(mesh,met,octree,&warn);

  if ( ns < 0 ) {
    fprintf(stderr,"  ## Unable to complete mesh. Exit program.\n");
    return(0);
  }

  if ( warn ) {
    fprintf(stderr,"  ## Error:");
    fprintf(stderr," unable to allocate a new point in last call of adpspl.\n");
    fprintf(stderr,"  ## Check the mesh size or ");
    fprintf(stderr,"increase the maximal authorized memory with the -m option.\n");
    fprintf(stderr,"  ## Uncomplete mesh. Exiting\n" );
    return(0);
  }

  /* renumerotation if available */
  if ( !_MMG5_scotchCall(mesh,met) )
    return(0);
  
  if(mesh->info.optimLES) {
    if(!_MMG5_optetLES(mesh,met,octree)) return(0);
  }
  else {
    if(!_MMG5_optet(mesh,met,octree)) return(0);
  }
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \return 0 if failed, 1 if success.
 *
 * Main adaptation routine.
 *
 */
int _MMG5_mmg3d1_delone(MMG5_pMesh mesh,MMG5_pSol met) {
  _MMG3D_pOctree octree = NULL;

  if ( abs(mesh->info.imprim) > 4 )
    fprintf(stdout,"  ** MESH ANALYSIS\n");

  if ( mesh->info.iso && !_MMG5_chkmani(mesh) ) {
    fprintf(stderr,"  ## Non orientable implicit surface. Exit program.\n");
    return(0);
  }

  /**--- stage 1: geometric mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** GEOMETRIC MESH\n");

  if ( !_MMG5_anatet(mesh,met,1,0) ) {
    fprintf(stderr,"  ## Unable to split mesh. Exiting.\n");
    if ( octree )
      /*free octree*/
      _MMG3D_freeOctree(mesh,&octree);
    return(0);
  }

#ifdef DEBUG
  _MMG3D_inqua(mesh,met0);
#endif

  /**--- stage 2: computational mesh */
  if ( abs(mesh->info.imprim) > 4 || mesh->info.ddebug )
    fprintf(stdout,"  ** COMPUTATIONAL MESH\n");

  /* define metric map */
  if ( !_MMG5_defsiz(mesh,met) ) {
    fprintf(stderr,"  ## Metric undefined. Exit program.\n");
    if ( octree )
      /*free octree*/
      _MMG3D_freeOctree(mesh,&octree);
    return(0);
  }

  if ( mesh->info.hgrad > 0. ) {
    if ( mesh->info.imprim )   fprintf(stdout,"\n  -- GRADATION : %8f\n",exp(mesh->info.hgrad));
    if ( !_MMG5_gradsiz(mesh,met) ) {
      fprintf(stderr,"  ## Gradation problem. Exit program.\n");
      if ( octree )
        /*free octree*/
        _MMG3D_freeOctree(mesh,&octree);
      return(0);
    }
  }

  if ( !_MMG5_anatet(mesh,met,2,0) ) {
    fprintf(stderr,"  ## Unable to split mesh. Exiting.\n");
    if ( octree )
      /*free octree*/
      _MMG3D_freeOctree(mesh,&octree);
    return(0);
  }

#ifdef DEBUG
  puts("---------------------------Fin anatet---------------------");
  _MMG3D_outqua(mesh,met);
#endif

  /* renumerotation if available */
  if ( !_MMG5_scotchCall(mesh,met) ) {
    if ( octree )
      /*free octree*/
      _MMG3D_freeOctree(mesh,&octree);
    return(0);
  }

  if ( mesh->info.octree > 0 ) {
    if ( !_MMG3D_initOctree(mesh,&octree,mesh->info.octree) ) {
      if ( octree )
        /*free octree*/
        _MMG3D_freeOctree(mesh,&octree);
    }
  }

  if ( !_MMG5_adptet_delone(mesh,met,octree) ) {
    fprintf(stderr,"  ## Unable to adapt. Exit program.\n");
    if ( octree )
      /*free octree*/
      _MMG3D_freeOctree(mesh,&octree);
    return(0);
  }

#ifdef DEBUG
  puts("---------------------Fin adptet-----------------");
  _MMG3D_outqua(mesh,met);
#endif
  /* in test phase: check if no element with 2 bdry faces */
  if ( !_MMG5_chkfemtopo(mesh) ) {
    fprintf(stderr,"  ## Topology of mesh unsuited for fem computations. Exit program.\n");
    if ( octree )
      /*free octree*/
      _MMG3D_freeOctree(mesh,&octree);
    return(0);
  }

  if ( mesh->info.iso && !_MMG5_chkmani(mesh) ) {
    fprintf(stderr,"  ## Non orientable implicit surface. Exit program.\n");
    if ( octree )
      /*free octree*/
      _MMG3D_freeOctree(mesh,&octree);
    return(0);
  }

  if ( octree )
    /*free octree*/
    _MMG3D_freeOctree(mesh,&octree);

  return(1);
}

#endif
