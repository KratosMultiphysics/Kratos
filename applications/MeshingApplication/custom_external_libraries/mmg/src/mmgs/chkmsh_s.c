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
 * \file mmgs/chkmsh_s.c
 * \brief Check the input mesh validity.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgs.h"


/**
 * \param mesh pointer toward the mesh structure.
 * \param severe level of performed check
 * \param base unused argument.
 * \return 0 if fail, 1 if success.
 *
 * Check the mesh validity
 *
 */
int MMG5_mmgsChkmsh(MMG5_pMesh mesh,int severe,int base) {
    MMG5_pPoint         ppt;
    MMG5_pTria          pt1,pt2;
    int                 adj,adj1,k,kk,l,nk,i,j,ip,lon,len;
    int                 *adja,*adjb,list[MMGS_LMAX+2];
    char                voy,voy1,i1,i2,j1,j2;
    static char         mmgErr0=0,mmgErr1=0,mmgErr2=0,mmgErr3=0,mmgErr4=0;
    static char         mmgErr5=0,mmgErr6=0,mmgErr7=0;

    for (k=1; k<=mesh->nt; k++) {
        pt1 = &mesh->tria[k];
        if ( !MG_EOK(pt1) )  continue;
        adja = &mesh->adja[3*(k-1)+1];

        for (i=0; i<3; i++) {
            if ( !adja[i] )  continue;
            i1  = MMG5_inxt2[i];
            i2  = MMG5_iprv2[i];
            adj = adja[i] / 3;
            voy = adja[i] % 3;
            if ( !adj && !(pt1->tag[i] & MG_GEO) ) {
              if ( !mmgErr0 ) {
                mmgErr0 = 1;
                fprintf(stderr,"  ## Error: %s: 0. at least 1 missing edge"
                        " tag (%d %d)\n",__func__,MMGS_indElt(mesh,k),
                        MMGS_indElt(mesh,adj));
                fprintf(stderr,"triangle %d: %d %d %d \n",MMGS_indElt(mesh,k),
                        MMGS_indPt(mesh,pt1->v[0]),MMGS_indPt(mesh,pt1->v[1]),
                        MMGS_indPt(mesh,pt1->v[2]));
                fprintf(stderr,"tag (%d): %d %d %d \n",MMGS_indElt(mesh,k),
                        pt1->tag[0],pt1->tag[1],pt1->tag[2]);
              }
              return 0;
            }
            if ( adj == k ) {
              if ( !mmgErr1 ) {
                mmgErr1 = 1;
                fprintf(stderr,"\n  ## Error: %s: 1. at least 1 wrong adjacency"
                        " (%d %d)\n",
                        __func__,MMGS_indElt(mesh,k),MMGS_indElt(mesh,adj));
                fprintf(stderr,"triangle %d: %d %d %d \n",MMGS_indElt(mesh,k),
                        MMGS_indPt(mesh,pt1->v[0]),MMGS_indPt(mesh,pt1->v[1]),
                        MMGS_indPt(mesh,pt1->v[2]));
                fprintf(stderr,"adj (%d): %d %d %d \n",MMGS_indElt(mesh,k),
                        MMGS_indElt(mesh,adja[0]/3),MMGS_indElt(mesh,adja[1]/3),
                        MMGS_indElt(mesh,adja[2]/3));
              }
              return 0;
            }
            pt2 = &mesh->tria[adj];
            if ( !MG_EOK(pt2) ) {
              if ( !mmgErr2 ) {
                mmgErr2 = 1;
                fprintf(stderr,"\n  ## Error: %s: 4. At least 1 invalid adjacent"
                        " (%d %d)\n",__func__,MMGS_indElt(mesh,adj),
                        MMGS_indElt(mesh,k));
                fprintf(stderr,"vertices of k   %d: %d %d %d\n",
                        MMGS_indElt(mesh,k),MMGS_indPt(mesh,pt1->v[0]),
                        MMGS_indPt(mesh,pt1->v[1]),MMGS_indPt(mesh,pt1->v[2]));
                fprintf(stderr,"vertices of adj %d: %d %d %d \n",
                        MMGS_indElt(mesh,adj),MMGS_indPt(mesh,pt2->v[0]),
                        MMGS_indPt(mesh,pt2->v[1]),MMGS_indPt(mesh,pt2->v[2]));
              }
              return 0;
            }
            if ( (pt1->tag[i] != pt2->tag[voy]) || (pt1->edg[i] != pt2->edg[voy] ) ) {
                fprintf(stderr,"\n  ## Error: %s: 3. at least 1 wrong"
                        " tag/ref (%d %d)"
                        "  %d - %d\n",__func__,MMGS_indElt(mesh,k),
                        MMGS_indElt(mesh,adj),pt1->tag[i],pt2->tag[voy]);
                return 0;
            }
            adjb = &mesh->adja[3*(adj-1)+1];
            adj1 = adjb[voy] / 3;
            voy1 = adjb[voy] % 3;
            if ( adj1 != k || voy1 != i ) {
              if ( !mmgErr3 ) {
                mmgErr3 = 1;
                fprintf(stderr,"\n  ## Error: %s: 2. at least 1 wrong"
                        " adjacency (%d %d)\n",__func__, MMGS_indElt(mesh,k),
                         MMGS_indElt(mesh,adj1));
                fprintf(stderr,"vertices of %d: %d %d %d \n",MMGS_indElt(mesh,k),
                        MMGS_indPt(mesh,pt1->v[0]),MMGS_indPt(mesh,pt1->v[1]),
                        MMGS_indPt(mesh,pt1->v[2]));
                fprintf(stderr,"vertices of adj %d: %d %d %d \n",
                        MMGS_indElt(mesh,adj),
                        MMGS_indPt(mesh,pt2->v[0]),MMGS_indPt(mesh,pt2->v[1]),
                        MMGS_indPt(mesh,pt2->v[2]));
                fprintf(stderr,"adj(%d): %d %d %d\n",MMGS_indElt(mesh,k),
                        MMGS_indElt(mesh,adja[0]/3),MMGS_indElt(mesh,adja[1]/3),
                        MMGS_indElt(mesh,adja[2]/3));
                fprintf(stderr,"adj(%d): %d %d %d\n",MMGS_indElt(mesh,adj),
                        MMGS_indElt(mesh,adjb[0]/3),MMGS_indElt(mesh,adjb[1]/3),
                        MMGS_indElt(mesh,adjb[2]/3));
              }
              return 0;
            }
            if ( !MS_SIN(pt1->tag[i]) ) {
                j1 = MMG5_inxt2[voy];
                j2 = MMG5_iprv2[voy];
                if ( pt2->v[j2] != pt1->v[i1] || pt2->v[j1] != pt1->v[i2] ) {
                  if ( !mmgErr4 ) {
                    mmgErr4 = 1;
                    fprintf(stderr,"\n  ## Error: %s: 8. at least 1 wrong"
                            " orientation (%d %d).\n",__func__,MMGS_indElt(mesh,k),
                            MMGS_indElt(mesh,adj));
                  }
                  return 0;
                }
            }
        }
    }

    if ( !severe )  return 1;

    for (k=1; k<=mesh->nt; k++) {
        pt1 = &mesh->tria[k];
        if ( !MG_EOK(pt1) )  continue;

        adja = &mesh->adja[3*(k-1)+1];
        for (i=0; i<3; i++) {
            if ( !adja[i] )  continue;

            ip  = pt1->v[i];
            ppt = &mesh->point[ip];
            if ( !MG_VOK(ppt) ) {
              if ( !mmgErr5 ) {
                mmgErr5 = 1;
                fprintf(stderr,"\n  ## Error: %s: 6. at least 1 unused"
                        " vertex (%d  %d)\n",__func__,
                        MMGS_indElt(mesh,k),MMGS_indPt(mesh,ip));
                fprintf(stderr,"%d %d %d\n",MMGS_indPt(mesh,pt1->v[0]),
                        MMGS_indPt(mesh,pt1->v[1]),MMGS_indPt(mesh,pt1->v[2]));
              }
              return 0;
            }
            else if ( MS_SIN(ppt->tag) )  continue;

            lon = boulet(mesh,k,i,list);
            if ( lon < 1 )  continue;
            for (l=0; l<lon; l++) {
                kk  = list[l] / 3;
                nk  = list[l] % 3;
                pt2 = &mesh->tria[kk];
                if ( pt2->v[nk] != ip ) {
                  if ( !mmgErr6 ) {
                    mmgErr6 = 1;
                    fprintf(stderr,"\n  ## Error: %s: 5. at least 1 wrong"
                            " ball (%d, %d).\n",__func__,MMGS_indPt(mesh,ip),
                            MMGS_indPt(mesh,pt2->v[nk]));
                  }
                  return 0;
                }
            }
            len = 0;
            for (kk=1; kk<=mesh->nt; kk++) {
                pt2 = &mesh->tria[kk];
                if ( !MG_EOK(pt2) )  continue;
                for (j=0; j<3; j++)
                    if ( pt2->v[j] == ip ) {
                        len++;
                        break;
                    }
            }
            if ( len != lon ) {
              if ( !mmgErr7 ) {
                mmgErr7 = 1;
                fprintf(stderr,"\n  ## Error: %s: 7. at least 1 incorrect"
                        " ball (%d: %d %d).\n",__func__,MMGS_indPt(mesh,ip),lon,len);
                ppt->tag |= MG_CRN + MG_REQ;
              }
              return 0;
            }
        }
    }

    return 1;
}

/**
 *
 * Check size prescriptions at point k.
 *
 * \warning not used
 */
int chkeigen(MMG5_pMesh mesh,MMG5_pSol met,int k,double lambda[3]) {
    MMG5_pPoint    p0;
    MMG5_pxPoint   go;
    int            ord;
    double         *m,*n,r[3][3],mr[6],mtan[3],vp[2][2];
    static char    mmgWarn=0;

    p0 = &mesh->point[k];
    assert( MG_VOK(p0) );

    m = &met->m[6*k];

    if ( MS_SIN(p0->tag) ) {
        lambda[0] = m[0];
        lambda[1] = m[3];
        lambda[2] = m[5];
    }
    else if ( p0->tag & MG_GEO ) {
        lambda[0] = m[0];
        lambda[1] = m[1];
        lambda[2] = m[2];
    }
    else {
        if ( p0->tag & MG_REF ) {
            go = &mesh->xpoint[p0->xp];
            n = &go->n1[0];
        }
        else
            n = &p0->n[0];

        if ( !MMG5_rotmatrix(n,r) )  return 0;
        MMG5_rmtr(r,m,mr);
        mtan[0] = mr[0];
        mtan[1] = mr[1];
        mtan[2] = mr[3];
        ord = MMG5_eigensym(mtan,lambda,vp);

        if ( !ord ) {
          if ( !mmgWarn ) {
            fprintf(stderr,"\n  ## Error: %s: at least 1 wrong matrix.\n",
                    __func__);
            mmgWarn = 1;
          }
          return 0;
        }
    }

    return 1;
}

/**
 * \return 0 if fail , 1 if success
 *
 * Check metric consistency.
 *
 * \warning not used.
 */
int chkmet(MMG5_pMesh mesh,MMG5_pSol met) {
    MMG5_pPoint    p0;
    MMG5_pxPoint   go;
    double         *m,*n,isqhmin,isqhmax,r[3][3],mr[6],mtan[3];
    double         vp[2][2],lambda[2];
    int            k;
    char           i;
    static char    mmgWarn0=0,mmgWarn1=0,mmgWarn2=0,mmgWarn3=0;

    isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
    isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

    /* First test : check wether metric at singular point has the suitable form */
    for(k=1; k<=mesh->np; k++) {
        p0 = &mesh->point[k];
        if ( !MG_VOK(p0) ) continue;

        if( MS_SIN(p0->tag) ) {
            m = &met->m[6*k];
            if( m[1] != 0.0 || m[2] != 0.0 || m[4] != 0.0 ){
              if ( !mmgWarn0 ) {
                mmgWarn0=1;
                fprintf(stderr,"\n   ## Error: %s: at least 1 wrong definition"
                        " of singular metric point (%d: "
                        "met %f %f %f %f %f %f).\n",__func__,
                        MMGS_indPt(mesh,k),m[0],m[1],m[2],m[3],m[4],m[5]);
              }
              return 0;
            }
            else if ( m[0] != m[3] || m[0] != m[5] || m[3] != m[5] )  {
              if ( !mmgWarn1 ) {
                mmgWarn1=1;
                fprintf(stderr,"\n  ## Error: %s: at least 1 wrong definition of"
                        " singular metric point (%d: "
                        "met %f %f %f %f %f %f).\n",__func__,
                        MMGS_indPt(mesh,k),m[0],m[1],m[2],m[3],m[4],m[5]);
              }
              return 0;
            }
        }
        else if ( p0->tag & MG_GEO ) {
            m = &met->m[6*k];
            for (i=0; i<3; i++) {
                if ( m[i] > isqhmin + 1.e-6 || m[i] < isqhmax - 1.e-6 ){
                  if ( !mmgWarn2 ) {
                    mmgWarn2=1;
                    fprintf(stderr,"\n  ## Error: %s: in definition of metric"
                            " at ridge point (%d: "
                            "met %f %f %f %f %f %f). \n",__func__,
                            MMGS_indPt(mesh,k),m[0],m[1],m[2],m[3],m[4],m[5]);
                    return 0;
                  }
                }
            }
        }
        else {
            m = &met->m[6*k];
            if ( MG_EDG(p0->tag) ) {
                go = &mesh->xpoint[p0->xp];
                n = &go->n1[0];
            }
            else{
                n = &p0->n[0];
            }

            /* Recovery of the eigenvalues of m */
            if ( !MMG5_rotmatrix(n,r) )  return 0;
            MMG5_rmtr(r,m,mr);
            mtan[0] = mr[0];
            mtan[1] = mr[1];
            mtan[2] = mr[3];
            MMG5_eigensym(mtan,lambda,vp);

            if ( lambda[0] > isqhmin + 1.e-6 || lambda[0] < isqhmax - 1.e-6 ){
              if ( !mmgWarn3 ) {
                mmgWarn3 = 1;
                fprintf(stderr,"\n  ## Error: %s: in definition of metric at"
                        " regular point (%d: "
                        "met %f %f %f %f %f %f).\n",__func__,MMGS_indPt(mesh,k),
                        m[0],m[1],m[2],m[3],m[4],m[5]);
                fprintf(stderr,"eigenvalue : %f \n",lambda[0]);
              }
              return 0;
            }

            if ( lambda[1] > isqhmin + 1.e-6 || lambda[1] < isqhmax - 1.e-6 ){
              if ( !mmgWarn3 ) {
                mmgWarn3 = 1;
                fprintf(stderr,"\n  ## Error: %s: in definition of metric at"
                        " regular point (%d: "
                        "met %f %f %f %f %f %f).\n",__func__,MMGS_indPt(mesh,k),
                        m[0],m[1],m[2],m[3],m[4],m[5]);
                fprintf(stderr,"eigenvalue : %f \n",lambda[1]);
              }
              return 0;
            }
        }
    }

    printf(" *** admissible metric.\n");

    return 1;
}

/**
 * \param mesh pointer toward the mesh
 *
 * \return 1 if success, 0 if fail.
 *
 * Check normal vectors consistency
 *
 * \warning unused
 */
int chknor(MMG5_pMesh mesh) {
    MMG5_pTria    pt;
    MMG5_pPoint   p0;
    MMG5_pxPoint    go;
    double   dd,ps,*n,nt[3];
    int      k;
    char     i;
    static char mmgWarn0=0, mmgWarn1=0;

    /* First test : check that all normal vectors at points are non 0 */
    for (k=1; k<=mesh->np; k++) {
        p0 = &mesh->point[k];
        if ( !MG_VOK(p0) ) continue;
        if ( MS_SIN(p0->tag) ) continue;
        if ( !(p0->tag & MG_GEO) ) continue;

        assert( p0->xp );
        go = &mesh->xpoint[p0->xp];
        n = &go->n1[0];

        dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if ( dd < 0.9 ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 non unitary normal"
                    " (point: %d normal n1 = %f %f %f). exit program\n",
                    __func__,MMGS_indPt(mesh,k),n[0],n[1],n[2]);
          }
          return 0;
        }

        n = &go->n2[0];

        dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if ( dd < 0.9 ) {
          if ( !mmgWarn0 ) {
            mmgWarn0 = 1;
            fprintf(stderr,"\n  ## Error: %s: at least 1 non unitary normal"
                    " (point: %d normal n2 = %f %f %f). exit program\n",
                    __func__,MMGS_indPt(mesh,k),n[0],n[1],n[2]);
          }
          return 0;
        }
    }

    /* Second test : check that all normal vectors at points are consistently oriented with
       respect to the underlying triangles */
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MG_EOK(pt) ) continue;

        MMG5_nortri(mesh,pt,nt);
        for (i=0; i<3; i++) {
            p0 = &mesh->point[pt->v[i]];
            if ( MS_SIN(p0->tag) ) continue;
            else if ( MG_EDG(p0->tag) ) {
                assert ( p0->xp );
                go = &mesh->xpoint[p0->xp];
                if ( p0->tag & MG_GEO ) {
                    n = &go->n1[0];
                    ps = n[0]*nt[0] + n[1]*nt[1] + n[2]*nt[2];
                    if ( ps < -0.99 ) {
                      if ( !mmgWarn1 ) {
                        mmgWarn1 = 1;
                        fprintf(stderr,"\n  ## Error: %s: at least 1"
                                " inconsistant normal (point %d in triangle %d):"
                                " ps = %f \n",__func__,MMGS_indPt(mesh,pt->v[i]),
                                MMGS_indElt(mesh,k),ps);
                      }
                      return 0;
                    }

                    n = &go->n2[0];
                    ps = n[0]*nt[0] + n[1]*nt[1] + n[2]*nt[2];
                    if ( ps < -0.99 ) {
                      if ( !mmgWarn1 ) {
                        mmgWarn1 = 1;
                        fprintf(stderr,"\n  ## Error: %s: at least 1"
                                " inconsistant normal (point %d in triangle %d):"
                                " ps = %f \n",__func__,MMGS_indPt(mesh,pt->v[i]),
                                MMGS_indElt(mesh,k),ps);
                      }
                      return 0;
                    }
                }
                else {
                    n = &go->n1[0];
                    ps = n[0]*nt[0] + n[1]*nt[1] + n[2]*nt[2];
                    if ( ps < -0.99 ) {
                      if ( !mmgWarn1 ) {
                        mmgWarn1 = 1;
                        fprintf(stderr,"\n  ## Error: %s: at least 1"
                                " inconsistant normal (point %d in triangle %d):"
                                " ps = %f \n",__func__,MMGS_indPt(mesh,pt->v[i]),
                                MMGS_indElt(mesh,k),ps);
                      }
                      return 0;
                    }
                }
            }
            else {
                n = &p0->n[0];
                ps = n[0]*nt[0] + n[1]*nt[1] + n[2]*nt[2];
                if ( ps < -0.99 ) {
                  if ( !mmgWarn1 ) {
                    mmgWarn1 = 1;
                    fprintf(stderr,"\n  ## Error: %s: at least 1"
                            " inconsistant normal (point %d in triangle %d):"
                            " ps = %f \n",__func__,MMGS_indPt(mesh,pt->v[i]),
                            MMGS_indElt(mesh,k),ps);
                  }
                  return 0;
                }
            }
        }
    }

    printf(" *** admissible normals.\n");

    return 1;
}
