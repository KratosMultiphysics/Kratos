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
int _MMG5_mmgsChkmsh(MMG5_pMesh mesh,int severe,int base) {
    MMG5_pPoint         ppt;
    MMG5_pTria          pt1,pt2;
    int                 adj,adj1,k,kk,l,nk,i,j,ip,lon,len;
    int                 *adja,*adjb,list[_MMG5_LMAX+2];
    char                voy,voy1,i1,i2,j1,j2;

    for (k=1; k<=mesh->nt; k++) {
        pt1 = &mesh->tria[k];
        if ( !MG_EOK(pt1) )  continue;
        adja = &mesh->adja[3*(k-1)+1];

        for (i=0; i<3; i++) {
            if ( !adja[i] )  continue;
            i1  = _MMG5_inxt2[i];
            i2  = _MMG5_iprv2[i];
            adj = adja[i] / 3;
            voy = adja[i] % 3;
            if ( !adj && !(pt1->tag[i] & MG_GEO) ) {
                fprintf(stderr,"  0. Missing edge tag %d %d\n",k,adj);
                fprintf(stderr,"triangle %d: %d %d %d \n",k,pt1->v[0],pt1->v[1],pt1->v[2]);
                fprintf(stderr,"tag (%d): %d %d %d \n",k,pt1->tag[0],pt1->tag[1],pt1->tag[2]);
                return(0);
            }
            if ( adj == k ) {
                fprintf(stderr,"  1. Wrong adjacency %d %d\n",k,adj);
                fprintf(stderr,"triangle %d: %d %d %d \n",k,pt1->v[0],pt1->v[1],pt1->v[2]);
                fprintf(stderr,"adj (%d): %d %d %d \n",k,adja[0]/3,adja[1]/3,adja[2]/3);
                return(0);
            }
            pt2 = &mesh->tria[adj];
            if ( !MG_EOK(pt2) ) {
                fprintf(stderr,"  4. Invalid adjacent %d %d\n",adj,k);
                fprintf(stderr,"vertices of k   %d: %d %d %d\n",
                       k,pt1->v[0],pt1->v[1],pt1->v[2]);
                fprintf(stderr,"vertices of adj %d: %d %d %d \n",
                       adj,pt2->v[0],pt2->v[1],pt2->v[2]);
                return(0);
            }
            if ( (pt1->tag[i] != pt2->tag[voy]) || (pt1->edg[i] != pt2->edg[voy] ) ) {
                fprintf(stderr,"  3. Wrong tag/ref %d %d  %d - %d\n",
                        k,adj,pt1->tag[i],pt2->tag[voy]);
                return(0);
            }
            adjb = &mesh->adja[3*(adj-1)+1];
            adj1 = adjb[voy] / 3;
            voy1 = adjb[voy] % 3;
            if ( adj1 != k || voy1 != i ) {
                fprintf(stderr,"  2. Wrong adjacency %d %d\n",k,adj1);
                fprintf(stderr,"vertices of %d: %d %d %d \n",k,
                       pt1->v[0],pt1->v[1],pt1->v[2]);
                fprintf(stderr,"vertices of adj %d: %d %d %d \n",adj,
                       pt2->v[0],pt2->v[1],pt2->v[2]);
                fprintf(stderr,"adj(%d): %d %d %d\n",k,adja[0]/3,adja[1]/3,adja[2]/3);
                fprintf(stderr,"adj(%d): %d %d %d\n",adj,adjb[0]/3,adjb[1]/3,adjb[2]/3);
                return(0);
            }
            if ( !MS_SIN(pt1->tag[i]) ) {
                j1 = _MMG5_inxt2[voy];
                j2 = _MMG5_iprv2[voy];
                if ( pt2->v[j2] != pt1->v[i1] || pt2->v[j1] != pt1->v[i2] ) {
                    fprintf(stderr,"  8. Wrong orientation %d %d\n",k,adj);
                    return(0);
                }
            }
        }
    }

    if ( !severe )  return(1);

    for (k=1; k<=mesh->nt; k++) {
        pt1 = &mesh->tria[k];
        if ( !MG_EOK(pt1) )  continue;

        adja = &mesh->adja[3*(k-1)+1];
        for (i=0; i<3; i++) {
            if ( !adja[i] )  continue;

            ip  = pt1->v[i];
            ppt = &mesh->point[ip];
            if ( !MG_VOK(ppt) ) {
                fprintf(stderr,"  6. Unused vertex %d  %d\n",k,ip);
                fprintf(stderr,"%d %d %d\n",pt1->v[0],pt1->v[1],pt1->v[2]);
                return(0);
            }
            else if ( MS_SIN(ppt->tag) )  continue;

            lon = boulet(mesh,k,i,list);
            if ( lon < 1 )  continue;
            for (l=0; l<lon; l++) {
                kk  = list[l] / 3;
                nk  = list[l] % 3;
                pt2 = &mesh->tria[kk];
                if ( pt2->v[nk] != ip ) {
                    fprintf(stderr,"  5. Wrong ball %d, %d\n",ip,pt2->v[nk]);
                    return(0);
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
                fprintf(stderr,"  7. Incorrect ball %d: %d %d\n",ip,lon,len);
                ppt->tag |= MG_CRN + MG_REQ;
            }
        }
    }

    return(1);
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

        if ( !_MMG5_rotmatrix(n,r) )  return(0);
        _MMG5_rmtr(r,m,mr);
        mtan[0] = mr[0];
        mtan[1] = mr[1];
        mtan[2] = mr[3];
        ord = _MMG5_eigensym(mtan,lambda,vp);

        if ( !ord ) {
            fprintf(stderr,"wrong matrix \n");
            exit(EXIT_FAILURE);
        }
    }

    return(1);
}

/**
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

    isqhmin = 1.0 / (mesh->info.hmin*mesh->info.hmin);
    isqhmax = 1.0 / (mesh->info.hmax*mesh->info.hmax);

    /* First test : check wether metric at singular point has the suitable form */
    for(k=1; k<=mesh->np; k++) {
        p0 = &mesh->point[k];
        if ( !MG_VOK(p0) ) continue;

        if( MS_SIN(p0->tag) ) {
            m = &met->m[6*k];
            if( m[1] != 0.0 || m[2] != 0.0 || m[4] != 0.0 ){
                fprintf(stderr,"   ### Error in definition of singular metric point %d,\
                     met %f %f %f %f %f %f  \n",k,m[0],m[1],m[2],m[3],m[4],m[5]);
                exit(EXIT_FAILURE);
            }
            else if ( m[0] != m[3] || m[0] != m[5] || m[3] != m[5] )  {
                fprintf(stderr,"   ### Error in definition of singular metric point %d,\
               met %f %f %f %f %f %f  \n",k,m[0],m[1],m[2],m[3],m[4],m[5]);
                exit(EXIT_FAILURE);
            }
        }
        else if ( p0->tag & MG_GEO ) {
            m = &met->m[6*k];
            for (i=0; i<3; i++) {
                if ( m[i] > isqhmin + 1.e-6 || m[i] < isqhmax - 1.e-6 ){
                    fprintf(stderr,"   ### Error in definition of metric at ridge point %d,\
                 met %f %f %f %f %f %f  \n",k,m[0],m[1],m[2],m[3],m[4],m[5]);
                    exit(EXIT_FAILURE);
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
            if ( !_MMG5_rotmatrix(n,r) )  return(0);
            _MMG5_rmtr(r,m,mr);
            mtan[0] = mr[0];
            mtan[1] = mr[1];
            mtan[2] = mr[3];
            _MMG5_eigensym(mtan,lambda,vp);

            if ( lambda[0] > isqhmin + 1.e-6 || lambda[0] < isqhmax - 1.e-6 ){
                fprintf(stderr,"   ### Error in definition of metric at regular point %d,\
               met %f %f %f %f %f %f  \n",k,m[0],m[1],m[2],m[3],m[4],m[5]);
                fprintf(stderr,"eigenvalue : %f \n",lambda[0]);
                exit(EXIT_FAILURE);
            }

            if ( lambda[1] > isqhmin + 1.e-6 || lambda[1] < isqhmax - 1.e-6 ){
                fprintf(stderr,"   ### Error in definition of metric at regular point %d,\
               met %f %f %f %f %f %f  \n",k,m[0],m[1],m[2],m[3],m[4],m[5]);
                fprintf(stderr,"eigenvalue : %f \n",lambda[1]);
                exit(EXIT_FAILURE);
            }
        }
    }

    printf(" *** admissible metric.\n");

    return(1);
}

/* Check normal vectors consistency */
int chknor(MMG5_pMesh mesh) {
    MMG5_pTria    pt;
    MMG5_pPoint   p0;
    MMG5_pxPoint    go;
    double   dd,ps,*n,nt[3];
    int      k;
    char     i;

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
            fprintf(stderr,"point : %d normal n1 = %f %f %f : exit program\n",k,n[0],n[1],n[2]);
            exit(EXIT_FAILURE);
        }

        n = &go->n2[0];

        dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if ( dd < 0.9 ) {
            fprintf(stderr,"point : %d normal n2 = %f %f %f : exit program\n",k,n[0],n[1],n[2]);
            exit(EXIT_FAILURE);
        }
    }

    /* Second test : check that all normal vectors at points are consistently oriented with
       respect to the underlying triangles */
    for (k=1; k<=mesh->nt; k++) {
        pt = &mesh->tria[k];
        if ( !MG_EOK(pt) ) continue;

        _MMG5_nortri(mesh,pt,nt);
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
                        fprintf(stderr,"point %d in triangle %d : inconsistant normal : ps = %f \n",pt->v[i],k,ps);
                        exit(EXIT_FAILURE);
                    }

                    n = &go->n2[0];
                    ps = n[0]*nt[0] + n[1]*nt[1] + n[2]*nt[2];
                    if ( ps < -0.99 ) {
                        fprintf(stderr,"point %d in triangle %d : inconsistant normal : ps = %f \n",pt->v[i],k,ps);
                        exit(EXIT_FAILURE);
                    }

                }
                else {
                    n = &go->n1[0];
                    ps = n[0]*nt[0] + n[1]*nt[1] + n[2]*nt[2];
                    if ( ps < -0.99 ) {
                        fprintf(stderr,"point %d in triangle %d : inconsistant normal : ps = %f \n",pt->v[i],k,ps);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            else {
                n = &p0->n[0];
                ps = n[0]*nt[0] + n[1]*nt[1] + n[2]*nt[2];
                if ( ps < -0.99 ) {
                    fprintf(stderr,"point %d in triangle %d : inconsistant normal : ps = %f \n",pt->v[i],k,ps);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    printf(" *** admissible normals.\n");

    return(1);
}
