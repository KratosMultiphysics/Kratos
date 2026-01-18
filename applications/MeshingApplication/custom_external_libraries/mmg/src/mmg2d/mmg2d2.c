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
 * \file mmg2d/mmg2d2.c
 * \brief Mesh generation functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "libmmg2d_private.h"
#include "mmgexterns_private.h"

/**
 * \param mesh pointer to the mesh structure.
 * \return 1 if success, 0 if fail.
 *
 * Remove the bounding box triangles.
 *
 */
int MMG2D_removeBBtriangles(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_int        ip1,ip2,ip3,ip4,k,iadr,*adja,iadr2,*adja2,iel,nd;
  int8_t          i,ii;
  static int8_t   mmgWarn0=0;

  /* Bounding Box vertices */
  ip1 = mesh->np-3;
  ip2 = mesh->np-2;
  ip3 = mesh->np-1;
  ip4 = mesh->np;

  nd = 0;
  for(k=1; k<=mesh->nt; k++) {
    pt  = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    if ( pt->base < 0 ) {
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for(i=0; i<3; i++) {
        if ( !adja[i] ) continue;
        iel = adja[i] / 3;
        ii = adja[i] % 3;
        iadr2 = 3*(iel-1) + 1;
        adja2 = &mesh->adja[iadr2];
        adja2 [ii] = 0;
      }
      MMG2D_delElt(mesh,k);
      continue;
    }
    else if ( !pt->base ) {
      if ( !mmgWarn0 ) {
        mmgWarn0 = 1;
        fprintf(stderr,"\n  ## Warning: %s: at least 1 undetermined"
                " triangle.\n",__func__);
      }
      nd++;
    }
  }

  if ( !nd ) {
    MMG2D_delPt(mesh,ip1);
    MMG2D_delPt(mesh,ip2);
    MMG2D_delPt(mesh,ip3);
    MMG2D_delPt(mesh,ip4);
  }
  else {
    fprintf(stderr,"\n  ## Error: %s: procedure failed :"
            " %" MMG5_PRId " indetermined triangles.\n",__func__,nd);
    return 0;
  }
  return 1;
}

/* Set tag to triangles in the case where there are no constrained edge
   in the supplied mesh: in = base ; out = -base ; undetermined = 0*/
int MMG2D_settagtriangles(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria        pt;
  int               iter,maxiter;
  MMG5_int          base,nd,k,ip1,ip2,ip3,ip4;

  /*BB vertex*/
  ip1=(mesh->np-3);
  ip2=(mesh->np-2);
  ip3=(mesh->np-1);
  ip4=(mesh->np);

  base = ++mesh->base;
  iter    = 0;
  maxiter = 3;
  do {
    nd = 0;
    for(k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;
      if ( !MMG2D_findtrianglestate(mesh,k,ip1,ip2,ip3,ip4,base) ) nd++ ;
    }

    if(mesh->info.ddebug) {
      printf(" ** how many undetermined triangles ? %" MMG5_PRId "\n",nd);
    }
  }
  while (nd && ++iter<maxiter);

  return 1;
}

/* Find out whether triangle pt is inside or outside (i.e. contains bb points or not) */
/* Return <0 value if triangle outside ; > 0 if triangle inside */
MMG5_int MMG2D_findtrianglestate(MMG5_pMesh mesh,MMG5_int k,MMG5_int ip1,MMG5_int ip2,MMG5_int ip3,MMG5_int ip4,MMG5_int base) {
  MMG5_pTria       pt;
  int              nb;
  int8_t           i;

  pt = &mesh->tria[k];

  /* Count how many vertices of pt are vertices of the boundary box */
  nb = 0;
  for(i=0; i<3; i++)
    if ( pt->v[i] == ip1 || pt->v[i] == ip2 || pt->v[i] == ip3 || pt->v[i] == ip4 ) nb++;

  /* Triangle to be deleted */
  if ( nb ) {
    pt->base = -base;
    pt->ref = 3;
    return -base;
  }
  else {
    pt->base = base;
    return base;
  }
}

/**
 * \param mesh pointer to mesh
 * \param k point index
 * \return index of one elt containing k or 0 (if no elt is found)
 *
 * Return the index of one triangle containing k (performing exhausting search
 * if needed).
 */
static inline
MMG5_int MMG2D_findTria_exhaust(MMG5_pMesh mesh,MMG5_int k) {
  MMG5_pPoint ppt = &mesh->point[k];
  static int8_t mmgWarn0=0;

  /* Find the triangle lel of the mesh containing ppt */
  MMG5_int lel = MMG2D_findTria(mesh,k);

  /* Exhaustive search if not found */
  if ( !lel ) {
    if ( mesh->info.ddebug ) {
      printf(" ** exhaustive search of point location.\n");
    }

    int kk;
    for(kk=1; kk<=mesh->nt; kk++) {
      lel = MMG2D_isInTriangle(mesh,kk,&ppt->c[0]);
      if ( lel ) break;
    }

    if ( kk>mesh->nt ) {
      if ( !mmgWarn0 ) {
        mmgWarn0 = 1;
        fprintf(stderr,"\n  ## Error: %s: unable to find triangle"
                " for at least vertex %" MMG5_PRId ".\n",__func__,k);
      }
      return 0;
    }
  }
  return lel;
}


/**
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the solution structure
 * \return  0 if fail.
 *
 * Insertion of the list of points inside the mesh
 * (Vertices mesh->np - 3, 2, 1, 0 are the vertices of the BB and have already been inserted)
 *
 */
int MMG2D_insertpointdelone(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pPoint   ppt;
  int           lon;
  MMG5_int      k,kk,ns,nus,nu,nud;
  int           iter,maxiter;
  static int8_t mmgWarn0=0,mmgWarn1=0,mmgWarn2=0;
  MMG5_int      list[MMG5_TRIA_LMAX];
  const int flag=-10;

  for(k=1; k<=mesh->np-4; k++) {
    ppt = &mesh->point[k];
    ppt->flag	= flag;
  }
  iter = 0;
  maxiter = 10;

	do {
    ns = nus = 0;
    nu = nud = 0;
    mmgWarn1 = mmgWarn2 = 0;
    for(k=1; k<=mesh->np-4; k++) {
      ppt = &mesh->point[k];
		  if(ppt->flag != flag) continue;
			nus++;

      list[0] = MMG2D_findTria_exhaust(mesh,k);
      if ( !list[0] ) {
        return 0;
      }

      /* Create the cavity of point k starting from list[0] */
      lon = MMG2D_cavity(mesh,sol,k,list);

      if ( lon < 1 ) {
        nu++;
        if ( !mmgWarn1 ) {
          mmgWarn1 = 1;
          if ( mesh->info.imprim > 6 || mesh->info.ddebug )
	          fprintf(stderr,"\n  ## Warning: %s: unable to insert "
		          "at least 1 vertex. (%" MMG5_PRId ")\n",__func__,k);
        }
        continue;
      } else {
				if(!MMG2D_delone(mesh,sol,k,list,lon)) {
			    if ( abs(mesh->info.imprim) > 4) {
            nud++;
            if ( !mmgWarn2 ) {
              mmgWarn2 = 1;
              if(mesh->info.imprim > 6 || mesh->info.ddebug)
	             	fprintf(stderr,"\n  ## Warning: %s: unable to"
			            " insert at least 1 point with Delaunay (%" MMG5_PRId ")\n",__func__,k);
            }
          }
        } else {
          ppt->flag = 0;
          ns++;
        }
      }
    }

    if ( abs(mesh->info.imprim) > 4)
      fprintf(stdout,"     %8" MMG5_PRId " vertex inserted %8" MMG5_PRId " not inserted\n",ns,nu+nud);
    if ( mesh->info.imprim >6 || mesh->info.ddebug )
      fprintf(stdout,"     unable to insert %8" MMG5_PRId " vertex : cavity %8" MMG5_PRId " -- delaunay %8" MMG5_PRId " \n",nu+nud,nu,nud);
  } while (ns && ++iter<maxiter);

	if(MMG5_abs(nus-ns)) {
    if ( mesh->info.imprim > 6 || mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: unable to"
              " insert %8" MMG5_PRId " point with Delaunay \n",__func__,MMG5_abs(nus-ns));
      fprintf(stdout,"     try to insert with splitbar\n");
    }
    mmgWarn2 = 0;
    nus = ns = 0;
    /*try to insert using splitbar*/
    for(k=1; k<=mesh->np-4; k++) {
      ppt = &mesh->point[k];
		  if(ppt->flag != flag) continue;
			nus++;

      list[0] = MMG2D_findTria_exhaust(mesh,k);
      if ( !list[0] ) {
        return 0;
      }

      if(!MMG2D_splitbar(mesh,list[0],k)) {
        if ( !mmgWarn2 ) {
          mmgWarn2 = 1;
          if ( mesh->info.imprim >6 || mesh->info.ddebug )
            fprintf(stderr,"\n  ## Warning: %s: unable to"
                    " insert at least 1 point with splitbar (%" MMG5_PRId ")\n",__func__,k);
        }
      } else {
        ns++;
      }
    }
    if ( MMG5_abs(nus-ns) ) {
      fprintf(stderr,"  ## Warning: %s: %" MMG5_PRId " point(s) not "
            "inserted. Check your output mesh\n",__func__,MMG5_abs(nus-ns));
      return 0;
    }
  }
	return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 *
 * \return 0 if fail, 1 if success.
 *
 * Put different references on different subdomains
 *
 */
int MMG2D_markSD(MMG5_pMesh mesh) {
  MMG5_pTria   pt,pt1;
  MMG5_pEdge   ped;
  MMG5_pPoint  ppt;
  MMG5_int     k,l,iadr,*adja,ped0,ped1,ipil,ncurc,nref;
  MMG5_int     kinit,nt,nsd,ip1,ip2,ip3,ip4,ned,iel;
  int          voy;
  int8_t       i,i1,i2;
  MMG5_int     *list;

  /* Reset flag field for triangles */
  for(k=1 ; k<=mesh->nt ; k++)
    mesh->tria[k].flag = mesh->mark;

  MMG5_SAFE_CALLOC(list,mesh->nt,MMG5_int,return 0);
  kinit = 0;
  nref  = 0;
  ip1   =  mesh->np;

  /* Catch first triangle with vertex ip1 */
  for(k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    pt->flag = mesh->mark;
    pt->ref  = 0;
    if ( (!kinit) && ( pt->v[0]==ip1 || pt->v[1]==ip1 || pt->v[2]==ip1) ) kinit = k;
  }

  /* Travel mesh by adjacencies to set references on triangles as long as no boundary is met */
  do {
    nref++;
    list[0] = kinit;
    ipil = 0;
    ncurc = 0;
    do {
      k = list[ipil];
      pt = &mesh->tria[k];
      pt->ref = nref;
      iadr = 3*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for(i=0; i<3; i++) {
        iel = adja[i] / 3;
        pt1 = &mesh->tria[iel];

        if( !iel || pt1->ref == nref ) continue;

        i1 = MMG5_inxt2[i];
        i2 = MMG5_iprv2[i];
        ped0 = pt->v[i1];
        ped1 = pt->v[i2];

        /* WARNING: exhaustive search among edges, to be optimized with a hashing structure */
        for(l=1; l<=mesh->na; l++) {
          ped = &mesh->edge[l];
          if( ( ped->a == ped0 && ped->b == ped1 ) || ( ped->b == ped0 && ped->a == ped1 ) ) break;
        }
        if ( l <= mesh->na ) continue;

        pt1->ref = nref;
        list[++ncurc] = iel;
      }
      ++ipil ;
    }
    while ( ipil <= ncurc );

    kinit = 0;
    for(k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;
      pt->flag = mesh->mark;
      if ( !kinit && !(mesh->tria[k].ref) ) kinit = k;
    }
  }
  while ( kinit );

  if ( mesh->info.imprim > 0  ) {
    /* nref - 1 subdomains because Bounding Box triangles have been counted */
    fprintf(stdout,"     %8" MMG5_PRId " sub-domains\n",nref-1);
  }

  MMG5_SAFE_FREE(list);

  /* Remove BB triangles and vertices */
  /*BB vertex*/
  ip1=(mesh->np-3);
  ip2=(mesh->np-2);
  ip3=(mesh->np-1);
  ip4=(mesh->np);

  /* Case when there are inner and outer triangles */
  if ( nref != 1 ) {
    nt = mesh->nt;
    for(k=1; k<=nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      if ( pt->ref != 1 ) continue;
      /*update adjacencies*/
      iadr = 3*(k-1)+1;
      adja = &mesh->adja[iadr];
      for(i=0; i<3; i++) {
        if ( !adja[i] ) continue;
        iel = adja[i] / 3;
        voy = adja[i] % 3;
        (&mesh->adja[3*(iel-1)+1])[voy] = 0;
      }
      MMG2D_delElt(mesh,k);
    }
  }
  /* Remove only the triangles containing one of the BB vertex*/
  else {
    nt = mesh->nt;
    for(k=1 ; k<=nt ; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;

      if( !(pt->v[0]==ip1 || pt->v[1]==ip1 || pt->v[2]==ip1 ||
           pt->v[0]==ip2 || pt->v[1]==ip2 || pt->v[2]==ip2 ||
           pt->v[0]==ip3 || pt->v[1]==ip3 || pt->v[2]==ip3 ||
           pt->v[0]==ip4 || pt->v[1]==ip4 || pt->v[2]==ip4 ) ) continue;

      /*update adjacencies*/
      iadr = 3*(k-1)+1;
      adja = &mesh->adja[iadr];
      for(i=0 ; i<3 ; i++) {
        if(!adja[i]) continue;
        iel = adja[i]/3;
        voy = adja[i]%3;
        (&mesh->adja[3*(iel-1)+1])[voy] = 0;
      }
      MMG2D_delElt(mesh,k);
    }
  }

  MMG2D_delPt(mesh,ip1);
  MMG2D_delPt(mesh,ip2);
  MMG2D_delPt(mesh,ip3);
  MMG2D_delPt(mesh,ip4);

  if(mesh->info.nsd) {
    nsd = mesh->info.nsd;
    nt = mesh->nt;
    for(k=1 ; k<=nt ; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;
      if ( pt->ref == nsd ) continue;
      MMG2D_delElt(mesh,k);
    }
  }

  /* Remove vertex*/
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for (i=0; i<3; i++) {
      ppt = &mesh->point[ pt->v[i] ];
      ppt->tag &= ~MG_NUL;
    }
  }
  /* Remove edge*/
  ned = mesh->na;
  for (k=1; k<=ned; k++) {
    ped = &mesh->edge[k];
    if ( !ped->a )  continue;
    ppt = &mesh->point[ ped->a ];
    if ( !MG_VOK(ppt) ) {
      MMG5_delEdge(mesh,k);
      continue;
    }
    ppt = &mesh->point[ ped->b ];
    if ( !MG_VOK(ppt) ) {
      MMG5_delEdge(mesh,k);
      continue;
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the sol structure.
 * \return 0 if fail, 1 if success.
 *
 * Mesh triangulation.
 *
 **/
int MMG2D_mmg2d2(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria pt;
  double     c[2];
  MMG5_int   ip1,ip2,ip3,ip4,jel,kel,iadr,*adja;

  mesh->base = 0;
  assert ( !mesh->nt );

  /* Create the 4 vertices of the bounding box */
  /* Bottom left corner */
  c[0] = -0.5; //mesh->info.min[0] - 1.;
  c[1] = -0.5; // mesh->info.min[1] - 1.;
  ip1 = MMG2D_newPt(mesh,c,0);
  if ( !ip1 ) {
    /* reallocation of point table */
    MMG2D_POINT_REALLOC(mesh,sol,ip1,mesh->gap,
                         fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                 " a new point.\n",__func__);
                         MMG5_INCREASE_MEM_MESSAGE();return 0;,
                         c,0);
  }

  /* Top left corner */
  c[0] = -0.5; //mesh->info.min[0] - 1.;
  c[1] =  MMG2D_PRECI / mesh->info.delta *(mesh->info.max[1]-mesh->info.min[1]) + 0.5;//mesh->info.max[1] + 1.;
  ip2 = MMG2D_newPt(mesh,c,0);
  if ( !ip2 ) {
    /* reallocation of point table */
    MMG2D_POINT_REALLOC(mesh,sol,ip2,mesh->gap,
                         fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                 " a new point.\n",__func__);
                         MMG5_INCREASE_MEM_MESSAGE();return 0;,
                         c,0);
  }

  /* Bottom right corner */
  c[0] =  MMG2D_PRECI / mesh->info.delta *(mesh->info.max[0]-mesh->info.min[0]) + 0.5;//mesh->info.max[0] + 1.;
  c[1] = -0.5;//mesh->info.min[1] - 1.;
  ip3 = MMG2D_newPt(mesh,c,0);
  if ( !ip3 ) {
    /* reallocation of point table */
    MMG2D_POINT_REALLOC(mesh,sol,ip3,mesh->gap,
                         fprintf(stderr,"\n  ## Error: %s: unable to allocate a "
                                 " new point.\n",__func__);
                         MMG5_INCREASE_MEM_MESSAGE(); return 0;,
                         c,0);
  }

  /* Top right corner */
  c[0] =  MMG2D_PRECI / mesh->info.delta *(mesh->info.max[0]-mesh->info.min[0]) + 0.5;//mesh->info.max[0] + 1.;
  c[1] = MMG2D_PRECI / mesh->info.delta *(mesh->info.max[1]-mesh->info.min[1]) + 0.5;//mesh->info.max[1] + 1.;
  ip4 = MMG2D_newPt(mesh,c,0);
  if ( !ip4 ) {
    /* reallocation of point table */
    MMG2D_POINT_REALLOC(mesh,sol,ip4,mesh->gap,
                         fprintf(stderr,"\n  ## Error: %s: unable to allocate a"
                                 " new point.\n",__func__);
                         MMG5_INCREASE_MEM_MESSAGE(); return 0;,
                         c,0);
  }

  assert ( ip1 == mesh->np-3 );
  assert ( ip2 == mesh->np-2 );
  assert ( ip3 == mesh->np-1 );
  assert ( ip4 == mesh->np );

  /* Create the first two triangles in the mesh and the adjacency relations */
  jel  = MMG2D_newElt(mesh);
  if ( !jel ) {
    MMG2D_TRIA_REALLOC(mesh,jel,mesh->gap,
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate a"
                               " new element.\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");return 0);
  }
  pt       = &mesh->tria[jel];
  pt->v[0] = ip1;
  pt->v[1] = ip4;
  pt->v[2] = ip2;
  pt->base = mesh->base;

  kel  = MMG2D_newElt(mesh);
  if ( !kel ) {
    MMG2D_TRIA_REALLOC(mesh,kel,mesh->gap,
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                               " a new element.\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();
                       printf("  Exit program.\n");return 0);
  }
  pt   = &mesh->tria[kel];
  pt->v[0] = ip1;
  pt->v[1] = ip3;
  pt->v[2] = ip4;
  pt->base = mesh->base;

  iadr = 3*(jel-1) + 1;
  adja = &mesh->adja[iadr];
  adja[2] = 3*kel + 1;

  iadr = 3*(kel-1) + 1;
  adja = &mesh->adja[iadr];
  adja[1] = 3*jel + 2;

  /* Insertion of vertices in the mesh */
  if ( !MMG2D_insertpointdelone(mesh,sol) ) return 0;

  if ( mesh->info.imprim > 0 )
    fprintf(stdout,"     Insertion succeed\n");

  /* Enforcement of the boundary edges */
  if ( !MMG2D_bdryenforcement(mesh,sol) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to enforce the boundaries.\n",
      __func__);
    return 0;
  }

  if(mesh->info.ddebug)
    if ( !MMG5_chkmsh(mesh,1,0) ) return 0;

  /* Mark SubDomains and remove the bounding box triangles */
  if ( mesh->na ) {
   if ( ! MMG2D_markSD(mesh) ) return 0;
  }
  else {
    /* Tag triangles : in = base ; out = -base ; Undetermined = 0*/
    if ( !MMG2D_settagtriangles(mesh,sol) ) return 0;
    if ( !MMG2D_removeBBtriangles(mesh) ) return 0;
  }

  return 1;
}
