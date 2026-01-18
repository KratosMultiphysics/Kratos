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
 * \file mmg/mmg2.c
 * \author Algiane Froehly (Bx INP/Inria/UBordeaux)
 * \author Luca Cirrottola (Inria)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * Common functions for ls discretization.
 *
 */

#include "mmgcommon_private.h"
#include "mmgexterns_private.h"

/**
 * \param pim    multimaterials inverse data table.
 * \param ref    material reference.
 * \return the key of the material in the lookup/hash table.
 *
 * Compute key for the material in the hash table.
 */
static MMG5_int MMG5_InvMat_key(MMG5_pInvMat pim,int ref) {
  return (ref - pim->offset);
}

/**
 * \param k      index of the material.
 * \param attr   the attribute of the material (nosplit/split/plus/minus).
 * \return the code to be stored in the inverse data table entry.
 *
 * Compute the material entry in the inverse data table, storing the index of
 * the parent material and the attribute of the child material.
 */
static int MMG5_InvMat_code(int k,int attr) {
  return 4*(k+1)+attr;
}

/**
 * \param pim multimaterials inverse data table.
 * \param ref    material reference.
 * \return Index of the material.
 *
 * Get index of the parent material from lookup table.
 */
static int MMG5_InvMat_getIndex(MMG5_pInvMat pim,int ref) {
  int key = MMG5_InvMat_key(pim,ref);
  /* The parent index is stored as 4*(k+1) */
  return (pim->lookup[key] / 4 - 1);
}

/**
 * \param mesh   pointer to the mesh structure.
 * \param pim    multimaterials inverse data table.
 * \param ref    material reference.
 * \return the nosplit/split/plus/minus attribute of the material.
 *
 * Get attribute of the child material (nosplit/split/plus/minus) from lookup
 * table.
 */
static int MMG5_InvMat_getAttrib(MMG5_pInvMat pim,int ref) {
  int key = MMG5_InvMat_key(pim,ref);
  /* The nosplit/split/plus/minus attribute is stored as the rest of the
   * integer division. */
  return (pim->lookup[key] % 4);
}

/**
 * \param pim    multimaterials inverse data table.
 * \param key    material key in the inverse data table.
 * \return 0 if an entry for the material already exists, 1 otherwise.
 *
 * Check if a material reference already exists in the material lookup table.
 */
static int MMG5_InvMat_check(MMG5_pInvMat pim,int key) {
  return pim->lookup[key] ? 0 : 1;
}

/**
 * \param pim    multimaterials inverse data table.
 * \param ref    reference of the material..
 * \param k      index of the material in the direct table.
 * \return 0 if an entry for the material already exists, 1 otherwise.
 *
 * Raise an error if trying to overwrite a reference entry in the material
 * lookup table.
 */
static void MMG5_InvMat_error(MMG5_pInvMat pim,int ref,int k) {
  fprintf(stderr,"\n   ## Warning: Overwrite material reference %d"
    " (from LSReferences line %d) with another entry from LSReferences line %d."
    ,ref,MMG5_InvMat_getIndex(pim,ref)+1,k+1);
  fprintf(stderr,"\n               Check your LSReferences table: if possible,"
          " each material reference should be unique,\n"
          "                if not possible, you may"
          " encounter unexpected issues (wrong domain mapping or erroneous"
          " detection of non-manifold level-set)!\n");
}

/**
 * \param mesh   pointer to the mesh structure.
 * \param pim    multimaterials inverse data table.
 * \param k      index of the material in the input table.
 *
 * Set materials lookup table entry.
 */
static int MMG5_InvMat_set(MMG5_pMesh mesh,MMG5_pInvMat pim,int k) {
  MMG5_pMat pm;
  int       key;

  /* Get material */
  pm = &mesh->info.mat[k];

  /** Store the dosplit attribute of the parent material */
  key = MMG5_InvMat_key(pim,pm->ref);
  if( !MMG5_InvMat_check(pim,key) ) {
    MMG5_InvMat_error(pim,pm->ref,k);
  }
  pim->lookup[key] = MMG5_InvMat_code(k,pm->dospl);

  /** Store the child material sign with the parent material index (in the
   *  lookup table).
   *  1) 0 is a legit material index, so store the parent as 4*(k+1).
   *  2) If a child material has the same reference as the parent, this
   *     effectively overwrites the result of the previous instruction.
   *  3) No different child materials are allowed to have the same reference,
   *     and this must have already been checked. */
  if( pm->dospl ) {
    key = MMG5_InvMat_key(pim,pm->rin);
    if( !MMG5_InvMat_check(pim,key) ) {
      MMG5_InvMat_error(pim,pm->rin,k);
    }
    pim->lookup[key] = MMG5_InvMat_code(k,MG_MINUS);

    key = MMG5_InvMat_key(pim,pm->rex);
    if( !MMG5_InvMat_check(pim,key) ) {
      MMG5_InvMat_error(pim,pm->rex,k);
    }
    pim->lookup[key] = MMG5_InvMat_code(k,MG_PLUS);
  }

  return 1;
}

/**
 * \param mesh   pointer to the mesh structure.
 * \param pim    multimaterials inverse data table.
 * \param ref    material reference.
 * \param pref   pointer to the parent material reference.
 * \return 1 if found, 0 if not found.
 *
 * Get reference of the parent material in multimaterial mode.
 *.
 */
static int MMG5_InvMat_getParent(MMG5_pMesh mesh,MMG5_pInvMat pim,MMG5_int ref,MMG5_int *pref) {
  MMG5_pMat pm;
  int       k;

  /* The parent index */
  k = MMG5_InvMat_getIndex(pim,ref);

  /* Material not found in the table */
  if( k == -1 ) {
    fprintf(stderr,"\n  ## Warning: %s: material %" MMG5_PRId " not found in table.\n",
            __func__,ref);
    fprintf(stderr,"              Please ensure that you provide all mesh"
            " references in the material map\n"
            "              (that is, the whole list of"
            " surface materials in lssurf mode,\n"
            "              and the whole list of domain"
            " materials in ls mode).\n" );
    return 0;
  }

  /* Get the material in the lookup table and return the parent reference */
  pm = &mesh->info.mat[k];
  *pref = pm->ref;
  return 1;
}

/**
 * \param mesh pointer to the mesh
 * \param ref  final reference for which we are searching the initial one
 * \param pref pointer to the reference of the parent material.
 * \return 1 if found, 0 otherwise.
 *
 * Retrieve the starting domain reference (parent material) associated to the
 * split reference ref. Allow the call in non-multimaterial mode.
 *
 */
int MMG5_getStartRef(MMG5_pMesh mesh,MMG5_int ref,MMG5_int *pref) {
  MMG5_pInvMat pim;

  /* No multi-materials nor single material reference preservation */
  if( !mesh->info.nmat ) {
    *pref = 0;
    return 1;
  }

  /* Get parent of material */
  pim = &mesh->info.invmat;

  /* Return 0 if the material does not exist, 1 otherwise */
  if( !MMG5_InvMat_getParent(mesh,pim,ref,pref) )
    return 0;
  else
    return 1;
}

/**
 * \param mesh   pointer to the mesh structure.
 * \param pim    multimaterials inverse data table.
 *
 * Print materials lookup table.
 */
static void MMG5_InvMat_print(MMG5_pMesh mesh,MMG5_pInvMat pim) {
  MMG5_int ref,pref;

  /* Scan all references in the table limits, some may not exist */
  for( ref = pim->offset; ref < pim->offset + pim->size; ref++ ) {
    if( !MMG5_InvMat_getParent(mesh,pim,ref,&pref) ) continue;
    printf("%" MMG5_PRId " (%" MMG5_PRId "): %" MMG5_PRId " %d\n",ref,
           MMG5_InvMat_key(pim,ref),pref,
           MMG5_InvMat_getAttrib(pim,ref));
  }
}

/**
 * \param mesh   pointer to the mesh structure.
 * \return 1 if success, 0 if fail.
 *
 *
 * Initialize handling of multimaterial mode.
 *
 * An indexed table of materials has been provided by the MMG5_Mat datatype in
 * the form:
 *
 *   index | dospl       | ref       | rin       | rex
 *   -------------------------------------------------------
 *   0     | dospl_0     | ref_0     | rin_0     | rex_0
 *   ...   | ...         | ...       | ...       |
 *   k     | dospl_k     | ref_k     | rin_k     | rex_k
 *   ...   | ...         | ...       | ...       |
 *   n-1   | dospl_{n-1} | ref_{n-1} | rin_{n-1} | rex_{n-1}
 *
 * where dospl is the split/preserve attribute of the material, and rin,rex
 * are its child materials (if dospl). Viceversa, ref is the parent material
 * for rin,rex.
 *
 * Here a lookup table for material references is built through trivial hashing
 * for all references (both parent and child materials) with the key:
 *
 *   key = ref - ref_min,    ref_min = min_{k = 0,n-1} (ref_k, rin_k, rex_k)
 *
 * For all references, it is important to store
 * - the index k of the row in the original table,
 * - the characteristic attribute (parent-split, parent-preserve, child-interior,
 *   child-exterior) of the material.
 *
 * Since
 *   dospl = 0 or 1,  MG_PLUS = 2, and MG_MINUS = 3
 *
 * we can store the attribute as dospl (for the parent material) or MG_MINUS/
 * MG_PLUS (for the child materials), with its value ranging between 0 and 3.
 *
 * A convenient entry to store both the index and the attribute in the lookup
 * table is thus:
 *
 *   entry = 4*(index+1) + attribute
 *
 * leading to a lookup table in the form (the key ordering is only an example):
 *
 *   key   | entry
 *   ------------------------
 *   ...   | ...
 *   ref_k | 4*(k+1)+dospl_k
 *   ...   | ...
 *   rin_k | 4*(k+1)+MG_MINUS
 *   ...   | ...
 *   rex_k | 4*(k+1)+MG_PLUS
 *   ...   | ...
 *
 * where the index and the attribute of the material can be retrieved as
 *
 *   index     = entry / 4 -1
 *   attribute = entry % 4
 *
 *
 * What if two materials have the same reference?
 *  - child references should be distinct (the entry in the lookup table would
 *    be overwritten),
 *  - a parent material can have itself as child (a positive attribute would say
 *    it should be split, the attribute value would say if it is interior or
 *    exterior).
 * Thus, each of the maps parent->child_in and parent->child_ext is injective,
 * but a fixed point is allowed.
 *
 * Why child materials should be different?
 *  - because on failure it is necessary to recover the parent references and
 *    apply them on the tetrahedra to restore the starting materials
 *    distribution.
 *
 */
int MMG5_MultiMat_init(MMG5_pMesh mesh) {
  MMG5_pMat    pm;
  MMG5_pInvMat pim;
  int          k;
  MMG5_int     refmax,refmin;

  /* Nothing to do if no multi-material option */
  if( !mesh->info.nmat ) return 1;

  /* Error if all the materials have not been set */
  if( mesh->info.nmati < mesh->info.nmat ) {
    fprintf(stderr,"\n ## Error: %s: Only %d materials out of %d have been set.\n",
        __func__,mesh->info.nmati,mesh->info.nmat);
    return 0;
  }

  /* Get pointer to the structure */
  pim = &mesh->info.invmat;

  /* Initialize the max and min reference */
  refmax = 0;

  refmin = MMG5_INTMAX;

  /* Look for the max/min reference provided in material table */
  for( k = 0; k < mesh->info.nmat; k++ ) {
    pm = &mesh->info.mat[k];
    /* Update max and min val for original ref */
    if( pm->ref > refmax ) refmax = pm->ref;
    if( pm->ref < refmin ) refmin = pm->ref;
    if( !pm->dospl ) continue;
    /* Update max and min val with interior ref */
    if( pm->rin > refmax ) refmax = pm->rin;
    if( pm->rin < refmin ) refmin = pm->rin;
    /* Update max and min val with exterior ref */
    if( pm->rex > refmax ) refmax = pm->rex;
    if( pm->rex < refmin ) refmin = pm->rex;
  }

  /* Look for the max/min reference of tetra, triangles and edges provided
   * inside the mesh (worst case to avoid memory error when checking the
   * the inverse map). Looking at vertices is useless as
   * we will never check for the mapping of reference of vertices */
  for ( k=1; k<=mesh->ne; ++k ) {
    if( mesh->tetra[k].ref > refmax ) refmax = mesh->tetra[k].ref;
    if( mesh->tetra[k].ref < refmin ) refmin = mesh->tetra[k].ref;
  }
  for ( k=1; k<=mesh->nt; ++k ) {
    if( mesh->tria[k].ref > refmax ) refmax = mesh->tria[k].ref;
    if( mesh->tria[k].ref < refmin ) refmin = mesh->tria[k].ref;
  }
  for ( k=1; k<=mesh->na; ++k ) {
    if( mesh->edge[k].ref > refmax ) refmax = mesh->edge[k].ref;
    if( mesh->edge[k].ref < refmin ) refmin = mesh->edge[k].ref;
  }

  /* Get span of the lookup table */
  pim->offset = refmin;
  pim->size   = refmax - refmin + 1;
  assert( pim->size > 0 );

  /* Allocate lookup table */
  MMG5_ADD_MEM(mesh,pim->size*sizeof(int),"materials lookup table",return 0);
  MMG5_SAFE_CALLOC(pim->lookup,pim->size,int,return 0);

  /* Fill lookup table */
  for( k = 0; k < mesh->info.nmat; k++ ) {
    if( !MMG5_InvMat_set(mesh,pim,k) )
      return 0;
  }

  // MMG5_InvMat_print(mesh,pim);
  return 1;
}

/**
 * \param mesh   pointer to the mesh structure.
 * \param ref    initial reference.
 * \param refint internal reference after ls discretization.
 * \param refext external reference after ls discretization.
 * \return 1 if entity can be splitted, 0 if cannot be splitted.
 *
 * Identify whether an entity with reference ref should be split, and the
 * labels of the resulting entities.
 *
 */
int MMG5_isSplit(MMG5_pMesh mesh,MMG5_int ref,MMG5_int *refint,MMG5_int *refext) {
  MMG5_pInvMat pim;
  MMG5_pMat    pm;
  int          k;

  /* Default case: split with references MG_MINUS, MG_PLUS */
  if( !mesh->info.nmat ) {
    *refint = MG_MINUS;
    *refext = MG_PLUS;
    return 1;
  }

  /* Check in the info->mat table if reference ref is supplied by the user */
  pim = &mesh->info.invmat;
  k = MMG5_InvMat_getIndex(pim,ref);

  assert( k != -1 );
  pm = &mesh->info.mat[k];

  if ( !pm->dospl ) {
    return 0;
  } else {
    *refint = pm->rin;
    *refext = pm->rex;
    return 1;
  }
}

/**
 * \param mesh   pointer to the mesh structure.
 * \param ref    initial reference.
 * \return 1 if entity cannot be split, 0 if can be split.
 *
 * Identify whether an entity with reference ref should not be split.
 *
 */
int MMG5_isNotSplit(MMG5_pMesh mesh,MMG5_int ref) {
  MMG5_pInvMat pim;

  /* Split material by default if not in multi-material mode */
  if( !mesh->info.nmat ) return 0;

  /* Look in the table otherwise */
  pim = &mesh->info.invmat;
  if( !MMG5_InvMat_getAttrib(pim,ref) )
    return 0;
  else
    return 1;

}

/**
 * \param mesh   pointer to the mesh structure.
 * \param ref0   reference of the first tetrahedron sharing the face.
 * \param ref1   reference of the second tetrahedron sharing the face..
 * \return 1 if face is on the discrete level set, 0 if not.
 *
 * Identify whether a face is on the discrete level set or not.
 *
 */
int MMG5_isLevelSet(MMG5_pMesh mesh,MMG5_int ref0,MMG5_int ref1) {
  MMG5_pInvMat pim;
  int8_t       found0,found1;

  /* Check whether multimaterial case or not */
  if( mesh->info.nmat ) {
    /* Retrieve levelset information from the lookup table */
    pim = &mesh->info.invmat;
    found0 = MMG5_InvMat_getAttrib(pim,ref0);
    found1 = MMG5_InvMat_getAttrib(pim,ref1);

    if( (found0+found1) == (MG_MINUS+MG_PLUS) ) return 1;
    else return 0;

  } else {
    /* Single material, check references directly */
    if( ( ref0 == MG_MINUS && ref1 == MG_PLUS ) ||
        ( ref1 == MG_MINUS && ref0 == MG_PLUS ) ) return 1;
    else return 0;
  }
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set function.
 * \return 1 if success, 0 if fail.
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 */
int MMG5_snpval_ls(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTria       pt,pt1;
  MMG5_pPoint      p0;
  double           v1,v2,*tmp;
  MMG5_int         k,kk,iel,ns,nc,ip,ip1,ip2,npl,nmn;
  int              ilist;
  int8_t           i,j,j1,j2;
  MMG5_int         list[MMG5_TRIA_LMAX+2];

  /* Allocate memory for tmp */
  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(tmp,mesh->npmax+1,double,return 0);

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Snap values of sol that are close to 0 to 0 exactly */
  ns = nc = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if ( fabs(sol->m[k]) < MMG5_EPS ) {
      tmp[k] = sol->m[k];
      p0->flag = 1;
      sol->m[k] = 0.0;
      ns++;
    }
  }

  /* Check that the snapping process has not led to a nonmanifold situation */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      ip1 = pt->v[MMG5_inxt2[i]];
      ip2 = pt->v[MMG5_iprv2[i]];

      p0 = &mesh->point[ip];
      v1 = sol->m[ip1];
      v2 = sol->m[ip2];

      /* Catch a snapped point by a triangle where there is a sign change: use
       * the same convention than in ismaniball to evaluate sign changes. If
       * travelled in direct sense from a triangle, an edge is considered
       * without sign change if first vertex is 0. It has a sign change if
       * second vertex is 0 or if we have 2 vertices with different signs
       * (otherwise a 0 vertex leads to count 2 sign changes instead of one). */
      int smsgn = ((fabs(v2) < MMG5_EPS) || MG_SMSGN(v1,v2)) ? 1 : 0;
      if ( p0->flag && !smsgn ) {
        if ( !MMG5_ismaniball(mesh,sol,k,i) ) {
          if ( tmp[ip] < 0.0 )
            sol->m[ip] = -100.0*MMG5_EPS;
          else
            sol->m[ip] = 100.0*MMG5_EPS;
          nc++;
        }
        p0->flag = 0;
      }
    }
  }

  /* Check that the ls function does not show isolated spots with 0 values (without sign changes) */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      if ( fabs(sol->m[ip]) >= MMG5_EPS ) continue;
      npl = nmn = 0;
      int8_t opn; //unused
      ilist = MMG5_boulet(mesh,k,i,list,1,&opn);
      for(kk=0; kk<ilist; kk++) {
        iel = list[kk] / 3;
        j = list[kk] % 3;
        j1 = MMG5_inxt2[j];
        j2 = MMG5_iprv2[i];
        pt1 = &mesh->tria[iel];
        ip1 = pt1->v[j1];
        ip2 = pt1->v[j2];
        if ( sol->m[ip1] >= MMG5_EPS ) npl = 1;
        else if ( sol->m[ip1] <= -MMG5_EPS ) nmn = 1;

        if ( sol->m[ip2] >= MMG5_EPS ) npl = 1;
        else if ( sol->m[ip2] <= -MMG5_EPS ) nmn = 1;
      }

      if ( npl == 1 && nmn == 0 )
        sol->m[ip] = 100.0*MMG5_EPS;
      else if ( npl == 0 && nmn == 1 )
        sol->m[ip] = -100.0*MMG5_EPS;
    }
  }

  if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+nc > 0 )
    fprintf(stdout,"     %8" MMG5_PRId " points snapped, %" MMG5_PRId " corrected\n",ns,nc);

  /* memory free */
  MMG5_DEL_MEM(mesh,tmp);

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set values.
 * \param start index of the starting tria
 * \param istart local index (inside the tria \a start) of the vertex that we check.
 * \return 1 if success, 0 if fail
 *
 * Check whether snapping the value of vertex \a istart of \a start to 0 exactly
 * leads to a non manifold situation.
 *
 * \warning: we assume that the triangle \a start has vertex \a istart
 * with value 0 and the other two with changing values.
 *
 */
int MMG5_ismaniball(MMG5_pMesh mesh, MMG5_pSol sol, MMG5_int start, int8_t istart) {
  MMG5_pTria       pt;
  double           v1, v2;
  MMG5_int         refstart,*adja,k,ip1,ip2,end1;
  int8_t           i,i1,smsgn;
  static int8_t    mmgWarn=0;

  k = start;
  refstart = mesh->tria[k].ref;
  i = MMG5_inxt2[istart];

  /* First loop: stop if an external boundary, or a change in signs (or a 0) is met
     recall that MG_SMGSGN(a,b) = 1 provided a*b >0 */
  do{
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i] / 3;
    i1 = adja[i] % 3;
    i = MMG5_iprv2[i1];

    if ( k==0 ) break;

    pt = &mesh->tria[k];

    ip1 = pt->v[i1];
    ip2 = pt->v[i];

    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    if ( (fabs(v1) < MMG5_EPS) && (fabs(v2) < MMG5_EPS) ) {
      /* Do not authorize a snap that leads to a triangle with only 0 vertices */
      return 0;
    }

    /* Authorize change of references only provided the boundary reference is mesh->info.isoref */
    if ( pt->ref != refstart && pt->edg[i1] != mesh->info.isoref ) {
      smsgn = 0;
      k = 0;
    } else {
      /* Evaluation of sign change using following convention: If
       * travelled in direct sense from a triangle, an edge is considered
       * without sign change if first vertex is 0. It has a sign change if
       * second vertex is 0 or if we have 2 vertices with different signs
       * (otherwise a 0 vertex leads to count 2 sign changes instead of one). */
      smsgn = (fabs(v1) < MMG5_EPS) || ( (fabs(v2) > MMG5_EPS) && MG_SMSGN(v1,v2) ) ? 1 : 0;
    }
  }
  while ( smsgn && (k != start) );

  if ( k==start ) {
    /* Complete ball has been travelled without crossing a boundary or finding a
     * sign change: we are in the special case where v1 = v2 = v[istart] = 0 in
     * tria start. In this case, test MG_SMSGN(v1,v2) returns 0 while smsgn is
     * computed to 1, which is non consistent.  */
    assert ( smsgn );
    return 0;
  }

  end1 = k;
  k = start;
  i = MMG5_iprv2[istart];

  /* Second loop: same travel in the opposite sense */
  do{
    adja = &mesh->adja[3*(k-1)+1];
    k = adja[i] / 3;
    i1 = adja[i] % 3;
    i = MMG5_inxt2[i1];

    if ( k==0 ) break;

    pt = &mesh->tria[k];
    ip1 = pt->v[i1];
    ip2 = pt->v[i];

    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    if ( (fabs(v1) < MMG5_EPS) && (fabs(v2) < MMG5_EPS) ) {
      /* Do not authorize a snap that leads to a triangle with only 0 vertices */
      return 0;
    }

    if ( pt->ref != refstart && pt->edg[i1] != mesh->info.isoref ) {
      smsgn = 0;
      k = 0;
    } else {
      /* Evaluation of sign change using following convention: If
       * travelled in undirect sense from a triangle, an edge is considered
       * without sign change if second vertex is 0. It has a sign change if
       * first vertex is 0 or if we have 2 vertices with different signs
       * (it allows to evaluate the same splitted edges than the first loop). */
      smsgn = (fabs(v2) < MMG5_EPS) || ( (fabs(v1) > MMG5_EPS) && MG_SMSGN(v1,v2) ) ? 1 : 0;
    }
  }
  while ( smsgn && (k != start) );

  assert ( k!=start );

  /* If first stop was due to an external boundary, the second one must too
     (k==end1==0); else, the final triangle for the first travel must be that of
     the second one */
  if ( k != end1 ) {
    if ( !mmgWarn ) {
      mmgWarn = 1;
      fprintf(stderr,"\n  ## Warning: %s: unsnap at least 1 point "
              "(point %" MMG5_PRId " in tri %" MMG5_PRId ").\n",__func__,
              MMG5_indPt(mesh,mesh->tria[start].v[istart]),MMG5_indElt(mesh,start));
    }
    return 0;
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param ip0 First vertex of the triangle
 * \param ip1 Second vertex of the triangle
 * \param ip2 Third vertex of the triangle
 * \return area of the triangle
 *
 * Calculate the area of a triangle given by its vertices
 *
 **/
static inline
double MMG5_voltri(MMG5_pMesh mesh,MMG5_int ip0,MMG5_int ip1,MMG5_int ip2) {
  MMG5_pPoint    p0,p1,p2;
  double         vol;

  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  p2 = &mesh->point[ip2];

  vol = (p1->c[0]-p0->c[0])*(p2->c[1]-p0->c[1]) - (p1->c[1]-p0->c[1])*(p2->c[0]-p0->c[0]);
  vol = 0.5*fabs(vol);

  return vol;
}

/**
 * \param mesh pointer to the mesh structure
 * \param sol pointer to the ls function
 * \param k index of the triangle
 * \param pm 1 for computation of positive subdomain, -1 for negative one
 *
 * \return volfrac
 *
 * Calculate the area of the positive (if pm == 1) or negative (if pm == -1) subdomain
 * inside triangle k defined by the ls function in sol
 *
 **/
static inline
double MMG5_vfrac(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int k,int pm) {
  MMG5_pTria    pt;
  MMG5_pPoint   ppt[3];
  double        v[3],vfp,vfm,lam,area,eps,o1[2],o2[2];
  MMG5_int      ip[3],nplus,nminus,nzero;
  int8_t        i,i0,i1,i2,imin1,iplus1,iz;

  eps = MMG5_EPS*MMG5_EPS;
  pt = &mesh->tria[k];

  ip[0] = pt->v[0];
  ip[1] = pt->v[1];
  ip[2] = pt->v[2];

  ppt[0] = &mesh->point[ip[0]];
  ppt[1] = &mesh->point[ip[1]];
  ppt[2] = &mesh->point[ip[2]];

  v[0] = sol->m[ip[0]];
  v[1] = sol->m[ip[1]];
  v[2] = sol->m[ip[2]];

  /* Identify number of zero, positive and negative vertices, and corresponding indices */
  nplus = nminus = nzero = 0;
  imin1 = iplus1 = iz = -1;

  for (i=0; i<3; i++) {
    if ( fabs(v[i]) < eps ) {
      nzero++;
      if ( iz < 0 ) iz = i;
    }
    else if ( v[i] >= eps ) {
      nplus++;
      if ( iplus1 < 0 ) iplus1 = i;
    }
    else {
      nminus++;
      if ( imin1 < 0 ) imin1 = i;
    }
  }

  /* Degenerate case */
  if ( nzero == 3 ) return 0.0;

  /* Whole triangle is positive */
  if ( nminus == 0 ) {
    vfp = (ppt[1]->c[0]-ppt[0]->c[0])*(ppt[2]->c[1]-ppt[0]->c[1]) - (ppt[1]->c[1]-ppt[0]->c[1])*(ppt[2]->c[0]-ppt[0]->c[0]);
    vfp = 0.5*fabs(vfp);
    if ( pm == 1 ) return vfp;
    else           return 0.0;
  }

  /* Whole triangle is negative */
  if ( nplus == 0 ) {
    vfm = (ppt[1]->c[0]-ppt[0]->c[0])*(ppt[2]->c[1]-ppt[0]->c[1]) - (ppt[1]->c[1]-ppt[0]->c[1])*(ppt[2]->c[0]-ppt[0]->c[0]);
    vfm = 0.5*fabs(vfm);
    if ( pm == -1 ) return vfm;
    else            return 0.0;
  }

  /* Exactly one vertex is negative */
  if ( nminus == 1 ) {
    i0 = imin1;
    i1 = MMG5_inxt2[i0];
    i2 = MMG5_iprv2[i0];

    lam = v[i0] / (v[i0]-v[i1]);
    o1[0] = ppt[i0]->c[0] + lam*(ppt[i1]->c[0]-ppt[i0]->c[0]);
    o1[1] = ppt[i0]->c[1] + lam*(ppt[i1]->c[1]-ppt[i0]->c[1]);

    lam = v[i0] / (v[i0]-v[i2]);
    o2[0] = ppt[i0]->c[0] + lam*(ppt[i2]->c[0]-ppt[i0]->c[0]);
    o2[1] = ppt[i0]->c[1] + lam*(ppt[i2]->c[1]-ppt[i0]->c[1]);

    vfm = (o1[0]-ppt[i0]->c[0])*(o2[1]-ppt[i0]->c[1]) - (o1[1]-ppt[i0]->c[1])*(o2[0]-ppt[i0]->c[0]);
    vfm = 0.5*fabs(vfm);

    if ( pm == -1 ) return vfm;
    else {
      area = (ppt[1]->c[0]-ppt[0]->c[0])*(ppt[2]->c[1]-ppt[0]->c[1]) - (ppt[1]->c[1]-ppt[0]->c[1])*(ppt[2]->c[0]-ppt[0]->c[0]);
      area = 0.5*fabs(area);
      vfp = area-vfm;
      return vfp;
    }
  }

  /* Exactly one vertex is positive */
  if ( nplus == 1 ) {
    i0 = iplus1;
    i1 = MMG5_inxt2[i0];
    i2 = MMG5_iprv2[i0];

    lam = v[i0] / (v[i0]-v[i1]);
    o1[0] = ppt[i0]->c[0] + lam*(ppt[i1]->c[0]-ppt[i0]->c[0]);
    o1[1] = ppt[i0]->c[1] + lam*(ppt[i1]->c[1]-ppt[i0]->c[1]);

    lam = v[i0] / (v[i0]-v[i2]);
    o2[0] = ppt[i0]->c[0] + lam*(ppt[i2]->c[0]-ppt[i0]->c[0]);
    o2[1] = ppt[i0]->c[1] + lam*(ppt[i2]->c[1]-ppt[i0]->c[1]);

    vfp = (o1[0]-ppt[i0]->c[0])*(o2[1]-ppt[i0]->c[1]) - (o1[1]-ppt[i0]->c[1])*(o2[0]-ppt[i0]->c[0]);
    vfp = 0.5*fabs(vfp);

    if ( pm == 1 ) return vfp;
    else {
      area = (ppt[1]->c[0]-ppt[0]->c[0])*(ppt[2]->c[1]-ppt[0]->c[1]) - (ppt[1]->c[1]-ppt[0]->c[1])*(ppt[2]->c[0]-ppt[0]->c[0]);
      area = 0.5*fabs(area);
      vfm = area-vfp;
      return vfm;
    }
  }

  /* Should not pass here */
  return 0.0;
}

/** Return 1 if reference ref is in the base reference table, 0 otherwise */
static inline
int MMG5_isbr(MMG5_pMesh mesh,MMG5_int ref) {
  MMG5_int k;

  for(k=0; k<mesh->info.nbr; k++)
    if ( ref == mesh->info.br[k] ) return 1;

  return 0;
}

/**
 * \param mesh pointer to the mesh
 * \param sol pointer to the level-set
 *
 * \return 1 if success, 0 otherwise
 *
 * Removal of small parasitic components (bubbles of material, etc) with volume less than
 * mesh->info.rmc * volume of the mesh.
 *
 */
int MMG5_rmc(MMG5_pMesh mesh, MMG5_pSol sol){
  MMG5_pTria     pt,pt1,pt2;
  double         volc,voltot,v0,v1,v2;
  MMG5_int       k,kk,l,ll,ncp,ncm,ip0,ip1,ip2,cur,ipile,*pile,*adja,base;
  int8_t         i,i1,i2,onbr;

  ncp = 0;
  ncm = 0;

  /* Erase triangle flags */
  for (k=1; k<=mesh->nt; k++) mesh->tria[k].flag = 0;

  /* Calculate volume of the total mesh */
  voltot = 0.0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];
    voltot += MMG5_voltri(mesh,ip0,ip1,ip2);
  }

  /* Memory allocation for pile */
  MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_int),"temporary table",
               printf("  Exit program.\n");
               return 0);
  MMG5_SAFE_CALLOC(pile,mesh->nt+1,MMG5_int,return 0);

  /* Investigate only positive connected components */
  base = ++mesh->base;

  for (k=1; k<=mesh->nt; k++) {
    ipile = 0;
    volc = 0.0;
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    if ( pt->flag == base ) continue;

    /* Checks signs of the LS function at the 3 vertices of pt */
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];

    v0 = sol->m[ip0];
    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    if ( v0 <= 0.0 && v1 <= 0.0 && v2 <= 0.0 ) continue;

    /* Add triangle to pile if one vertex is > 0 */
    pt->flag = base;
    pile[ipile] = k;
    ipile++;
    if ( ipile > mesh->nt ) {
      fprintf(stderr,"\n  ## Problem in length of pile; function rmc.\n"
              " Check that the level-set intersect the mesh.\n"
              " Exit program.\n");

      return 0;
    }

    /* Pile up all the positive connected component attached to the first triangle */
    cur = 0;
    do {
      kk = pile[cur];
      pt1 = &mesh->tria[kk];

      /* Add local volume fraction of the positive subdomain to volc */
      volc += MMG5_vfrac(mesh,sol,kk,1);

      /* Add adjacent triangles to kk via positive vertices to the pile, if need be */
      adja = &mesh->adja[3*(kk-1)+1];
      for (i=0; i<3; i++) {
        ip0 = pt1->v[i];
        if ( sol->m[ip0] <= 0.0 ) continue;

        i1 = MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];

        /* First neighbor of positive vertex i */
        ll = adja[i1] / 3;
        if ( ll ) {
          pt2 = &mesh->tria[ll];
          if ( pt2->flag != base ) {
            pt2->flag = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->nt ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }

        /* Second neighbor of positive vertex i */
        ll = adja[i2] / 3;
        if ( ll ) {
          pt2 = &mesh->tria[ll];
          if ( pt2->flag != base ) {
            pt2->flag = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->nt ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }
      }
    }
    while ( ++cur < ipile );

    /* Remove connected component if its volume is too small */
    if ( volc < mesh->info.rmc*voltot ) {
      for (l=0; l<ipile; l++) {
        pt1 = &mesh->tria[pile[l]];
        for (i=0; i<3; i++) {
          ip0 = pt1->v[i];
          if ( sol->m[ip0] > 0.0 ) sol->m[ip0] = -100*MMG5_EPS;
        }
      }
      ncp++;
    }

  }

  /* Investigate only negative connected components */
  base = ++mesh->base;

  for (k=1; k<=mesh->nt; k++) {
    ipile = 0;
    volc = 0.0;
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    if ( pt->flag == base ) continue;

    /* Checks signs of the LS function at the 3 vertices of pt */
    ip0 = pt->v[0];
    ip1 = pt->v[1];
    ip2 = pt->v[2];

    v0 = sol->m[ip0];
    v1 = sol->m[ip1];
    v2 = sol->m[ip2];

    if ( v0 >= 0.0 && v1 >= 0.0 && v2 >= 0.0 ) continue;

    /* Pile up all the negative connected component attached to the first triangle */
    pt->flag = base;
    pile[ipile] = k;
    ipile++;
    if ( ipile > mesh->nt ) {
      fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
      return 0;
    }

    cur = 0;
    do {
      kk = pile[cur];
      pt1 = &mesh->tria[kk];

      /* Add local volume fraction of the negative subdomain to volc */
      volc += MMG5_vfrac(mesh,sol,kk,-1);

      /* Add adjacent triangles to kk via negative vertices to the pile, if need be */
      adja = &mesh->adja[3*(kk-1)+1];
      for (i=0; i<3; i++) {
        ip0 = pt1->v[i];
        if ( sol->m[ip0] >= 0.0 ) continue;

        i1= MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];

        /* First neighbor of negative vertex i */
        ll = adja[i1] / 3;
        if ( ll ) {
          pt2 = &mesh->tria[ll];
          if ( pt2->flag != base ) {
            pt2->flag = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->nt ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }

        /* Second neighbor of negative vertex i */
        ll = adja[i2] / 3;
        if ( ll ) {
          pt2 = &mesh->tria[ll];
          if ( pt2->flag != base ) {
            pt2->flag = base;
            pile[ipile] = ll;
            ipile++;
            if ( ipile > mesh->nt ) {
              fprintf(stderr,"\n  ## Problem in length of pile; function rmc. Exit program.\n");
              return 0;
            }
          }
        }

      }
    }
    while ( ++cur < ipile );

    /* Remove connected component if its volume is too small */
    if ( volc < mesh->info.rmc*voltot ) {
      for (l=0; l<ipile; l++) {
        pt1 = &mesh->tria[pile[l]];
        for (i=0; i<3; i++) {
          ip0 = pt1->v[i];
          if ( sol->m[ip0] < 0.0 ) sol->m[ip0] = 100*MMG5_EPS;
        }
      }
      ncm++;
    }

    /* Remove connected component if it is not attached to one base reference */
    if ( mesh->info.nbr ) {
      onbr = 0;
      for (l=0; l<ipile; l++) {
        pt1 = &mesh->tria[pile[l]];
        for (i=0; i<3; i++) {
          if ( MMG5_isbr(mesh,pt1->edg[i]) ) {
            i1 = MMG5_inxt2[i];
            i2 = MMG5_inxt2[i1];
            ip1 = pt1->v[i1];
            if ( sol->m[ip1] < 0.0 )  {
              onbr = 1;
              break;
            }
            ip2 = pt1->v[i2];
            if ( sol->m[ip2] < 0.0 )  {
              onbr = 1;
              break;
            }
          }
        }
        if ( onbr ) break;
      }

      if ( !onbr ) {
        for (l=0; l<ipile; l++) {
          pt1 = &mesh->tria[pile[l]];
          for (i=0; i<3; i++) {
            ip0 = pt1->v[i];
            if ( sol->m[ip0] < 0.0 ) sol->m[ip0] = 100*MMG5_EPS;
          }
        }
        ncm++;
      }
    }

  }

  /* Erase triangle flags */
  for (k=1; k<=mesh->nt; k++) mesh->tria[k].flag = 0;

  /* Release memory */
  MMG5_DEL_MEM(mesh,pile);

  if ( mesh->info.imprim > 0 || mesh->info.ddebug ) {
    printf("\n  *** Removed %" MMG5_PRId " positive parasitic bubbles and %" MMG5_PRId " negative parasitic bubbles\n",ncp,ncm);
  }

  return(1);
}

/**
 * \param mesh pointer to the mesh
 *
 * Reset mesh->info.isoref vertex and edge references to 0.
 *
 */
int MMG5_resetRef_ls(MMG5_pMesh mesh) {
  MMG5_pTria      pt;
  MMG5_pPoint     p0;
  MMG5_int        ref,k;
  int8_t          i;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;

    for (i=0; i<3; i++) {
      p0 = &mesh->point[pt->v[i]];
      if ( pt->edg[i] == mesh->info.isoref ) pt->edg[i] = 0;
      if ( p0->ref == mesh->info.isoref ) p0->ref = 0;
    }
  }

  /* Reset the triangle references to their initial distribution */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !pt->v[0] ) continue;
    if( !MMG5_getStartRef(mesh,pt->ref,&ref) ) return 0;
    pt->ref = ref;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param sol pointer to the level-set values.
 * \return 1.
 *
 * Set references to tris according to the sign of the level set function.
 *
 */
int MMG5_setref_ls(MMG5_pMesh mesh, MMG5_pSol sol) {
  MMG5_pTria    pt;
  double        v,v1;
  int           ier;
  MMG5_int      k,ip,ip1,ref,refint,refext;
  int8_t        i,i1,i2,nmn,npl,nz;

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    ref = pt->ref;
    nmn = npl = nz = 0;
    for (i=0; i<3; i++) {
      ip = pt->v[i];
      v = sol->m[ip];

      if ( v > 0.0 )
        npl++;
      else if ( v < 0.0 )
        nmn++;
      else
        nz++;
    }

    assert(nz < 3);

    /* Keep the initial triangle references of the mesh if iso==2, set
     * positive and negative ls refs otherwise */
    if ( mesh->info.iso != 2 ) {

      /* find if current reference should be splitted and the new positive and negative refs */
      ier = MMG5_isSplit(mesh,ref,&refint,&refext);
      if ( ier ) {
        if ( npl ) {
          assert( !nmn );
          pt->ref = refext;
        }
        else {
          assert ( nmn );
          pt->ref = refint;
        }
      }
    }

    /* Set mesh->info.isoref ref at ls edges and at the points of these edges */
    if ( nz == 2 ) {
      for (i=0; i<3; i++) {
        ip  = pt->v[MMG5_inxt2[i]];
        ip1 = pt->v[MMG5_iprv2[i]];
        v   = sol->m[ip];
        v1  = sol->m[ip1];
        if ( v == 0.0 && v1 == 0.0) {
          pt->edg[i]  = mesh->info.isoref;
          pt->tag[i] |= MG_REF;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[pt->v[i1]].ref = mesh->info.isoref;
          mesh->point[pt->v[i2]].ref = mesh->info.isoref;
        }
      }
    }

  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param start index of starting tria.
 * \param istart local index of point that we check (in tria \a start)
 * \return 1 if the ball is manifold, 0 otherwise.
 *
 * Check whether the ball of vertex i in tria start is manifold;
 *
 * \warning inxt[i] is an edge belonging to the implicit boundary.
 *
 */
int MMG5_chkmaniball(MMG5_pMesh mesh, MMG5_int start, int8_t istart) {
  MMG5_int           refstart,*adja,k;
  int8_t             i,i1;

  k = start;
  i = istart;

  i1 = MMG5_iprv2[i];


  MMG5_pTria pt = &mesh->tria[start];
  assert( MG_EDG(pt->tag[i1]) && (pt->edg[i1]==mesh->info.isoref) );

  /** Step 1: Travel while another part of the implicit boundary is not met */
  refstart = pt->ref;
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;

    if ( !k ) break;

    if ( mesh->info.iso !=2 ) {
      /** Normal or multi-material mode: check for change in triangle references */
      if ( mesh->tria[k].ref != refstart) break;
    }
    else {
      /** Input reference preservation mode (mmgs --keep-ref option): Check if
       * we cross an isoref edge */
      if ( mesh->tria[k].edg[i]==mesh->info.isoref ) break;
    }
    i = MMG5_inxt2[i];
  }
  while ( k!=start );

  assert(k!=start); //unexpected case

  /** Step 2: Check why the loop has ended */
  /** a./ Case where a boundary is hit: travel in the other sense from start, and
      make sure that a boundary is hit too */
  if ( k == 0 ) {
    k = start;
    i = istart;

    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_iprv2[i];
    k = adja[i1] / 3;
    i = adja[i1] % 3;
    i = MMG5_iprv2[i];

    /** Ball is manifold if tested point is connected to two external edges */
    if ( k == 0 ) return 1;

    do {
      adja = &mesh->adja[3*(k-1)+1];
      i1 = MMG5_iprv2[i];

      k = adja[i1] / 3;
      i = adja[i1] % 3;

      if ( !k ) break;

      if ( mesh->info.iso !=2 ) {
        /* Normal or multi-material mode: check for change in triangle references */
        if ( mesh->tria[k].ref == refstart) break;
      }
      else {
        /* Input reference preservation mode (mmgs --keep-ref option): Check if
         * we cross an isoref edge */
        if ( mesh->tria[k].edg[i]==mesh->info.isoref ) break;
      }
      i = MMG5_iprv2[i];
    }
    while ( k!=start );

    assert(k!=start); //unexpected case

    return !k;
  }

  /** b./ General case: go on travelling until another implicit boundary is met */
  i = MMG5_inxt2[i];
  do {
    adja = &mesh->adja[3*(k-1)+1];
    i1 = MMG5_inxt2[i];

    k = adja[i1] / 3;
    i = adja[i1] % 3;

    if ( !k ) break;

    if ( mesh->info.iso !=2 ) {
      /* Check tria ref change */
      if ( mesh->tria[k].ref == refstart) break;
    }
    else {
      /* Check if we cross an isoref edge */
      if ( mesh->tria[k].edg[i]==mesh->info.isoref ) break;
    }

    i = MMG5_inxt2[i];
  }
  while ( k!=start );

  /** Ball is non-manifold if at least 3 boundary segments meeting at p */
  if ( k != start )
    return 0;

  return 1;
}

/**
 * \param mesh pointer to the mesh.
 * \return 1 if the mesh is manifold, 0 otherwise.
 *
 * Check whether the resulting two subdomains occupying mesh are manifold.
 *
 */
int MMG5_chkmanimesh(MMG5_pMesh mesh) {
  MMG5_pTria      pt,pt1;
  MMG5_int        *adja,k;
  MMG5_int        cnt,iel;
  int8_t          i,i1;
  static int8_t   mmgWarn = 0;

  /** First check: check whether one triangle in the mesh has 3 boundary faces */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    adja = &mesh->adja[3*(k-1)+1];
    cnt = 0;
    for (i=0; i<3; i++) {
      iel = adja[i] / 3;

      if (!iel ) {
        cnt++;
        continue;
      }
      else {
        if ( mesh->info.iso !=2 ) {
          /* Multi-material mode may lead to have only 1 isoref edge around a
             point (due to nosplit option): check tria ref change */
          pt1 = &mesh->tria[iel];
          if ( pt1->ref != pt->ref ) cnt++;
        }
        else {
          /* keep-ref mode: check if isoref edge */
          if ( pt->edg[i] == mesh->info.isoref ) cnt++;
        }
      }
    }
    if( cnt == 3 ) {
      if ( !mmgWarn ) {
        mmgWarn = 1;
        fprintf(stderr,"\n  ## Warning: %s: at least 1 triangle with 3 boundary"
                " edges.\n",__func__);
      }
    }
  }

  /** Second check: check whether the configuration is manifold in the ball of
     each point; each vertex on the implicit boundary is caught in such a way
     that i1 inxt[i1] is one edge of the implicit boundary */
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for (i=0; i<3; i++) {
      adja = &mesh->adja[3*(k-1)+1];
      iel = adja[i] / 3;

      if (! iel ) continue;

      if ( mesh->info.iso !=2 ) {
        /* Check change of tria ref */
        pt1 = &mesh->tria[iel];
        if ( pt->ref == pt1->ref || pt->edg[i]!= mesh->info.isoref ) continue;
      }
      else {
        /* Check isoref edge only */
        if ( pt->edg[i] != mesh->info.isoref ) continue;
      }

      i1 = MMG5_inxt2[i];
      if ( !MMG5_chkmaniball(mesh,k,i1) ) {
        fprintf(stderr,"   *** Topological problem\n");
        fprintf(stderr,"       non manifold curve at point %" MMG5_PRId "\n",pt->v[i1]);
        fprintf(stderr,"       non manifold curve at tria %" MMG5_PRId " (ip %d)\n", MMG5_indElt(mesh,k),i1);
        return 0;
      }
    }
  }

  if ( mesh->info.imprim > 0 || mesh->info.ddebug )
    fprintf(stdout,"  *** Manifold implicit surface.\n");

  return 1;
}
