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
 * \brief Functions used in mmg applications
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 04 2015
 * \copyright GNU Lesser General Public License.
 **/

#include "mmgcommon_private.h"

/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param node pointer toward a MMG5_iNode (cell for linked list)
 * \return 1 if we can alloc the node \a node, 0 otherwise.
 *
 * Node allocation.
 *
 */
static inline
int MMG5_Alloc_inode( MMG5_pMesh mesh, MMG5_iNode **node ) {

  MMG5_ADD_MEM(mesh,sizeof(MMG5_iNode),"boundary reference node",
                return 0;);

  MMG5_SAFE_MALLOC(*node,1,MMG5_iNode,return 0);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param liLi pointer toward the address of the root of the linked list.
 * \param val value to add to the linked list.
 * \return 1 if the node is inserted, 0 if the node is not inserted, -1 if fail.
 *
 * Add a node with value \a val to a sorted linked list with unique entries.
 *
 * \remark as the linked list had unique entries, we don't insert a node if it
 * exists.
 *
 */
int MMG5_Add_inode( MMG5_pMesh mesh, MMG5_iNode **liLi, int val ) {
  MMG5_iNode  *newNode, *cur;

  cur = *liLi;

  /* Travel through the linked list and search if the value val exist or, if
   * not, where to insert it */
  if ( cur ) {
    if ( val < (*liLi)->val ) {
      /* Add a value at the list head */
      if ( !MMG5_Alloc_inode(mesh,&newNode) ) return -1;

      newNode->val = val;
      newNode->nxt = (*liLi);

      (*liLi) = newNode;

      return 1;

    }
    else if (val == (*liLi)->val ) return 0;

    while ( cur->nxt && ( val >= (cur->nxt)->val) )
      cur = cur->nxt;

    if ( val == cur->val ) return 0;

    if ( !MMG5_Alloc_inode(mesh,&newNode) ) return -1;

    newNode->val = val;
    newNode->nxt = cur->nxt;
    cur->nxt = newNode;
  }
  else {
    if ( !MMG5_Alloc_inode(mesh,&newNode) ) return -1;

    newNode->val = val;
    newNode->nxt = NULL;

    *liLi = newNode;
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure (for count of used memory).
 * \param liLi pointer toward the root of the linked list.
 *
 * Free the memory used by the linked list whose root is \a liLi.
 *
 */
void MMG5_Free_ilinkedList( MMG5_pMesh mesh, MMG5_iNode *liLi ) {
  MMG5_iNode *cur,*nxt;

  cur = liLi;
  while (cur) {
    nxt = cur;
    cur = cur->nxt;

    MMG5_DEL_MEM(mesh,nxt);
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdryRefs pointer toward the list of the boundary references.
 * \return npar, the number of local parameters at triangles if success,
 * 0 otherwise.
 *
 * Count the local default values at triangles and fill the list of the boundary
 * references.
 *
 */
int MMG5_countLocalParamAtTri( MMG5_pMesh mesh,MMG5_iNode **bdryRefs) {
  MMG5_int    k;
  int         npar,ier;

  /** Count the number of different boundary references and list it */
  (*bdryRefs) = NULL;

  k = mesh->nt? mesh->tria[1].ref : 0;

  /* Try to alloc the first node */
  ier = MMG5_Add_inode( mesh, bdryRefs, k );
  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate the first boundary"
           " reference node.\n",__func__);
    return 0;
  }
  else {
    assert(ier);
    npar = 1;
  }

  for ( k=1; k<=mesh->nt; ++k ) {
    ier = MMG5_Add_inode( mesh, bdryRefs, mesh->tria[k].ref );

    if ( ier < 0 ) {
      printf("  ## Warning: %s: unable to list the tria references."
             " Uncomplete parameters file.\n",__func__ );
      break;
    }
    else if ( ier ) ++npar;
  }

  return npar;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param bdryRefs pointer toward the list of the boundary references.
 * \param out pointer toward the file in which to write.
 * \return 1 if success, 0 otherwise.
 *
 * Write the local default values at triangles in the parameter file.
 *
 */
int MMG5_writeLocalParamAtTri( MMG5_pMesh mesh, MMG5_iNode *bdryRefs,
                                FILE *out ) {
  MMG5_iNode *cur;

  cur = bdryRefs;
  while( cur ) {
    fprintf(out,"%"MMG5_PRId" Triangle %e %e %e \n",cur->val,
            mesh->info.hmin, mesh->info.hmax,mesh->info.hausd);
    cur = cur->nxt;
  }

  MMG5_Free_ilinkedList(mesh,bdryRefs);

  return 1;
}
