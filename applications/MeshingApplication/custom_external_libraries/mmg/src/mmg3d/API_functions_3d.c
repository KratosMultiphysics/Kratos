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
 * \file mmg3d/API_functions_3d.c
 * \brief C API functions definitions for MMG3D library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 * \warning Use the MMG3D_ prefix: MMG5_ prefix will became obsolete soon...
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmg3d/libmmg3d.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * C API for MMG3D library.
 *
 */

#include "libmmg3d.h"
#include "libmmg3d_private.h"
#include "inlined_functions_3d_private.h"
#include "PRoctree_3d_private.h"

int MMG3D_Init_mesh(const int starter,...) {
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = MMG3D_Init_mesh_var(argptr);

  va_end(argptr);

  return ier;
}
void MMG3D_Init_fileNames(MMG5_pMesh mesh,MMG5_pSol sol
  ) {

  MMG5_Init_fileNames(mesh,sol);
  return;
}

int MMG3D_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin) {

  return MMG5_Set_inputMeshName(mesh,meshin);
}

int MMG3D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin) {
  return MMG5_Set_inputSolName(mesh,sol,solin);
}

int MMG3D_Set_inputParamName(MMG5_pMesh mesh, const char* fparamin) {
  return MMG5_Set_inputParamName(mesh,fparamin);
}

int MMG3D_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout) {

  return MMG5_Set_outputMeshName(mesh,meshout);
}

int MMG3D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout) {
  return MMG5_Set_outputSolName(mesh,sol,solout);
}

void MMG3D_Init_parameters(MMG5_pMesh mesh) {

  /* Init common parameters for mmgs and mmg3d. */
  MMG5_Init_parameters(mesh);

  /* default values for integers */
  mesh->info.lag      = MMG5_LAG;
  mesh->info.setfem   = MMG5_FEM;
  mesh->info.optim    = MMG5_OFF;
  mesh->info.optimLES = MMG5_OFF;
  /* [0/1]    ,avoid/allow surface modifications */
  mesh->info.nosurf   =  MMG5_OFF;
#ifdef USE_SCOTCH
   /* [1/0]    , Turn on/off the renumbering using SCOTCH */
  mesh->info.renum    = MMG5_ON;
#else
   /* [0]    , Turn on/off the renumbering using SCOTCH */
  mesh->info.renum    = MMG5_OFF;
#endif

  /* default values for doubles */
  /* level set value */
  mesh->info.ls       = MMG5_LS;

#ifndef MMG_PATTERN
  mesh->info.PROctree = MMG5_PROCTREE;
#endif
}

int MMG3D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, MMG5_int np, int typSol) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && sol->m )
    fprintf(stderr,"\n  ## Warning: %s: old solution deletion.\n",__func__);

  if ( typEntity != MMG5_Vertex ) {
    fprintf(stderr,"\n  ## Error: %s: mmg3d need a solution imposed on vertices.\n",
            __func__);
    return 0;
  }

  sol->type = typSol;

  if ( typSol == MMG5_Scalar ) {
    sol->size = 1;
  }
  else if ( typSol == MMG5_Vector ) {
    sol->size = 3;
  }
  else if ( typSol == MMG5_Tensor ) {
    sol->size = 6;
    /* User will provide its own metric: classical storage at ridges */
    mesh->info.metRidTyp = 0;
  }
  else {
    fprintf(stderr,"\n  ## Error: %s: type of solution not yet implemented.\n",
            __func__);
    return 0;
  }

  sol->dim = 3;
  if ( np ) {
    sol->np  = np;
    sol->npi = np;
    if ( sol->m )
      MMG5_DEL_MEM(mesh,sol->m);

    sol->npmax = mesh->npmax;
    MMG5_ADD_MEM(mesh,(sol->size*(sol->npmax+1))*sizeof(double),"initial solution",
                  return 0);
    MMG5_SAFE_CALLOC(sol->m,(sol->size*(sol->npmax+1)),double,return 0);
  }
  return 1;
}

int MMG3D_Set_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol *sol,int nsols,
                                 MMG5_int nentities, int *typSol) {
  MMG5_pSol psl;
  char      data[18];
  int       j;

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && mesh->nsols ) {
    if ( *sol ) {
      fprintf(stderr,"\n  ## Warning: %s: old solutions array deletion.\n",
              __func__);
      MMG5_DEL_MEM(mesh,*sol);
    }
  }

  /** Sol tab allocation */
  mesh->nsols = nsols;

  MMG5_ADD_MEM(mesh,nsols*sizeof(MMG5_Sol),"solutions array",
                return 0);
  MMG5_SAFE_CALLOC(*sol,nsols,MMG5_Sol,return 0);

  for ( j=0; j<nsols; ++j ) {
    psl = *sol + j;
    psl->ver = 2;

    /* Give an arbitrary name to the solution */
    sprintf(data,"sol_%d",j);
    if ( !MMG3D_Set_inputSolName(mesh,psl,data) ) {
      return 0;
    }
    /* Give an arbitrary name to the solution */
    sprintf(data,"sol_%d.o",j);
    if ( !MMG3D_Set_outputSolName(mesh,psl,data) ) {
      return 0;
    }

    if ( !MMG3D_Set_solSize(mesh,psl,MMG5_Vertex,mesh->np,typSol[j]) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to set the size of the"
              " solution num %d.\n",__func__,j);
      return 0;
    }
  }
  return 1;
}

/**
 * \param mesh pointer to the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nprism number of prisms.
 * \param nt number of triangles.
 * \param nquad number of quads.
 * \param na number of edges.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Check the input mesh size and assign their values to the mesh.
 *
 */
int MMG3D_setMeshSize_initData(MMG5_pMesh mesh, MMG5_int np, MMG5_int ne, MMG5_int nprism,
                               MMG5_int nt, MMG5_int nquad, MMG5_int na ) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->tria || mesh->tetra || mesh->edge) )
    fprintf(stderr,"\n  ## Warning: %s: old mesh deletion.\n",__func__);

  if ( !np ) {
    fprintf(stderr,"  ** MISSING DATA:\n");
    fprintf(stderr,"     Your mesh must contains at least points.\n");
    return 0;
  }
  if ( !ne && (mesh->info.imprim > 4 || mesh->info.ddebug) ) {
    fprintf(stderr,"  ** WARNING:\n");
    fprintf(stderr,"     Your mesh don't contains tetrahedra.\n");
  }
  if ( mesh->point )
    MMG5_DEL_MEM(mesh,mesh->point);
  if ( mesh->tetra )
    MMG5_DEL_MEM(mesh,mesh->tetra);
  if ( mesh->prism )
    MMG5_DEL_MEM(mesh,mesh->prism);
  if ( mesh->tria )
    MMG5_DEL_MEM(mesh,mesh->tria);
  if ( mesh->quadra )
    MMG5_DEL_MEM(mesh,mesh->quadra);
  if ( mesh->edge )
    MMG5_DEL_MEM(mesh,mesh->edge);

  mesh->np  = np;
  mesh->ne  = ne;
  mesh->nt  = nt;
  mesh->na  = na;
  mesh->nprism = nprism;
  mesh->nquad  = nquad;

  mesh->npi = mesh->np;
  mesh->nei = mesh->ne;
  mesh->nti = mesh->nt;
  mesh->nai = mesh->na;

  return 1;
}

int MMG3D_Set_meshSize(MMG5_pMesh mesh, MMG5_int np, MMG5_int ne, MMG5_int nprism,
                       MMG5_int nt, MMG5_int nquad, MMG5_int na ) {

  /* Check input data and set mesh->ne/na/np/nt to the suitable values */
  if ( !MMG3D_setMeshSize_initData(mesh,np,ne,nprism,nt,nquad,na) )
    return 0;

  /* Check the -m option */
  if( mesh->info.mem > 0) {
    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne)) {
      if ( !MMG3D_memOption(mesh) )  return 0;
    } else if(mesh->info.mem < 39) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory  %d\n",__func__,
              mesh->info.mem);
      return 0;
    }
  } else {
    if ( !MMG3D_memOption(mesh) )  return 0;
  }

  /* Mesh allocation and linkage */
  if ( !MMG3D_setMeshSize_alloc( mesh ) ) return 0;

  return 1;
}

int MMG3D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, MMG5_int* np, int* typSol) {

  if ( typEntity != NULL )
    *typEntity = MMG5_Vertex;

  if ( typSol != NULL ) {
    if ( sol->size == 1 )
      *typSol    = MMG5_Scalar;
    else if ( sol->size == 3 )
      *typSol    = MMG5_Vector;
    else if ( sol->size == 6 )
      *typSol    = MMG5_Tensor;
    else
      *typSol    = MMG5_Notype;
  }

  assert( (!sol->np) || (sol->np == mesh->np));

  if ( np != NULL )
    *np = sol->np;

  return 1;
}

int MMG3D_Get_solsAtVerticesSize(MMG5_pMesh mesh, MMG5_pSol *sol, int *nsols,
                                 MMG5_int* np, int* typSol) {
  MMG5_pSol      psl;
  MMG5_int       j;

  if ( !mesh ) {
    fprintf(stderr,"\n  ## Error: %s: your mesh structure must be allocated"
            " and filled\n",__func__);
    return 0;
  }

  if ( nsols != NULL )
    *nsols = mesh->nsols;

  for ( j=0; j<mesh->nsols; ++j ) {
    psl = *sol + j;

    if ( typSol != NULL ) {
      typSol[j]    = psl->type;
    }

    assert( (!psl->np) || (psl->np == mesh->np));
  }
  if ( np != NULL )
    *np = mesh->np;

  return 1;
}

int MMG3D_Get_meshSize(MMG5_pMesh mesh, MMG5_int* np, MMG5_int* ne, MMG5_int* nprism,
                       MMG5_int* nt, MMG5_int * nquad, MMG5_int* na) {

  if ( np != NULL )
    *np = mesh->np;
  if ( ne != NULL )
    *ne = mesh->ne;
  if ( nprism != NULL )
    *nprism = mesh->nprism;
  if ( nt != NULL )
    *nt = mesh->nt;
  if ( nquad != NULL )
    *nquad = mesh->nquad;
  if ( na != NULL )
    *na = mesh->na;

  return 1;
}

int MMG3D_Set_vertex(MMG5_pMesh mesh, double c0, double c1, double c2, MMG5_int ref, MMG5_int pos) {

  if ( !mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: you must set the number of points with the",
            __func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting vertices in mesh.\n");
    return 0;
  }

  if ( pos > mesh->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new point.\n",__func__);
    fprintf(stderr,"    max number of points: %" MMG5_PRId "\n",mesh->npmax);
    MMG5_INCREASE_MEM_MESSAGE();
    return 0;
  }

  if ( pos > mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new vertex at position %" MMG5_PRId ".",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of vertices: %" MMG5_PRId "\n",mesh->np);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the vertex.\n");
    return 0;
  }

  mesh->point[pos].c[0] = c0;
  mesh->point[pos].c[1] = c1;
  mesh->point[pos].c[2] = c2;
  mesh->point[pos].ref  = ref;
  mesh->point[pos].tag  = MG_NUL;
  mesh->point[pos].flag = 0;
  mesh->point[pos].tmp = 0;

  return 1;
}

int MMG3D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
                     int* isCorner, int* isRequired) {

  if ( mesh->npi == mesh->np ) {
    mesh->npi = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_vertex function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",mesh->np);
    }
  }

  mesh->npi++;

  if ( mesh->npi > mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get point.\n",__func__);
    fprintf(stderr,"     The number of call of MMG3D_Get_vertex function");
    fprintf(stderr," can not exceed the number of points: %" MMG5_PRId "\n ",mesh->np);
    return 0;
  }

  return MMG3D_GetByIdx_vertex( mesh,c0,c1,c2,ref,isCorner,isRequired,mesh->npi);
}

int MMG3D_GetByIdx_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, MMG5_int* ref,
                          int* isCorner, int* isRequired, MMG5_int idx) {

  if ( idx < 1 || idx > mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get point at position %" MMG5_PRId ".\n",
            __func__,idx);
    fprintf(stderr,"     Your vertices numbering goes from 1 to %" MMG5_PRId "\n",mesh->np);
    return 0;
  }

  *c0  = mesh->point[idx].c[0];
  *c1  = mesh->point[idx].c[1];
  *c2  = mesh->point[idx].c[2];
  if ( ref != NULL )
    *ref = mesh->point[idx].ref;

  if ( isCorner != NULL ) {
    if ( mesh->point[idx].tag & MG_CRN )
      *isCorner = 1;
    else
      *isCorner = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->point[idx].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}

int  MMG3D_Set_vertices(MMG5_pMesh mesh, double *vertices,MMG5_int *refs) {

  MMG5_pPoint ppt;
  MMG5_int    i,j;

  /*coordinates vertices*/
  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*3;
    ppt->c[0]  = vertices[j];
    ppt->c[1]  = vertices[j+1];
    ppt->c[2]  = vertices[j+2];

    ppt->tag = MG_NUL;
    ppt->flag = 0;
    ppt->tmp = 0;

    if ( refs != NULL )
      ppt->ref   = refs[i-1];
  }

  return 1;
}


int  MMG3D_Get_vertices(MMG5_pMesh mesh, double* vertices, MMG5_int* refs,
                        int* areCorners, int* areRequired) {
  MMG5_pPoint ppt;
  MMG5_int    i,j;

  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*3;
    vertices[j] = ppt->c[0];
    vertices[j+1] = ppt->c[1];
    vertices[j+2] = ppt->c[2];

    j = i-1;
    if ( refs != NULL )
      refs[j] = ppt->ref;

    if ( areCorners !=NULL ) {
      if ( ppt->tag & MG_CRN )
        areCorners[j] = 1;
      else
        areCorners[j] = 0;
    }

    if ( areRequired != NULL ) {
      if ( ppt->tag & MG_REQ )
        areRequired[j] = 1;
      else
        areRequired[j] = 0;
    }
  }

  return 1;
}

int MMG3D_Set_tetrahedron(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2, MMG5_int v3, MMG5_int ref, MMG5_int pos) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  double      vol;
  MMG5_int    aux;
  int         j,ip;

  if ( !mesh->ne ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of elements with the",
            __func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting elements in mesh\n");
    return 0;
  }

  if ( pos > mesh->nemax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new element.\n",
            __func__);
    fprintf(stderr,"    max number of element: %" MMG5_PRId "\n",mesh->nemax);
    MMG5_INCREASE_MEM_MESSAGE();
    return 0;
  }

  if ( pos > mesh->ne ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new tetrahedron at position %" MMG5_PRId ".",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of tetrahedron: %" MMG5_PRId "\n",mesh->ne);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the tetrahedron.\n");
    return 0;
  }

  pt = &mesh->tetra[pos];
  pt->v[0] = v0;
  pt->v[1] = v1;
  pt->v[2] = v2;
  pt->v[3] = v3;
  pt->ref  = MMG5_abs(ref);

  mesh->point[pt->v[0]].tag &= ~MG_NUL;
  mesh->point[pt->v[1]].tag &= ~MG_NUL;
  mesh->point[pt->v[2]].tag &= ~MG_NUL;
  mesh->point[pt->v[3]].tag &= ~MG_NUL;

  vol = MMG5_orvol(mesh->point,pt->v);
  if ( fabs(vol) <= MMG5_EPSD2 ) {
    fprintf(stderr,"\n  ## Error: %s: tetrahedron %" MMG5_PRId " has volume null.\n",
            __func__,pos);
    for ( ip=0; ip<4; ip++ ) {
      ppt = &mesh->point[pt->v[ip]];
      for ( j=0; j<3; j++ ) {
        if ( fabs(ppt->c[j])>0. ) {
          fprintf(stderr," Check that you don't have a sliver tetrahedron.\n");
          return 0;
        }
      }
    }
    fprintf(stderr,"  All vertices have zero coordinates.");
    fprintf(stderr," Check that you have set the vertices before the tetrahedra.\n");
    return 0;
  }
  else if ( vol < 0.0 ) {
    /* Possibly switch 2 vertices number so that each tet is positively oriented */
    aux = pt->v[2];
    pt->v[2] = pt->v[3];
    pt->v[3] = aux;
    /* mesh->xt temporary used to count reoriented tetra */
    mesh->xt++;
  }

  return 1;
}

int MMG3D_Get_tetrahedron(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* v3,
                          MMG5_int* ref, int* isRequired) {

  if ( mesh->nei == mesh->ne ) {
    mesh->nei = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of"
              " tetrahedra.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_tetrahedron function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of tetrahedron: %" MMG5_PRId "\n ",mesh->ne);
    }
  }

  mesh->nei++;

  if ( mesh->nei > mesh->ne ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get tetra.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_tetrahedron function");
    fprintf(stderr," can not exceed the number of tetra: %" MMG5_PRId "\n ",mesh->ne);
    return 0;
  }

  *v0  = mesh->tetra[mesh->nei].v[0];
  *v1  = mesh->tetra[mesh->nei].v[1];
  *v2  = mesh->tetra[mesh->nei].v[2];
  *v3  = mesh->tetra[mesh->nei].v[3];
  if ( ref != NULL ) {
    *ref = mesh->tetra[mesh->nei].ref;
  }

  if ( isRequired != NULL ) {
    if ( mesh->tetra[mesh->nei].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}

int  MMG3D_Set_tetrahedra(MMG5_pMesh mesh, MMG5_int *tetra, MMG5_int *refs) {
  MMG5_pPoint ppt;
  MMG5_pTetra pt;
  double      vol;
  int         ip;
  MMG5_int    aux,i,j;

  mesh->xp = 0;
  for (i=1;i<=mesh->ne;i++)
  {
    j = (i-1)*4;
    pt = &mesh->tetra[i];
    pt->v[0]  = tetra[j];
    pt->v[1]  = tetra[j+1];
    pt->v[2]  = tetra[j+2];
    pt->v[3]  = tetra[j+3];

    if ( refs != NULL )
      pt->ref   = MMG5_abs(refs[i-1]);

    mesh->point[pt->v[0]].tag &= ~MG_NUL;
    mesh->point[pt->v[1]].tag &= ~MG_NUL;
    mesh->point[pt->v[2]].tag &= ~MG_NUL;
    mesh->point[pt->v[3]].tag &= ~MG_NUL;

    vol = MMG5_orvol(mesh->point,pt->v);

    if ( fabs(vol) <= MMG5_EPSD2 ) {
      fprintf(stderr,"\n  ## Error: %s: tetrahedron %" MMG5_PRId " has volume null.\n",
              __func__,i);

      for ( ip=0; ip<4; ip++ ) {
        ppt = &mesh->point[pt->v[ip]];
        for ( j=0; j<3; j++ ) {
          if ( fabs(ppt->c[j])>0. ) {
            fprintf(stderr," Check that you don't have a sliver tetrahedron.\n");
            return 0;
          }
        }
      }

      fprintf(stderr,"  All vertices have zero coordinates.");
      fprintf(stderr," Check that you have set the vertices before the tetrahedra.\n");
      return 0;
    }
    else if ( vol < 0.0 ) {
      /* Possibly switch 2 vertices number so that each tet is positively oriented */
      aux = pt->v[2];
      pt->v[2] = pt->v[3];
      pt->v[3] = aux;

      ++mesh->xp;
    }
  }

  return 1;
}

int  MMG3D_Get_tetrahedra(MMG5_pMesh mesh, MMG5_int *tetra, MMG5_int *refs, int * areRequired) {
  MMG5_pTetra pt;
  MMG5_int    i, j;

  for (i=1;i<=mesh->ne;i++)
  {
    j = (i-1)*4;
    pt = &mesh->tetra[i];
    tetra[j]   = pt->v[0];
    tetra[j+1] = pt->v[1];
    tetra[j+2] = pt->v[2];
    tetra[j+3] = pt->v[3];
    if ( refs!=NULL )
      refs[i-1]  = pt->ref ;
    if ( areRequired != NULL ) {
      if ( pt->tag & MG_REQ )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}

int MMG3D_Set_prism(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2,
                    MMG5_int v3, MMG5_int v4, MMG5_int v5, MMG5_int ref, MMG5_int pos) {
  MMG5_pPrism pp;

  if ( !mesh->nprism ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of prisms with the",
            __func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting elements in mesh\n");
    return 0;
  }

  if ( pos > mesh->nprism ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new prism at position %" MMG5_PRId ".",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of prism: %" MMG5_PRId "\n",mesh->nprism);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the prism.\n");
    return 0;
  }

  pp = &mesh->prism[pos];
  pp->v[0] = v0;
  pp->v[1] = v1;
  pp->v[2] = v2;
  pp->v[3] = v3;
  pp->v[4] = v4;
  pp->v[5] = v5;
  pp->ref  = ref;

  mesh->point[pp->v[0]].tag &= ~MG_NUL;
  mesh->point[pp->v[1]].tag &= ~MG_NUL;
  mesh->point[pp->v[2]].tag &= ~MG_NUL;
  mesh->point[pp->v[3]].tag &= ~MG_NUL;
  mesh->point[pp->v[4]].tag &= ~MG_NUL;
  mesh->point[pp->v[5]].tag &= ~MG_NUL;


  return 1;
}

int MMG3D_Get_prism(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* v3,
                    MMG5_int* v4, MMG5_int* v5, MMG5_int* ref, int* isRequired) {
  static MMG5_int npri = 0;

  if ( npri == mesh->nprism ) {
    npri = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of prisms.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_prism function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of prisms: %" MMG5_PRId "\n ",mesh->nprism);
    }
  }

  ++npri;

  if ( npri > mesh->nprism ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get prism.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_prism function");
    fprintf(stderr," can not exceed the number of prism: %" MMG5_PRId "\n ",mesh->nprism);
    return 0;
  }

  *v0  = mesh->prism[npri].v[0];
  *v1  = mesh->prism[npri].v[1];
  *v2  = mesh->prism[npri].v[2];
  *v3  = mesh->prism[npri].v[3];
  *v4  = mesh->prism[npri].v[4];
  *v5  = mesh->prism[npri].v[5];

  if ( ref != NULL ) {
    *ref = mesh->prism[npri].ref;
  }

  if ( isRequired != NULL ) {
    if ( mesh->prism[npri].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}

int  MMG3D_Set_prisms(MMG5_pMesh mesh, MMG5_int *prisms, MMG5_int *refs) {
  MMG5_pPrism pp;
  MMG5_int    j,i;

  for (i=1;i<=mesh->nprism;i++)
  {
    j = (i-1)*6;
    pp = &mesh->prism[i];
    pp->v[0]  = prisms[j];
    pp->v[1]  = prisms[j+1];
    pp->v[2]  = prisms[j+2];
    pp->v[3]  = prisms[j+3];
    pp->v[4]  = prisms[j+4];
    pp->v[5]  = prisms[j+5];

    if ( refs != NULL )
      pp->ref   = refs[i-1];

    mesh->point[pp->v[0]].tag &= ~MG_NUL;
    mesh->point[pp->v[1]].tag &= ~MG_NUL;
    mesh->point[pp->v[2]].tag &= ~MG_NUL;
    mesh->point[pp->v[3]].tag &= ~MG_NUL;
    mesh->point[pp->v[4]].tag &= ~MG_NUL;
    mesh->point[pp->v[5]].tag &= ~MG_NUL;

  }

  return 1;
}

int  MMG3D_Get_prisms(MMG5_pMesh mesh, MMG5_int *prisms, MMG5_int *refs, int * areRequired) {
  MMG5_pPrism pp;
  MMG5_int    j,i;

  for (i=1;i<=mesh->nprism;i++)
  {
    j = (i-1)*6;
    pp = &mesh->prism[i];
    prisms[j]   = pp->v[0];
    prisms[j+2] = pp->v[1];
    prisms[j+1] = pp->v[2];
    prisms[j+3] = pp->v[3];
    prisms[j+4] = pp->v[4];
    prisms[j+5] = pp->v[5];

    if ( refs!=NULL )
      refs[i-1]  = pp->ref ;
    if ( areRequired != NULL ) {
      if ( pp->tag & MG_REQ )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}



int MMG3D_Set_triangle(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2, MMG5_int ref,MMG5_int pos) {

  if ( !mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of triangles"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting triangles in mesh\n");
    return 0;
  }

  if ( pos > mesh->ntmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new triangle.\n",
            __func__);
    fprintf(stderr,"    max number of triangle: %" MMG5_PRId "\n",mesh->ntmax);
    MMG5_INCREASE_MEM_MESSAGE();
    return 0;
  }

  if ( pos > mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new triangle at"
            " position %" MMG5_PRId ".",__func__,pos);
    fprintf(stderr," Overflow of the given number of triangles: %" MMG5_PRId "\n",mesh->nt);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the triangle.\n");
    return 0;
  }

  mesh->tria[pos].v[0] = v0;
  mesh->tria[pos].v[1] = v1;
  mesh->tria[pos].v[2] = v2;
  mesh->tria[pos].ref  = ref;

  return 1;
}

int MMG3D_Get_triangle(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* ref
                       ,int* isRequired) {
  MMG5_pTria  ptt;

  if ( mesh->nti == mesh->nt ) {
    mesh->nti = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of"
              " triangles.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_triangle function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of triangles: %" MMG5_PRId "\n ",mesh->nt);
    }
  }

  mesh->nti++;

  if ( mesh->nti > mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get triangle.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_triangle function");
    fprintf(stderr," can not exceed the number of triangles: %" MMG5_PRId "\n ",mesh->nt);
    return 0;
  }

  ptt = &mesh->tria[mesh->nti];
  *v0  = ptt->v[0];
  *v1  = ptt->v[1];
  *v2  = ptt->v[2];
  if ( ref != NULL )
    *ref = ptt->ref;

  if ( isRequired != NULL ) {
    if ( (ptt->tag[0] & MG_REQ) && (ptt->tag[1] & MG_REQ) &&
         (ptt->tag[2] & MG_REQ) )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}
int  MMG3D_Set_triangles(MMG5_pMesh mesh, MMG5_int *tria, MMG5_int *refs) {

  MMG5_pTria ptt;
  MMG5_int   i, j;

  for (i=1;i<=mesh->nt;i++)
  {
    j = (i-1)*3;
    ptt = &mesh->tria[i];
    ptt->v[0] = tria[j]  ;
    ptt->v[1] = tria[j+1];
    ptt->v[2] = tria[j+2];
    if ( refs != NULL )
      ptt->ref  = refs[i-1];
  }
  return 1;
}

int  MMG3D_Get_triangles(MMG5_pMesh mesh, MMG5_int *tria, MMG5_int *refs, int *areRequired) {
  MMG5_pTria ptt;
  MMG5_int   i, j;

  for (i=1;i<=mesh->nt;i++)
  {
    j = (i-1)*3;
    ptt = &mesh->tria[i];
    tria[j]   = ptt->v[0];
    tria[j+1] = ptt->v[1];
    tria[j+2] = ptt->v[2];

    if ( refs!=NULL )
      refs[i-1]  = ptt->ref ;
    if ( areRequired != NULL ) {
      if ( (ptt->tag[0] & MG_REQ) && (ptt->tag[1] & MG_REQ) &&
           (ptt->tag[2] & MG_REQ) )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}

int MMG3D_Set_quadrilateral(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2, MMG5_int v3,
                            MMG5_int ref,MMG5_int pos) {

  if ( !mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of quadrilaterals"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting quadrilaterals in mesh\n");
    return 0;
  }

  if ( pos > mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new quadrilateral"
            " at position %" MMG5_PRId ".",__func__,pos);
    fprintf(stderr," Overflow of the given number of quadrilaterals: %" MMG5_PRId "\n",mesh->nquad);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the quadrilateral.\n");
    return 0;
  }

  mesh->quadra[pos].v[0] = v0;
  mesh->quadra[pos].v[1] = v1;
  mesh->quadra[pos].v[2] = v2;
  mesh->quadra[pos].v[3] = v3;
  mesh->quadra[pos].ref  = ref;

  return 1;
}

int MMG3D_Get_quadrilateral(MMG5_pMesh mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* v2, MMG5_int* v3,
                            MMG5_int* ref,int* isRequired) {
  MMG5_pQuad       pq;
  static MMG5_int  nqi = 0;

  if ( nqi == mesh->nquad ) {
    nqi = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter"
              " of quadrilaterals.\n",__func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_quadrilateral function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of quadrilaterals: %" MMG5_PRId "\n ",mesh->nquad);
    }
  }

  nqi++;

  if ( nqi > mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get quadrilateral.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_quadrilateral function");
    fprintf(stderr," can not exceed the number of quadrilaterals: %" MMG5_PRId "\n ",mesh->nquad);
    return 0;
  }

  pq = &mesh->quadra[nqi];
  *v0  = pq->v[0];
  *v1  = pq->v[1];
  *v2  = pq->v[2];
  *v3  = pq->v[3];
  if ( ref != NULL )
    *ref = pq->ref;

  if ( isRequired != NULL ) {
    if ( (pq->tag[0] & MG_REQ) && (pq->tag[1] & MG_REQ) &&
         (pq->tag[2] & MG_REQ) && (pq->tag[3] & MG_REQ))
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}

int  MMG3D_Set_quadrilaterals(MMG5_pMesh mesh, MMG5_int *quads, MMG5_int *refs) {
  MMG5_pQuad  pq;
  MMG5_int    j,i;

  for (i=1;i<=mesh->nquad;i++)
  {
    j = (i-1)*4;
    pq = &mesh->quadra[i];
    pq->v[0] = quads[j]  ;
    pq->v[1] = quads[j+1];
    pq->v[2] = quads[j+2];
    pq->v[3] = quads[j+3];
    if ( refs != NULL )
      pq->ref  = refs[i-1];
  }
  return 1;
}

int  MMG3D_Get_quadrilaterals(MMG5_pMesh mesh, MMG5_int *quads, MMG5_int *refs, int *areRequired) {
  MMG5_pQuad  pq;
   MMG5_int   j,i;

  for (i=1;i<=mesh->nquad;i++)
  {
    j = (i-1)*4;
    pq = &mesh->quadra[i];
    quads[j]   = pq->v[0];
    quads[j+1] = pq->v[1];
    quads[j+2] = pq->v[2];
    quads[j+3] = pq->v[3];

    if ( refs!=NULL )
      refs[i-1]  = pq->ref ;
    if ( areRequired != NULL ) {
      if ( (pq->tag[0] & MG_REQ) && (pq->tag[1] & MG_REQ) &&
           (pq->tag[2] & MG_REQ) && (pq->tag[3] & MG_REQ) )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }
  return 1;
}

int MMG3D_Set_edge(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int ref, MMG5_int pos) {

  if ( !mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of edges with"
            " the",__func__);
    fprintf(stderr," MMG3D_Set_meshSize function before setting edges in mesh\n");
    return 0;
  }
  if ( pos > mesh->namax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to allocate a new edge.\n",
            __func__);
    fprintf(stderr,"    max number of edge: %" MMG5_PRId "\n",mesh->namax);
    MMG5_INCREASE_MEM_MESSAGE();
    return 0;
  }
  if ( pos > mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new edge at position %" MMG5_PRId ".",
            __func__,pos);
    fprintf(stderr," Overflow of the given number of edges: %" MMG5_PRId "\n",mesh->na);
    fprintf(stderr,"\n  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the edge.\n");
    return 0;
  }

  mesh->edge[pos].a = v0;
  mesh->edge[pos].b = v1;
  mesh->edge[pos].ref  = ref;
  mesh->edge[pos].tag |= MG_REF;

  return 1;
}

int MMG3D_Get_edge(MMG5_pMesh mesh, MMG5_int* e0, MMG5_int* e1, MMG5_int* ref
                   ,int* isRidge, int* isRequired) {

  if ( mesh->nai == mesh->na ) {
    mesh->nai = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of edges.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_edge function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of edges: %" MMG5_PRId "\n ",mesh->na);
    }
  }

  mesh->nai++;

  if ( mesh->nai > mesh->na ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get edge.\n",__func__);
    fprintf(stderr,"    The number of call of MMG3D_Get_edge function");
    fprintf(stderr," can not exceed the number of edges: %" MMG5_PRId "\n ",mesh->na);
    return 0;
  }

  *e0  = mesh->edge[mesh->nai].a;
  *e1  = mesh->edge[mesh->nai].b;
  if ( ref!=NULL )
    *ref = mesh->edge[mesh->nai].ref;

  if ( isRidge != NULL ) {
    if ( mesh->edge[mesh->nai].tag & MG_GEO )
      *isRidge = 1;
    else
      *isRidge = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->edge[mesh->nai].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return 1;
}

int MMG3D_Set_edges(MMG5_pMesh mesh, MMG5_int *edges, MMG5_int *refs) {
  MMG5_int i,j;

  for (i=1;i<=mesh->na;i++)
  {
    j = (i-1)*2;

    mesh->edge[i].a    = edges[j];
    mesh->edge[i].b    = edges[j+1];
    if ( refs != NULL )
      mesh->edge[i].ref  = refs[i-1];
    mesh->edge[i].tag |= MG_REF;
  }

  return 1;
}

int MMG3D_Get_edges(MMG5_pMesh mesh, MMG5_int* edges,MMG5_int *refs,int* areRidges,int* areRequired) {
  MMG5_int i,j;

  for (i=1;i<=mesh->na;i++)
  {
    j = (i-1)*2;
    edges[j]   = mesh->edge[i].a;
    edges[j+1] = mesh->edge[i].b;

    if ( refs!=NULL )
      refs[i-1] = mesh->edge[i].ref;

    if ( areRidges != NULL ) {
      if ( mesh->edge[i].tag & MG_GEO )
        areRidges[i-1] = 1;
      else
        areRidges[i-1] = 0;
    }

    if ( areRequired != NULL ) {
      if ( mesh->edge[i].tag & MG_REQ )
        areRequired[i-1] = 1;
      else
        areRequired[i-1] = 0;
    }
  }

  return 1;
}

int MMG3D_Set_corner(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_CRN;
  return 1;
}

int MMG3D_Unset_corner(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag &= ~MG_CRN;
  return 1;
}

int MMG3D_Set_requiredVertex(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_REQ;
  mesh->point[k].tag &= ~MG_NUL;
  return 1;
}

int MMG3D_Unset_requiredVertex(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag &= ~MG_REQ;
  return 1;
}

int MMG3D_Set_requiredTetrahedron(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->ne );
  mesh->tetra[k].tag |= MG_REQ;
  return 1;
}

int MMG3D_Unset_requiredTetrahedron(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->ne );
  mesh->tetra[k].tag &= ~MG_REQ;
  return 1;
}

int MMG3D_Set_requiredTetrahedra(MMG5_pMesh mesh, MMG5_int *reqIdx, MMG5_int nreq) {
  MMG5_int k;

  for ( k=0; k<nreq; ++k ){
    mesh->tetra[reqIdx[k]].tag |= MG_REQ;
  }

  return 1;
}

int MMG3D_Unset_requiredTetrahedra(MMG5_pMesh mesh, MMG5_int *reqIdx, MMG5_int nreq) {
  MMG5_int k;

  for ( k=0; k<nreq; ++k ){
    mesh->tetra[reqIdx[k]].tag &= ~MG_REQ;
  }

  return 1;
}

int MMG3D_Set_requiredTriangle(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->nt );
  mesh->tria[k].tag[0] |= MG_REQ;
  mesh->tria[k].tag[1] |= MG_REQ;
  mesh->tria[k].tag[2] |= MG_REQ;
  return 1;
}

int MMG3D_Unset_requiredTriangle(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->nt );
  mesh->tria[k].tag[0] &= ~MG_REQ;
  mesh->tria[k].tag[1] &= ~MG_REQ;
  mesh->tria[k].tag[2] &= ~MG_REQ;
  return 1;
}

int MMG3D_Set_requiredTriangles(MMG5_pMesh mesh, MMG5_int* reqIdx, MMG5_int nreq) {
  MMG5_int k;

  for ( k=0; k<nreq; ++k ){
    mesh->tria[reqIdx[k]].tag[0] |= MG_REQ;
    mesh->tria[reqIdx[k]].tag[1] |= MG_REQ;
    mesh->tria[reqIdx[k]].tag[2] |= MG_REQ;
  }
  return 1;
}

int MMG3D_Unset_requiredTriangles(MMG5_pMesh mesh, MMG5_int* reqIdx, MMG5_int nreq) {
  MMG5_int k;

  for ( k=0; k<nreq; ++k ){
    mesh->tria[reqIdx[k]].tag[0] &= ~MG_REQ;
    mesh->tria[reqIdx[k]].tag[1] &= ~MG_REQ;
    mesh->tria[reqIdx[k]].tag[2] &= ~MG_REQ;
  }
  return 1;
}

int MMG3D_Set_parallelTriangle(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->nt );
  mesh->tria[k].tag[0] |= MG_PARBDY;
  mesh->tria[k].tag[1] |= MG_PARBDY;
  mesh->tria[k].tag[2] |= MG_PARBDY;
  return 1;
}

int MMG3D_Unset_parallelTriangle(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->nt );
  mesh->tria[k].tag[0] &= ~MG_PARBDY;
  mesh->tria[k].tag[1] &= ~MG_PARBDY;
  mesh->tria[k].tag[2] &= ~MG_PARBDY;
  return 1;
}

int MMG3D_Set_parallelTriangles(MMG5_pMesh mesh, MMG5_int* parIdx, MMG5_int npar) {
  MMG5_int k;

  for ( k=0; k<npar; ++k ){
    mesh->tria[parIdx[k]].tag[0] |= MG_PARBDY;
    mesh->tria[parIdx[k]].tag[1] |= MG_PARBDY;
    mesh->tria[parIdx[k]].tag[2] |= MG_PARBDY;
  }
  return 1;
}

int MMG3D_Unset_parallelTriangles(MMG5_pMesh mesh, MMG5_int* parIdx, MMG5_int npar) {
  MMG5_int k;

  for ( k=0; k<npar; ++k ){
    mesh->tria[parIdx[k]].tag[0] &= ~MG_PARBDY;
    mesh->tria[parIdx[k]].tag[1] &= ~MG_PARBDY;
    mesh->tria[parIdx[k]].tag[2] &= ~MG_PARBDY;
  }
  return 1;
}

int MMG3D_Set_ridge(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag |= MG_GEO;
  return 1;
}

int MMG3D_Unset_ridge(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag &= ~MG_GEO;
  return 1;
}

int MMG3D_Set_requiredEdge(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag |= MG_REQ;
  return 1;
}

int MMG3D_Unset_requiredEdge(MMG5_pMesh mesh, MMG5_int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag &= ~MG_REQ;
  return 1;
}

int MMG3D_Set_normalAtVertex(MMG5_pMesh mesh, MMG5_int k, double n0, double n1, double n2) {

  assert ( k <= mesh->np );
  mesh->point[k].n[0] = n0;
  mesh->point[k].n[1] = n1;
  mesh->point[k].n[2] = n2;

  ++mesh->nc1;

  return 1;
}

int MMG3D_Get_normalAtVertex(MMG5_pMesh mesh, MMG5_int k, double *n0, double *n1, double *n2) {

  assert ( k <= mesh->np );
  (*n0) = mesh->point[k].n[0];
  (*n1) = mesh->point[k].n[1];
  (*n2) = mesh->point[k].n[2];

  return 1;
}

double MMG3D_Get_tetrahedronQuality(MMG5_pMesh mesh,MMG5_pSol met, MMG5_int k) {
  double      qual = 0.;
  MMG5_pTetra pt;

  if ( k < 1 || k > mesh->ne ) {
    fprintf(stderr,"\n  ## Error: %s: unable to access to tetra %" MMG5_PRId ".\n",
            __func__,k);
    fprintf(stderr,"     Tetra numbering goes from 1 to %" MMG5_PRId "\n",mesh->ne);
    return 0.;
  }
  pt = &mesh->tetra[k];
  assert ( MG_EOK(pt) );

  if ( (!met) || (!met->m) || met->size==1 ) {
    if ( mesh->info.optimLES) {
      /* Skewness */
      qual =  MMG3D_ALPHAD * MMG3D_caltetLES_iso(mesh,met,pt);
    } else {
      /* iso quality */
      qual =  MMG3D_ALPHAD * MMG5_caltet_iso(mesh,NULL,pt);
    }
  }
  else if ( !mesh->info.metRidTyp ) {
    qual =  MMG3D_ALPHAD * MMG5_caltet33_ani(mesh,met,pt);
  }
  else {
    qual = MMG3D_ALPHAD * MMG5_caltet_ani(mesh,met,pt);
  }

  return qual;
}

int MMG3D_Set_scalarSol(MMG5_pSol met, double s, MMG5_int pos) {

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }
  if ( pos < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return 0;
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    max number of solutions: %" MMG5_PRId "\n",met->npmax);
    return 0;
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution at"
            " position %" MMG5_PRId ".",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %" MMG5_PRId "\n",met->np);
    fprintf(stderr,"\n  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return 0;
  }

  met->m[pos] = s;
  return 1;
}


int MMG3D_Get_scalarSol(MMG5_pSol met, double* s) {

  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_scalarSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",__func__);
    fprintf(stderr,"     The number of call of MMG3D_Get_scalarSol function");
    fprintf(stderr," can not exceed the number of points: %" MMG5_PRId "\n ",met->np);
    return 0;
  }

  *s  = met->m[met->npi];

  return 1;
}

int MMG3D_Set_scalarSols(MMG5_pSol met, double *s ) {
  MMG5_int k;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of solution"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }

  for ( k=0; k<met->np; ++k )
    met->m[k+1] = s[k];

  return 1;
}

int MMG3D_Get_scalarSols(MMG5_pSol met, double* s) {
  MMG5_int k;

  for ( k=0; k<met->np; ++k )
    s[k]  = met->m[k+1];

  return 1;
}

int MMG3D_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz, MMG5_int pos) {

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of solution"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }
  if ( pos < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return 0;
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    max number of solutions: %" MMG5_PRId "\n",met->npmax);
    return 0;
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution at"
            " position %" MMG5_PRId ".",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %" MMG5_PRId "\n",met->np);
    fprintf(stderr,"\n  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return 0;
  }

  met->m[3*pos]   = vx;
  met->m[3*pos+1] = vy;
  met->m[3*pos+2] = vz;

  return 1;
}


int MMG3D_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz) {

  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_vectorSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",__func__);
    fprintf(stderr,"     The number of call of MMG3D_Get_vectorSol function");
    fprintf(stderr," can not exceed the number of points: %" MMG5_PRId "\n ",met->np);
    return 0;
  }

  *vx  = met->m[3*met->npi];
  *vy  = met->m[3*met->npi+1];
  *vz  = met->m[3*met->npi+2];

  return 1;
}

int MMG3D_Set_vectorSols(MMG5_pSol met, double *sols) {
  double   *m;
  MMG5_int k,j;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of solution"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j+3];
    m[0] = sols[j];
    m[1] = sols[j+1];
    m[2] = sols[j+2];
  }

  return 1;
}

int MMG3D_Get_vectorSols(MMG5_pSol met, double* sols) {
  double   *m;
  MMG5_int k, j;

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j+3];
    sols[j]   = m[0];
    sols[j+1] = m[1];
    sols[j+2] = m[2];
  }

  return 1;
}

int MMG3D_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                        double m22,double m23, double m33, MMG5_int pos) {

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of solution"
            " with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }
  if ( pos < 1 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",
            __func__);
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return 0;
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new solution.\n",__func__);
    fprintf(stderr,"    max number of solutions: %" MMG5_PRId "\n",met->npmax);
    return 0;
  }

  if ( pos > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: attempt to set new solution at "
            "position %" MMG5_PRId ".",__func__,pos);
    fprintf(stderr," Overflow of the given number of solutions: %" MMG5_PRId "\n",met->np);
    fprintf(stderr,"\n  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return 0;
  }

  met->m[6*pos]   = m11;
  met->m[6*pos+1] = m12;
  met->m[6*pos+2] = m13;
  met->m[6*pos+3] = m22;
  met->m[6*pos+4] = m23;
  met->m[6*pos+5] = m33;

  return 1;
}


int MMG3D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                        double *m22,double *m23, double *m33) {

  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the MMG3D_Get_tensorSol function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %" MMG5_PRId "\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get solution.\n",__func__);
    fprintf(stderr,"     The number of call of MMG3D_Get_tensorSol function");
    fprintf(stderr," can not exceed the number of points: %" MMG5_PRId "\n ",met->np);
    return 0;
  }

  *m11 = met->m[6*met->npi];
  *m12 = met->m[6*met->npi+1];
  *m13 = met->m[6*met->npi+2];
  *m22 = met->m[6*met->npi+3];
  *m23 = met->m[6*met->npi+4];
  *m33 = met->m[6*met->npi+5];

  return 1;
}

int MMG3D_Set_tensorSols(MMG5_pSol met, double *sols) {
  double   *m;
  MMG5_int k,j;

  if ( !met->np ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of"
            " solution with the",__func__);
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return 0;
  }

  for ( k=0; k<met->np; ++k ) {
    j = 6*k;
    m = &met->m[j+6];

    m[0] = sols[j];
    m[1] = sols[j+1];
    m[2] = sols[j+2];
    m[3] = sols[j+3];
    m[4] = sols[j+4];
    m[5] = sols[j+5];
  }
  return 1;
}

int MMG3D_Get_tensorSols(MMG5_pSol met, double *sols) {
  double   *m;
  MMG5_int k,j;

  for ( k=0; k<met->np; ++k ) {
    j = 6*k;
    m = &met->m[j+6];

    sols[j]   = m[0];
    sols[j+1] = m[1];
    sols[j+2] = m[2];
    sols[j+3] = m[3];
    sols[j+4] = m[4];
    sols[j+5] = m[5];
  }

  return 1;
}

int  MMG3D_Set_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double* s,MMG5_int pos) {
  MMG5_pSol psl;

  /* Warning: users give indices from 1 to nsols */
  psl = sol + (i-1);

  switch ( psl->type ) {
  case MMG5_Scalar:
    return MMG3D_Set_scalarSol(psl,s[0],pos);
    break;

  case MMG5_Vector:
    MMG3D_Set_vectorSol(psl,s[0],s[1],s[2],pos);
    break;

  case MMG5_Tensor:
    MMG3D_Set_tensorSol(psl,s[0],s[1],s[2],s[3],s[4],s[5],pos);
    break;

  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected type of solution: %s.\n",
            __func__,MMG5_Get_typeName(psl->type));
    return 0;
  }
  return 1;
}

int  MMG3D_Get_ithSol_inSolsAtVertices(MMG5_pSol sol,int i, double *s,MMG5_int pos) {
  MMG5_pSol psl;

  /* Warning: users give indices from 1 to nsols */
  psl = sol + (i-1);

  psl->npi = pos-1;

  switch ( psl->type ) {
  case MMG5_Scalar:
    return MMG3D_Get_scalarSol(psl,&s[0]);
    break;

  case MMG5_Vector:
    MMG3D_Get_vectorSol(psl,&s[0],&s[1],&s[2]);
    break;

  case MMG5_Tensor:
    MMG3D_Get_tensorSol(psl,&s[0],&s[1],&s[2],&s[3],&s[4],&s[5]);
    break;

  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected type of solution: %s\n",
            __func__,MMG5_Get_typeName(psl->type));
    return 0;
  }

  return 1;
}

int  MMG3D_Set_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double *s) {
  MMG5_pSol psl;

  /* Warning: users give indices from 1 to nsols */
  psl = sol + (i-1);

  switch ( psl->type ) {
  case MMG5_Scalar:
    return MMG3D_Set_scalarSols(psl,s);
    break;

  case MMG5_Vector:
    MMG3D_Set_vectorSols(psl,s);
    break;

  case MMG5_Tensor:
    MMG3D_Set_tensorSols(psl,s);
    break;

  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected type of solution: %s.\n",
            __func__,MMG5_Get_typeName(psl->type));
    return 0;
  }

  return 1;
}

int  MMG3D_Get_ithSols_inSolsAtVertices(MMG5_pSol sol,int i, double *s) {
  MMG5_pSol psl;

  /* Warning: users give indices from 1 to nsols */
  psl = sol + (i-1);

  switch ( psl->type ) {
  case MMG5_Scalar:
    return MMG3D_Get_scalarSols(psl,s);
    break;

  case MMG5_Vector:
    MMG3D_Get_vectorSols(psl,s);
    break;

  case MMG5_Tensor:
    MMG3D_Get_tensorSols(psl,s);
    break;

  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected type of solution: %s\n",
            __func__,MMG5_Get_typeName(psl->type));
    return 0;
  }

  return 1;
}

void MMG3D_Set_handGivenMesh(MMG5_pMesh mesh) {
  MMG5_int k, aux;

  /* Possibly switch 2 vertices number so that each tet is positively oriented */
  for (k=1; k<=mesh->ne; k++) {
    if ( MMG5_orvol(mesh->point,mesh->tetra[k].v) < 0.0 ) {
      /* mesh->xt temporary used to count reoriented tetra */
      mesh->xt++;
      aux = mesh->tetra[k].v[2];
      mesh->tetra[k].v[2] = mesh->tetra[k].v[3];
      mesh->tetra[k].v[3] = aux;
    }
  }
  return;
}

int MMG3D_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( (mesh->npi != mesh->np) || (mesh->nei != mesh->ne) ) {
    fprintf(stderr,"\n  ## Error: %s: if you don't use the MMG3D_loadMesh"
            " function,",__func__);
    fprintf(stderr," you must call the MMG3D_Set_meshSize function to have a");
    fprintf(stderr," valid mesh.\n");
    fprintf(stderr," Missing datas.\n");
    return 0;
  }

  if ( met->npi != met->np ) {
    fprintf(stderr,"\n  ## Error: %s: if you don't use the MMG3D_loadSol"
            " function,",__func__);
    fprintf(stderr," you must call the MMG3D_Set_solSize function to have a");
    fprintf(stderr," valid solution.\n");
    fprintf(stderr," Missing datas.\n");
    return 0;
  }

  /*  Check mesh data */
  if ( mesh->info.ddebug ) {
    if ( (!mesh->np) || (!mesh->point) ||
         (!mesh->ne) || (!mesh->tetra) ) {
      fprintf(stderr,"  ** MISSING DATA.\n");
      fprintf(stderr," Check that your mesh contains points and tetrahedra.\n");
      fprintf(stderr," Exit program.\n");
      return 0;
    }
  }

  if ( mesh->dim != 3 ) {
    fprintf(stderr,"  ** 3 DIMENSIONAL MESH NEEDED. Exit program.\n");
    return 0;
  }
  if ( met->dim != 3 ) {
    fprintf(stderr,"  ** WRONG DIMENSION FOR METRIC. Exit program.\n");
    return 0;
  }
  if ( !mesh->ver )  mesh->ver = 2;
  if ( !met ->ver )  met ->ver = 2;

  return 1;
}

static inline
int MMG3D_skipIso(MMG5_pMesh mesh) {
  MMG5_pTria  ptt,ptt1;
  MMG5_pEdge  pa,pa1;
  MMG5_int    k;

  if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
    fprintf(stderr,"\n  ## Warning: %s: skip of all entites with %"MMG5_PRId
            " reference.\n",__func__,mesh->info.isoref);

  /* Skip triangles with mesh->info.isoref refs */
  k = 1;
  do {
    ptt = &mesh->tria[k];
    if ( MMG5_abs(ptt->ref) != mesh->info.isoref ) continue;
    /* here ptt is the first tri of mesh->tria that we want to delete */
    do {
      ptt1 = &mesh->tria[mesh->nti];
    }
    while( (MMG5_abs(ptt1->ref) == mesh->info.isoref) && (k <= --mesh->nti) );

    if ( MMG5_abs(ptt1->ref) != mesh->info.isoref )
      /* ptt1 is the last tri of mesh->tria that we want to keep */
      memcpy(ptt,ptt1,sizeof(MMG5_Tria));
  } while( ++k <= mesh->nti );

  if ( mesh->nti < mesh->nt ) {
    if( !mesh->nti )
      MMG5_DEL_MEM(mesh,mesh->tria);
    else {
      MMG5_ADD_MEM(mesh,mesh->nti-mesh->nt,"triangles",return 0);
      MMG5_SAFE_RECALLOC(mesh->tria,mesh->nt+1,(mesh->nti+1),MMG5_Tria,
                          "triangles",return 0);
    }
    mesh->nt = mesh->nti;
  }

  /* Skip edges with mesh->info.isoref refs */
  if ( mesh->na ) {
    k = 1;
    do {
      pa = &mesh->edge[k];
      if ( MMG5_abs(pa->ref) != mesh->info.isoref ) {
        pa->ref = MMG5_abs(pa->ref);
        continue;
      }
      /* here pa is the first edge of mesh->edge that we want to delete */
      do {
        pa1 = &mesh->edge[mesh->nai];
      }
      while( (MMG5_abs(pa1->ref) == mesh->info.isoref) && (k <= --mesh->nai) );

      if ( MMG5_abs(pa1->ref) != mesh->info.isoref ) {
        /* pa1 is the last edge of mesh->edge that we want to keep */
        memcpy(pa,pa1,sizeof(MMG5_Edge));
        pa1->ref = MMG5_abs(pa1->ref);
      }
    } while( ++k <= mesh->nai );

    if ( mesh->nai < mesh->na ) {
      if( !mesh->nai )
        MMG5_DEL_MEM(mesh,mesh->edge);
      else {
        MMG5_ADD_MEM(mesh,mesh->nai-mesh->na,"Edges",return 0);
        MMG5_SAFE_RECALLOC(mesh->edge,mesh->na+1,(mesh->nai+1),MMG5_Edge,
                            "edges",return 0);
      }
      mesh->na = mesh->nai;
    }
  }
  assert ( mesh->nai == mesh->na );

  /* delete tetrahedra references */
  for (k=1; k<=mesh->ne; k++) {
    mesh->tetra[k].ref = 0;
  }
  return 1;
}

int MMG3D_Add_tetrahedron(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, MMG5_int v2, MMG5_int v3, MMG5_int ref) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  double      vol;
  int         j,ip;
  MMG5_int    aux,vv[4],iel;

  vv[0] = v0;
  vv[1] = v1;
  vv[2] = v2;
  vv[3] = v3;

  for ( j=0; j<4; ++j ) {
    if ( vv[j] > mesh->np ) {
      fprintf(stderr,"\n  ## Error: %s: vertex %" MMG5_PRId " doesn't exist in the mesh.\n",
              __func__,vv[j]);
      fprintf(stderr,"    Use the MMG3D_Add_vertex function to add it.\n");
      return 0;
    }
  }

  iel = MMG3D_newElt(mesh);
  if ( !iel ) {
    MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                        fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                " a new element.\n",__func__);
                        MMG5_INCREASE_MEM_MESSAGE();
                        fprintf(stderr,"  Exit program.\n");
                        return 0);
  }

  pt = &mesh->tetra[iel];
  pt->v[0] = v0;
  pt->v[1] = v1;
  pt->v[2] = v2;
  pt->v[3] = v3;
  pt->ref  = MMG5_abs(ref);

  mesh->point[pt->v[0]].tag &= ~MG_NUL;
  mesh->point[pt->v[1]].tag &= ~MG_NUL;
  mesh->point[pt->v[2]].tag &= ~MG_NUL;
  mesh->point[pt->v[3]].tag &= ~MG_NUL;

  vol = MMG5_orvol(mesh->point,pt->v);
  if ( fabs(vol) <= MMG5_EPSD2 ) {
    fprintf(stderr,"\n  ## Error: %s: tetrahedron %" MMG5_PRId ": null volume.\n",
            __func__,iel);
    for ( ip=0; ip<4; ip++ ) {
      ppt = &mesh->point[pt->v[ip]];
      for ( j=0; j<3; j++ ) {
        if ( fabs(ppt->c[j])>0. ) {
          fprintf(stderr," Check that you don't have a sliver tetrahedron.\n");
          return -iel;
        }
      }
    }
    fprintf(stderr,"  All vertices have zero coordinates.");
    fprintf(stderr," Check that you have set the vertices before the tetrahedra.\n");
    return -iel;
  }
  else if ( vol < 0.0 ) {
    /* Possibly switch 2 vertices number so that each tet is positively oriented */
    aux = pt->v[2];
    pt->v[2] = pt->v[3];
    pt->v[3] = aux;
    /* mesh->xt temporary used to count reoriented tetra */
    mesh->xt++;

    return -iel;
  }

  return iel;
}

MMG5_int MMG3D_Add_vertex(MMG5_pMesh mesh,double c0,double c1,double c2,MMG5_int ref) {
  double      c[3];
  MMG5_int    ip,klink;

  c[0] = c0;
  c[1] = c1;
  c[2] = c2;

  ip = MMG3D_newPt(mesh,c,0,1);
  if ( !ip ) {
    MMG5_TAB_RECALLOC(mesh,mesh->point,mesh->npmax,mesh->gap,MMG5_Point,
                       "larger point table",
                       fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                               " a new point\n",__func__);
                       MMG5_INCREASE_MEM_MESSAGE();return 0);

    mesh->npnil = mesh->np+1;
    for (klink=mesh->npnil; klink<mesh->npmax-1; klink++)
      mesh->point[klink].tmp  = klink+1;

    /* We try again to add the point */
    ip = MMG3D_newPt(mesh,c,0,1);
    if ( !ip ) {
      fprintf(stderr,"\n  ## Error: %s: unable to allocate"
              " a new point\n",__func__);
      MMG5_INCREASE_MEM_MESSAGE();
      return 0;
    }
  }
  return ip;
}

int MMG3D_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam,MMG5_int val){
  int k;

  switch ( iparam ) {
    /* Integer parameters */
  case MMG3D_IPARAM_verbose :
    mesh->info.imprim   = val;
    break;
  case MMG3D_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf(stderr,"\n  ## Warning: %s: maximal memory authorized must be"
              " strictly positive.\n",__func__);
      fprintf(stderr,"  Reset to default value.\n");
    }
    else
      mesh->info.mem      = val;
    if ( !MMG3D_memOption(mesh) )  return 0;
    break;
#ifndef MMG_PATTERN
  case MMG3D_IPARAM_octree :
    mesh->info.PROctree   = val;
    break;
#endif
  case MMG3D_IPARAM_debug :
    mesh->info.ddebug   = val;
    break;
  case MMG3D_IPARAM_angle :
    /* free table that may contains old ridges */
    if ( mesh->htab.geom )
      MMG5_DEL_MEM(mesh,mesh->htab.geom);
    if ( mesh->xpoint )
      MMG5_DEL_MEM(mesh,mesh->xpoint);
    if ( mesh->xtetra )
      MMG5_DEL_MEM(mesh,mesh->xtetra);
    if ( !val )
      mesh->info.dhd    = -1.;
    else {
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: angle detection parameter set"
                " to default value\n",__func__);
      mesh->info.dhd    = MMG5_ANGEDG;
    }
    break;
  case MMG3D_IPARAM_nofem :
    mesh->info.setfem = (val==1)? 0 : 1;
    break;
  case MMG3D_IPARAM_opnbdy :
    mesh->info.opnbdy = val;
    break;
  case MMG3D_IPARAM_iso :
    mesh->info.iso      = val;
    if ( mesh->info.iso )
      if ( mesh->nt && !MMG3D_skipIso(mesh) )
        return 0;
    break;
  case MMG3D_IPARAM_isoref :
    mesh->info.isoref   = val;
    break;
  case MMG3D_IPARAM_isosurf :
    mesh->info.isosurf = val;
    break;
  case MMG3D_IPARAM_lag :
#ifdef USE_ELAS
    if ( val < 0 || val > 2 )
      return 0;
    mesh->info.lag = val;
    /* No connectivity changes unless lag >= 2 */
    if ( val < 2 ) {
      if ( !MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_noinsert,1) )
        return 0;
    }
#else
    fprintf(stderr,"\n  ## Error: %s"
            " \"lagrangian motion\" option unavailable (-lag):\n"
            " set the USE_ELAS CMake's flag to ON when compiling the mmg3d"
            " library to enable this feature.\n",__func__);
    return 0;
#endif
    break;
  case MMG3D_IPARAM_numsubdomain :
    mesh->info.nsd = val;
    break;
  case MMG3D_IPARAM_optim :
    mesh->info.optim = val;
    break;
  case MMG3D_IPARAM_optimLES :
    mesh->info.optimLES = val;
    break;
  case MMG3D_IPARAM_noinsert :
    mesh->info.noinsert = val;
    break;
  case MMG3D_IPARAM_noswap :
    mesh->info.noswap   = val;
    break;
  case MMG3D_IPARAM_nomove :
    mesh->info.nomove   = val;
    break;
  case MMG3D_IPARAM_nosurf :
    mesh->info.nosurf   = val;
    break;
  case MMG3D_IPARAM_nreg :
    mesh->info.nreg     = val;
    break;
  case MMG3D_IPARAM_xreg :
    mesh->info.xreg     = val;
    break;
  case MMG3D_IPARAM_nosizreq :
    mesh->info.nosizreq = val;
    break;
  case MMG3D_IPARAM_numberOfLocalParam :
    if ( mesh->info.par ) {
      MMG5_DEL_MEM(mesh,mesh->info.par);
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: new local parameter values\n",__func__);
    }
    mesh->info.npar   = val;
    mesh->info.npari  = 0;
    mesh->info.parTyp = 0;

    MMG5_ADD_MEM(mesh,mesh->info.npar*sizeof(MMG5_Par),"parameters",
                  printf("  Exit program.\n");
                  return 0);
    MMG5_SAFE_CALLOC(mesh->info.par,mesh->info.npar,MMG5_Par,return 0);

    MMG5_int inival;
    inival = MMG5_INTMAX;

    for (k=0; k<mesh->info.npar; k++) {
      mesh->info.par[k].elt   = MMG5_Noentity;
      mesh->info.par[k].ref   = inival;
      mesh->info.par[k].hausd = mesh->info.hausd;
      mesh->info.par[k].hmin  = mesh->info.hmin;
      mesh->info.par[k].hmax  = mesh->info.hmax;
    }

    break;
  case MMG3D_IPARAM_numberOfLSBaseReferences :
    mesh->info.nbr = val;
    MMG5_ADD_MEM(mesh,mesh->info.nbr*sizeof(int),"References",
                 printf("  Exit program.\n");
                 return 0);
    MMG5_SAFE_CALLOC(mesh->info.br,mesh->info.nbr,MMG5_int,return 0);

    for (k=0; k<mesh->info.nbr; k++)
      mesh->info.br[k] = 0;

    break;

  case MMG3D_IPARAM_numberOfMat :
    if ( mesh->info.mat ) {
      MMG5_DEL_MEM(mesh,mesh->info.mat);
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stderr,"\n  ## Warning: %s: new multi materials values\n",__func__);
    }
    mesh->info.nmat   = val;
    mesh->info.nmati  = 0;

    MMG5_ADD_MEM(mesh,(mesh->info.nmat)*sizeof(MMG5_Mat),"multi material",
                 printf("  Exit program.\n");
                 return 0);
    MMG5_SAFE_CALLOC(mesh->info.mat,mesh->info.nmat,MMG5_Mat,return 0);
    for (k=0; k<mesh->info.nmat; k++) {
      mesh->info.mat[k].ref   = 0;
    }

    break;

#ifdef USE_SCOTCH
  case MMG3D_IPARAM_renum :
    mesh->info.renum    = val;
    break;
#endif
  case MMG3D_IPARAM_anisosize :
    mesh->info.ani = val;
    break;
  default :
    fprintf(stderr,"\n  ## Error: %s: unknown type of parameter\n",__func__);
    return 0;
  }
  /* other options */

  return 1;
}

int MMG3D_Get_iparameter(MMG5_pMesh mesh, MMG5_int iparam) {

  switch ( iparam ) {
    /* Integer parameters */
  case MMG3D_IPARAM_verbose :
    return  mesh->info.imprim;
    break;
  case MMG3D_IPARAM_mem :
    return  mesh->info.mem;
    break;
#ifndef MMG_PATTERN
  case MMG3D_IPARAM_octree :
    return  mesh->info.PROctree;
    break;
#endif
  case MMG3D_IPARAM_debug :
    return  mesh->info.ddebug;
    break;
  case MMG3D_IPARAM_angle :
    if ( mesh->info.dhd <= 0. ) {
      return  0;
    }
    else {
      return  1;
    }
    break;
  case MMG3D_IPARAM_iso :
    return  mesh->info.iso;
    break;
  case MMG3D_IPARAM_lag :
    return  mesh->info.lag;
    break;
  case MMG3D_IPARAM_noinsert :
    return  mesh->info.noinsert;
    break;
  case MMG3D_IPARAM_noswap :
    return  mesh->info.noswap;
    break;
  case MMG3D_IPARAM_nomove :
    return  mesh->info.nomove;
    break;
  case MMG3D_IPARAM_nosurf :
    return  mesh->info.nosurf;
    break;
  case MMG3D_IPARAM_nreg :
    return mesh->info.nreg;
    break;
  case MMG3D_IPARAM_xreg :
    return mesh->info.xreg;
    break;
  case MMG3D_IPARAM_numberOfLocalParam :
    return  mesh->info.npar;
    break;
  case MMG3D_IPARAM_numberOfMat :
    return  mesh->info.nmat;
    break;
#ifdef USE_SCOTCH
  case MMG3D_IPARAM_renum :
    return  mesh->info.renum;
    break;
#endif
  default :
    fprintf(stderr,"\n  ## Error: %s: unknown type of parameter\n",__func__);
    return 0;
  }
}

int MMG3D_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val){

  switch ( dparam ) {
    /* double parameters */
  case MMG3D_DPARAM_angleDetection :
    mesh->info.dhd = val;
    mesh->info.dhd = MG_MAX(0.0, MG_MIN(180.0,mesh->info.dhd));
    mesh->info.dhd = cos(mesh->info.dhd*M_PI/180.0);
    break;
  case MMG3D_DPARAM_hmin :
    mesh->info.sethmin  = 1;
    mesh->info.hmin     = val;
    if ( mesh->info.sethmax && ( mesh->info.hmin >=  mesh->info.hmax ) ) {
      fprintf(stderr,"\n  ## Warning: hmin value must be strictly lower than hmax one"
              " (hmin = %lf  hmax = %lf ).\n",mesh->info.hmin, mesh->info.hmax);
    }

    break;
  case MMG3D_DPARAM_hmax :
    mesh->info.sethmax  = 1;
    mesh->info.hmax     = val;
    if ( mesh->info.sethmin && ( mesh->info.hmin >=  mesh->info.hmax ) ) {
      fprintf(stderr,"\n  ## Warning: hmin value must be strictly lower than hmax one"
              " (hmin = %lf  hmax = %lf ).\n",mesh->info.hmin, mesh->info.hmax);
    }

    break;
  case MMG3D_DPARAM_hsiz :
    mesh->info.hsiz     = val;
    break;
  case MMG3D_DPARAM_hgrad :
    mesh->info.hgrad    = val;
    if ( mesh->info.hgrad <= 0.0 )
      mesh->info.hgrad = -1.0;
    else
      mesh->info.hgrad = log(mesh->info.hgrad);
    break;
  case MMG3D_DPARAM_hgradreq :
    mesh->info.hgradreq    = val;
    if ( mesh->info.hgradreq <= 0.0 )
      mesh->info.hgradreq = -1.0;
    else
      mesh->info.hgradreq = log(mesh->info.hgradreq);
    break;
  case MMG3D_DPARAM_hausd :
    if ( val <=0 ) {
      fprintf(stderr,"\n  ## Error: %s: hausdorff number must be strictly"
              " positive.\n",__func__);
      return 0;
    }
    else
      mesh->info.hausd    = val;
    break;
  case MMG3D_DPARAM_ls :
    mesh->info.ls       = val;
    break;
  case MMG3D_DPARAM_xreg :
    if (val < 0.0 || val > 1.0) {
      fprintf(stderr,"\n  ## Error: %s: Coordinate regularization parameter must be comprised between 0 and 1.\n",__func__);
    }
    else
      mesh->info.lxreg    = val;
    break;
  case MMG3D_DPARAM_rmc :
    if ( !val ) {
      /* Default value */
      mesh->info.rmc      = MMG3D_VOLFRAC;
    }
    else {
      /* User customized value */
      mesh->info.rmc      = val;
    }
    break;
  default :
    fprintf(stderr,"\n  ## Error: %s: unknown type of parameter\n", __func__);
    return 0;
  }
  return 1;
}

int MMG3D_Set_localParameter(MMG5_pMesh mesh,MMG5_pSol sol, int typ, MMG5_int ref,
                             double hmin,double hmax,double hausd){
  MMG5_pPar par;
  int k;

  if ( !mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of local"
            " parameters",__func__);
    fprintf(stderr," with the MMG3D_Set_iparameters function before setting");
    fprintf(stderr," values in local parameters structure. \n");
    return 0;
  }
  if ( mesh->info.npari >= mesh->info.npar ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new local parameter.\n",
            __func__);
    fprintf(stderr,"    max number of local parameters: %d\n",mesh->info.npar);
    return 0;
  }
  if ( typ != MMG5_Triangle && typ != MMG5_Tetrahedron ) {
    fprintf(stderr,"\n  ## Warning: %s: you must apply your local parameters",
            __func__);
    fprintf(stderr," on triangles (MMG5_Triangle or %d) or tetrahedron"
            " (MMG5_Tetrahedron or %d).\n",MMG5_Triangle,MMG5_Tetrahedron);
    fprintf(stderr,"\n  ## Unknown type of entity: ignored.\n");
    return 0;
  }
  if ( ref < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative references are not allowed.\n",
            __func__);
    return 0;
  }

  if ( hmin <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative hmin value is not allowed.\n",
            __func__);
    return 0;
  }
  if ( hmax <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative hmax value is not allowed.\n",
            __func__);
    return 0;
  }
  if ( hausd <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative hausd value is not allowed.\n",
            __func__);
    return 0;
  }

  for (k=0; k<mesh->info.npari; k++) {
    par = &mesh->info.par[k];

    if ( par->elt == typ && par->ref == ref ) {
      par->hausd = hausd;
      par->hmin  = hmin;
      par->hmax  = hmax;
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
        fprintf(stderr,"\n  ## Warning: %s: new parameters (hausd, hmin and hmax)",
                __func__);
        fprintf(stderr," for entities of type %d and of ref %" MMG5_PRId "\n",typ,ref);
      }
      return 1;
    }
  }

  mesh->info.par[mesh->info.npari].elt   = typ;
  mesh->info.par[mesh->info.npari].ref   = ref;
  mesh->info.par[mesh->info.npari].hmin  = hmin;
  mesh->info.par[mesh->info.npari].hmax  = hmax;
  mesh->info.par[mesh->info.npari].hausd = hausd;

  switch ( typ )
  {
  case ( MMG5_Triangle ):
    mesh->info.parTyp |= MG_Tria;
    break;
  case ( MMG5_Tetrahedron ):
    mesh->info.parTyp |= MG_Tetra;
    break;
  default:
    fprintf(stderr,"\n  ## Error: %s: unexpected entity type: %s.\n",
            __func__,MMG5_Get_entitiesName(typ));
    return 0;
  }

  mesh->info.npari++;

  return 1;
}

int MMG3D_Set_multiMat(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int ref,int split,MMG5_int rmin,MMG5_int rplus) {
  return MMG5_Set_multiMat(mesh,sol,ref,split,rmin,rplus);
}

int MMG3D_Set_lsBaseReference(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int br){
  return MMG5_Set_lsBaseReference(mesh,sol,br);
}


int MMG3D_Free_allSols(MMG5_pMesh mesh,MMG5_pSol *sol) {

  return MMG5_Free_allSols(mesh,sol);
}

int MMG3D_Free_all(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = MMG3D_Free_all_var(argptr);

  va_end(argptr);

  return ier;
}

int MMG3D_Free_structures(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = MMG3D_Free_structures_var(argptr);

  va_end(argptr);

  return ier;
}

int MMG3D_Free_names(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = MMG3D_Free_names_var(argptr);

  va_end(argptr);

  return ier;
}
