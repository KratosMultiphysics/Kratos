/* =============================================================================
**  This file is part of the mmg software package for the triangular
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
 * \file mmgs/API_functions_s.c
 * \brief C API functions definitions for MMGS library.
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmgs/libmmgs.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * C API for MMGS library.
 *
 */

#include "mmgs.h"

void MMGS_Init_mesh(const int starter,...) {
  va_list argptr;

  va_start(argptr, starter);

  _MMGS_Init_mesh_var(argptr);

  va_end(argptr);

  return;
}

void MMGS_Init_fileNames(MMG5_pMesh mesh,MMG5_pSol sol
  ) {

  MMG5_Init_fileNames(mesh,sol);
  return;
}

int MMGS_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin) {

  return(MMG5_Set_inputMeshName(mesh,meshin));
}

int MMGS_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin) {
  return(MMG5_Set_inputSolName(mesh,sol,solin));
}

int MMGS_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout) {

  return(MMG5_Set_outputMeshName(mesh,meshout));
}

int MMGS_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout) {
  return(MMG5_Set_outputSolName(mesh,sol,solout));
}
void MMGS_Init_parameters(MMG5_pMesh mesh) {

  /* Init common parameters for mmgs and mmgs. */
  _MMG5_Init_parameters(mesh);

  mesh->info.renum    = 0;   /* [0/1], Turn off/on the renumbering using SCOTCH; */

}

int MMGS_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && sol->m )
    fprintf(stdout,"  ## Warning: new solution\n");

  if ( typEntity != MMG5_Vertex ) {
    fprintf(stderr,"  ## Error: MMGS5 need a solution imposed on vertices\n");
    return(0);
  }
  if ( typSol == MMG5_Scalar ) {
    sol->size = 1;
  }
  else if ( typSol == MMG5_Vector ) {
    sol->size = 3;
  }
  else if ( typSol == MMG5_Tensor ) {
    sol->size = 6;
  }
  else {
    fprintf(stderr,"  ## Error: type of solution not yet implemented\n");
    return(0);
  }

  sol->dim = 3;
  if ( np ) {
    sol->np  = np;
    sol->npi = np;
    if ( sol->m )
      _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));

    sol->npmax = mesh->npmax;
    _MMG5_ADD_MEM(mesh,(sol->size*(sol->npmax+1))*sizeof(double),"initial solution",
                  fprintf(stderr,"  Exit program.\n");
                  exit(EXIT_FAILURE));
    _MMG5_SAFE_CALLOC(sol->m,(sol->size*(sol->npmax+1)),double);
  }
  return(1);
}

int MMGS_Set_meshSize(MMG5_pMesh mesh, int np, int nt, int na) {
  int k;

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->tria || mesh->edge) )
    fprintf(stdout,"  ## Warning: new mesh\n");

  mesh->np  = np;
  mesh->nt  = nt;
  mesh->na  = na;
  mesh->npi = mesh->np;
  mesh->nti = mesh->nt;
  mesh->nai = mesh->na;

  if ( !np || !nt ) {
    fprintf(stderr,"  ** MISSING DATA:\n");
    fprintf(stderr,"     Your mesh must contains at least points and triangles.\n");
    return(0);
  }

  if ( mesh->point )
    _MMG5_DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(MMG5_Point));
  if ( mesh->tria )
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->nt+1)*sizeof(MMG5_Tria));
  if ( mesh->edge )
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->na+1)*sizeof(MMG5_Edge));

  /*tester si -m defini : renvoie 0 si pas ok et met la taille min dans info.mem */
  if( mesh->info.mem > 0) {
    if ( mesh->npmax < mesh->np || mesh->ntmax < mesh->nt) {
      _MMGS_memOption(mesh);

      if ( mesh->npmax < mesh->np || mesh->ntmax < mesh->nt) {
        fprintf(stderr,"not enough memory: np : %d %d nt : %d %d \n"
                ,mesh->npmax,mesh->np, mesh->ntmax,mesh->nt);
        return(0);
      }
    } else if(mesh->info.mem < 39) {
      fprintf(stderr,"not enough memory  %d\n",mesh->info.mem);
      return(0);
    }
  } else {
    mesh->npmax = MG_MAX(1.5*mesh->np,_MMG5_NPMAX);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,_MMG5_NTMAX);

  }
  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                fprintf(stderr,"  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point);

  _MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",return(0));
  _MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria);


  mesh->namax = mesh->na;
  if ( mesh->na ) {
    _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"initial edges",return(0));
    _MMG5_SAFE_CALLOC(mesh->edge,(mesh->na+1),MMG5_Edge);
  }

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->nt + 1;
  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    mesh->point[k].tmp  = k+1;
  }
  for (k=mesh->nenil; k<mesh->ntmax-1; k++) {
    mesh->tria[k].v[2] = k+1;
  }

  /* stats */
  if ( abs(mesh->info.imprim) > 6 ) {
    fprintf(stdout,"     NUMBER OF VERTICES     %8d\n",mesh->np);
    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES        %8d\n",mesh->na);
    }
    fprintf(stdout,"     NUMBER OF TRIANGLES    %8d\n",mesh->nt);
  }
  return(1);
}

int MMGS_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np, int* typSol) {

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

  sol->npi = 0;

  return(1);
}

int MMGS_Get_meshSize(MMG5_pMesh mesh, int* np, int* nt, int* na) {

  if ( np != NULL )
    *np = mesh->np;
  if ( nt != NULL )
    *nt = mesh->nt;
  if ( na != NULL )
    *na = mesh->na;

  return(1);
}

int MMGS_Set_vertex(MMG5_pMesh mesh, double c0, double c1, double c2, int ref, int pos) {

  if ( !mesh->np ) {
    fprintf(stderr,"  ## Error: you must set the number of points with the");
    fprintf(stderr," MMGS_Set_meshSize function before setting vertices in mesh\n");
    return(0);
  }

  if ( pos > mesh->npmax ) {
    fprintf(stderr,"  ## Error: unable to allocate a new point.\n");
    fprintf(stderr,"    max number of points: %d\n",mesh->npmax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->np ) {
    fprintf(stderr,"  ## Error: attempt to set new vertex at position %d.",pos);
    fprintf(stderr," Overflow of the given number of vertices: %d\n",mesh->np);
    fprintf(stderr,"  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the vertex.\n");
    return(0);
  }

  mesh->point[pos].c[0] = c0;
  mesh->point[pos].c[1] = c1;
  mesh->point[pos].c[2] = c2;
  mesh->point[pos].ref  = ref;
  mesh->point[pos].tag  = MG_NUL;
  mesh->point[pos].flag = 0;
  mesh->point[pos].tmp = 0;

  return(1);
}

int  MMGS_Set_vertices(MMG5_pMesh mesh, double *vertices,int *refs) {

  MMG5_pPoint ppt;
  int i,j;

  /*coordinates vertices*/
  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*3;
    ppt->c[0]  = vertices[j];
    ppt->c[1]  = vertices[j+1];
    ppt->c[2]  = vertices[j+2];

    ppt->tag  = MG_NUL;
    ppt->flag = 0;
    ppt->tmp = 0;

    if ( refs != NULL )
      ppt->ref   = refs[i-1];
  }

  return 1;
}


int MMGS_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, double* c2, int* ref,
                    int* isCorner, int* isRequired) {

 if ( mesh->npi == mesh->np ) {
   mesh->npi = 0;
   if ( mesh->info.ddebug ) {
    fprintf(stdout,"  ## Warning: reset the internal counter of points.\n");
    fprintf(stdout,"     You must pass here exactly one time (the first time ");
    fprintf(stdout,"you call the MMGS_Get_vertex function).\n");
    fprintf(stdout,"     If not, the number of call of this function");
    fprintf(stdout," exceed the number of points: %d\n ",mesh->np);
   }
 }

  mesh->npi++;

  if ( mesh->npi > mesh->np ) {
    fprintf(stderr,"  ## Error: unable to get point.\n");
    fprintf(stderr,"     The number of call of MMGS_Get_vertex function");
    fprintf(stderr," can not exceed the number of points: %d\n ",mesh->np);
    return(0);
  }

  *c0  = mesh->point[mesh->npi].c[0];
  *c1  = mesh->point[mesh->npi].c[1];
  *c2  = mesh->point[mesh->npi].c[2];
  if ( ref != NULL )
    *ref = mesh->point[mesh->npi].ref;

  if ( isCorner != NULL ) {
    if ( mesh->point[mesh->npi].tag & MG_CRN )
      *isCorner = 1;
    else
      *isCorner = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->point[mesh->npi].tag & MG_REQ )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

int  MMGS_Get_vertices(MMG5_pMesh mesh, double* vertices, int* refs,
                        int* areCorners, int* areRequired) {
  MMG5_pPoint ppt;
  int i,j;

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

int MMGS_Set_triangle(MMG5_pMesh mesh, int v0, int v1, int v2, int ref,int pos) {

  if ( !mesh->nt ) {
    fprintf(stderr,"  ## Error: You must set the number of triangles with the");
    fprintf(stderr," MMGS_Set_meshSize function before setting triangles in mesh\n");
    return(0);
  }

  if ( pos > mesh->ntmax ) {
    fprintf(stderr,"  ## Error: unable to allocate a new triangle.\n");
    fprintf(stderr,"    max number of triangle: %d\n",mesh->ntmax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->nt ) {
    fprintf(stderr,"  ## Error: attempt to set new triangle at position %d.",pos);
    fprintf(stderr," Overflow of the given number of triangles: %d\n",mesh->nt);
    fprintf(stderr,"  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the triangle.\n");
    return(0);
  }

  mesh->tria[pos].v[0] = v0;
  mesh->tria[pos].v[1] = v1;
  mesh->tria[pos].v[2] = v2;
  mesh->tria[pos].ref  = ref;

  mesh->point[v0].tag &= ~MG_NUL;
  mesh->point[v1].tag &= ~MG_NUL;
  mesh->point[v2].tag &= ~MG_NUL;

  return(1);
}

int MMGS_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref
                      ,int* isRequired) {
  MMG5_pTria  ptt;

  if ( mesh->nti == mesh->nt ) {
    mesh->nti = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of triangles.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMGS_Get_triangle function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of triangles: %d\n ",mesh->nt);
    }
  }

  mesh->nti++;

  if ( mesh->nti > mesh->nt ) {
    fprintf(stderr,"  ## Error: unable to get triangle.\n");
    fprintf(stderr,"    The number of call of MMGS_Get_triangle function");
    fprintf(stderr," can not exceed the number of triangles: %d\n ",mesh->nt);
    return(0);
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

  return(1);
}

int  MMGS_Set_triangles(MMG5_pMesh mesh, int *tria, int *refs) {
  MMG5_pTria ptt;
  int         i, j;

   for (i=1;i<=mesh->nt;i++)
   {
      j = (i-1)*3;
      ptt = &mesh->tria[i];
      ptt->v[0] = tria[j]  ;
      ptt->v[1] = tria[j+2];
      ptt->v[2] = tria[j+1];

      mesh->point[ptt->v[0]].tag &= ~MG_NUL;
      mesh->point[ptt->v[1]].tag &= ~MG_NUL;
      mesh->point[ptt->v[2]].tag &= ~MG_NUL;

      if ( refs != NULL )
        ptt->ref  = refs[i-1];
   }
   return 1;
}

int  MMGS_Get_triangles(MMG5_pMesh mesh, int *tria, int *refs, int *areRequired) {
  MMG5_pTria ptt;
  int         i, j;

   for (i=1;i<=mesh->nt;i++)
   {
      j = (i-1)*3;
      ptt = &mesh->tria[i];
      tria[j]   = ptt->v[0];
      tria[j+2] = ptt->v[1];
      tria[j+1] = ptt->v[2];

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


int MMGS_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref, int pos) {

  if ( !mesh->na ) {
    fprintf(stderr,"  ## Error: You must set the number of edges with the");
    fprintf(stderr," MMGS_Set_meshSize function before setting edges in mesh\n");
    return(0);
  }
  if ( pos > mesh->namax ) {
    fprintf(stderr,"  ## Error: unable to allocate a new edge.\n");
    fprintf(stderr,"    max number of edge: %d\n",mesh->namax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }
  if ( pos > mesh->na ) {
    fprintf(stderr,"  ## Error: attempt to set new edge at position %d.",pos);
    fprintf(stderr," Overflow of the given number of edges: %d\n",mesh->na);
    fprintf(stderr,"  ## Check the mesh size, its compactness or the position");
    fprintf(stderr," of the edge.\n");
    return(0);
  }

  mesh->edge[pos].a = v0;
  mesh->edge[pos].b = v1;
  mesh->edge[pos].ref  = ref;
  mesh->edge[pos].tag |= MG_REF;

  return(1);
}

int MMGS_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref
                  ,int* isRidge, int* isRequired) {

  if ( mesh->nai == mesh->na ) {
    mesh->nai = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of edges.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMGS_Get_edge function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of edges: %d\n ",mesh->na);
    }
  }

  mesh->nai++;

  if ( mesh->nai > mesh->na ) {
    fprintf(stderr,"  ## Error: unable to get edge.\n");
    fprintf(stderr,"    The number of call of MMGS_Get_edge function");
    fprintf(stderr," can not exceed the number of edges: %d\n ",mesh->na);
    return(0);
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

  return(1);
}

int MMGS_Set_corner(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_CRN;
  return(1);
}

int MMGS_Set_requiredVertex(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->np );
  mesh->point[k].tag |= MG_REQ;
  return(1);
}

int MMGS_Set_requiredTriangle(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->nt );
  mesh->tria[k].tag[0] |= MG_REQ;
  mesh->tria[k].tag[1] |= MG_REQ;
  mesh->tria[k].tag[2] |= MG_REQ;
  return(1);
}

int MMGS_Set_ridge(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag |= MG_GEO;
  return(1);
}

int MMGS_Set_requiredEdge(MMG5_pMesh mesh, int k) {
  assert ( k <= mesh->na );
  mesh->edge[k].tag |= MG_REQ;
  return(1);
}

int MMGS_Set_normalAtVertex(MMG5_pMesh mesh, int k, double n0, double n1, double n2) {

  assert ( k <= mesh->np );
  mesh->point[k].n[0] = n0;
  mesh->point[k].n[1] = n1;
  mesh->point[k].n[2] = n2;

  ++mesh->nc1;

  return(1);
}

int MMGS_Get_normalAtVertex(MMG5_pMesh mesh, int k, double *n0, double *n1, double *n2) {

  assert ( k <= mesh->np );
  (*n0) = mesh->point[k].n[0];
  (*n1) = mesh->point[k].n[1];
  (*n2) = mesh->point[k].n[2];

  return(1);
}

int MMGS_Set_scalarSol(MMG5_pSol met, double s, int pos) {

  if ( !met->np ) {
    fprintf(stderr,"  ## Error: You must set the number of solution with the");
    fprintf(stderr," MMGS_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }
  if ( pos < 1 ) {
    fprintf(stderr,"  ## Error: unable to set a new solution.\n");
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return(0);
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"  ## Error: unable to set a new solution.\n");
    fprintf(stderr,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stderr,"  ## Error: attempt to set new solution at position %d.",pos);
    fprintf(stderr," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stderr,"  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return(0);
  }

  met->m[pos] = s;
  return(1);
}

int MMGS_Get_scalarSol(MMG5_pSol met, double* s) {
  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of points.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMGS_Get_scalarSol function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of points: %d\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"  ## Error: unable to get solution.\n");
    fprintf(stderr,"     The number of call of MMGS_Get_scalarSol function");
    fprintf(stderr," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  *s  = met->m[met->npi];

  return(1);
}

int MMGS_Set_scalarSols(MMG5_pSol met, double *s ) {
  int k;

  if ( !met->np ) {
    fprintf(stderr,"  ## Error: You must set the number of solution with the");
    fprintf(stderr," MMGS_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }

  for ( k=0; k<met->np; ++k )
    met->m[k+1] = s[k];

  return(1);
}

int MMGS_Get_scalarSols(MMG5_pSol met, double* s) {
  int k;

  for ( k=0; k<met->np; ++k )
    s[k]  = met->m[k+1];

  return(1);
}


int MMGS_Set_vectorSol(MMG5_pSol met, double vx,double vy, double vz, int pos) {

  if ( !met->np ) {
    fprintf(stderr,"  ## Error: You must set the number of solution with the");
    fprintf(stderr," MMGS_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }
  if ( pos < 1 ) {
    fprintf(stderr,"  ## Error: unable to set a new solution.\n");
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return(0);
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"  ## Error: unable to set a new solution.\n");
    fprintf(stderr,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stderr,"  ## Error: attempt to set new solution at position %d.",pos);
    fprintf(stderr," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stderr,"  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return(0);
  }

  met->m[3*pos]   = vx;
  met->m[3*pos+1] = vy;
  met->m[3*pos+2] = vz;

  return(1);
}

int MMGS_Get_vectorSol(MMG5_pSol met, double* vx, double* vy, double* vz) {
  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of points.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMGS_Get_vectorSol function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of points: %d\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"  ## Error: unable to get solution.\n");
    fprintf(stderr,"     The number of call of MMGS_Get_vectorSol function");
    fprintf(stderr," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  *vx  = met->m[3*met->npi];
  *vy  = met->m[3*met->npi+1];
  *vz  = met->m[3*met->npi+2];

  return(1);
}


int MMGS_Set_vectorSols(MMG5_pSol met, double *sols) {
  double *m;
  int k,j;

  if ( !met->np ) {
    fprintf(stderr,"  ## Error: You must set the number of solution with the");
    fprintf(stderr," MMGS_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j+3];
    m[0] = sols[j];
    m[1] = sols[j+1];
    m[2] = sols[j+2];
  }

  return(1);
}

int MMGS_Get_vectorSols(MMG5_pSol met, double* sols) {
  double *m;
  int k, j;

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j+3];
    sols[j]   = m[0];
    sols[j+1] = m[1];
    sols[j+2] = m[2];
  }

  return(1);
}

int MMGS_Set_tensorSol(MMG5_pSol met, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos) {

  if ( !met->np ) {
    fprintf(stderr,"  ## Error: You must set the number of solution with the");
    fprintf(stderr," MMGS_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }
  if ( pos < 1 ) {
    fprintf(stderr,"  ## Error: unable to set a new solution.\n");
    fprintf(stderr,"    Minimal index of the solution position must be 1.\n");
    return(0);
  }
  if ( pos >= met->npmax ) {
    fprintf(stderr,"  ## Error: unable to set a new solution.\n");
    fprintf(stderr,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stderr,"  ## Error: attempt to set new solution at position %d.",pos);
    fprintf(stderr," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stderr,"  ## Check the solution size, its compactness or the position");
    fprintf(stderr," of the solution.\n");
    return(0);
  }

  met->m[6*pos]   = m11;
  met->m[6*pos+1] = m12;
  met->m[6*pos+2] = m13;
  met->m[6*pos+3] = m22;
  met->m[6*pos+4] = m23;
  met->m[6*pos+5] = m33;

  return(1);
}

int MMGS_Get_tensorSol(MMG5_pSol met, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33) {

  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of points.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMGS_Get_tensorSol function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of points: %d\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stderr,"  ## Error: unable to get solution.\n");
    fprintf(stderr,"     The number of call of MMGS_Get_tensorSol function");
    fprintf(stderr," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  *m11 = met->m[6*met->npi];
  *m12 = met->m[6*met->npi+1];
  *m13 = met->m[6*met->npi+2];
  *m22 = met->m[6*met->npi+3];
  *m13 = met->m[6*met->npi+4];
  *m33 = met->m[6*met->npi+5];

  return(1);
}


int MMGS_Set_tensorSols(MMG5_pSol met, double *sols) {
  double *m;
  int k,j;

  if ( !met->np ) {
    fprintf(stderr,"  ## Error: You must set the number of solution with the");
    fprintf(stderr," MMGS_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
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
  return(1);
}

int MMGS_Get_tensorSols(MMG5_pSol met, double *sols) {
  double *m;
  int k,j;

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

  return(1);
}


int MMGS_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( (mesh->npi != mesh->np) || (mesh->nti != mesh->nt) ) {
    fprintf(stderr,"  ## Error: if you don't use the MMGS_loadMesh function,");
    fprintf(stderr," you must call the MMGS_Set_meshSize function to have a");
    fprintf(stderr," valid mesh.\n");
    fprintf(stderr," Missing datas.\n");
    return(0);
  }

  if ( met->npi != met->np ) {
    fprintf(stderr,"  ## Error: if you don't use the MMGS_loadSol function,");
    fprintf(stderr," you must call the MMGS_Set_solSize function to have a");
    fprintf(stderr," valid solution.\n");
    fprintf(stderr," Missing datas.\n");
    return(0);
  }

  /*  Check mesh data */
  if ( mesh->info.ddebug ) {
    if ( (!mesh->np) || (!mesh->point) ||
         (!mesh->nt) || (!mesh->tria) ) {
      fprintf(stderr,"  ** MISSING DATA.\n");
      fprintf(stderr," Check that your mesh contains points and triangles.\n");
      fprintf(stderr," Exit program.\n");
      return(0);
    }
  }

  if ( mesh->dim != 3 ) {
    fprintf(stderr,"  ** 3 DIMENSIONAL MESH NEEDED. Exit program.\n");
    return(0);
  }
  if ( met->dim != 3 ) {
    fprintf(stderr,"  ** WRONG DIMENSION FOR METRIC. Exit program.\n");
    return(0);
  }
  if ( !mesh->ver )  mesh->ver = 2;
  if ( !met ->ver )  met ->ver = 2;

  return(1);
}

int MMGS_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam, int val){
  int k;

  switch ( iparam ) {
    /* Integer parameters */
  case MMGS_IPARAM_verbose :
    mesh->info.imprim   = val;
    break;
  case MMGS_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf(stdout,"  ## Warning: maximal memory authorized must be strictly positive.\n");
      fprintf(stdout,"  Reset to default value.\n");
    }
    else
      mesh->info.mem      = val;
    _MMGS_memOption(mesh);
    if(mesh->np && (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt)) {
      return(0);
    } else if(mesh->info.mem < 39)
      return(0);
    break;
  case MMGS_IPARAM_debug :
    mesh->info.ddebug   = val;
    break;
  case MMGS_IPARAM_angle :
    if ( mesh->xpoint )
      _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
    if ( !val )
      mesh->info.dhd    = -1.;
    else {
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stdout,"  ## Warning: angle detection parameter set to default value\n");
      mesh->info.dhd    = _MMG5_ANGEDG;
    }
    break;
  case MMGS_IPARAM_iso :
    if ( !mesh->info.iso )
      mesh->info.iso      = val;
    break;
  case MMGS_IPARAM_keepRef :
    if ( val )
      mesh->info.iso      = 2;
    break;
  case MMGS_IPARAM_noinsert :
    mesh->info.noinsert = val;
    break;
  case MMGS_IPARAM_noswap :
    mesh->info.noswap   = val;
    break;
  case MMGS_IPARAM_nomove :
    mesh->info.nomove   = val;
    break;
  case MMGS_IPARAM_nreg :
    mesh->info.nreg     = val;
    break;
  case MMGS_IPARAM_numberOfLocalParam :
    if ( mesh->info.par ) {
      _MMG5_DEL_MEM(mesh,mesh->info.par,mesh->info.npar*sizeof(MMG5_Par));
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stdout,"  ## Warning: new local parameter values\n");
    }
    mesh->info.npar   = val;
    mesh->info.npari  = 0;
    mesh->info.parTyp = 0;

    _MMG5_ADD_MEM(mesh,mesh->info.npar*sizeof(MMG5_Par),"parameters",
                  fprintf(stderr,"  Exit program.\n");
                  exit(EXIT_FAILURE));
    _MMG5_SAFE_CALLOC(mesh->info.par,mesh->info.npar,MMG5_Par);

    for (k=0; k<mesh->info.npar; k++) {
      mesh->info.par[k].elt   = MMG5_Noentity;
      mesh->info.par[k].ref   = INT_MAX;
      mesh->info.par[k].hausd = mesh->info.hausd;
      mesh->info.par[k].hmin  = mesh->info.hmin;
      mesh->info.par[k].hmax  = mesh->info.hmax;
    }

    break;
#ifdef USE_SCOTCH
  case MMGS_IPARAM_renum :
    mesh->info.renum    = val;
    break;
#endif
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    return(0);
  }
  /* other options */
  mesh->info.fem      = 0;
  return(1);
}

int MMGS_Get_iparameter(MMG5_pMesh mesh, int iparam) {

  switch ( iparam ) {
    /* Integer parameters */
  case MMGS_IPARAM_verbose :
    return ( mesh->info.imprim );
    break;
  case MMGS_IPARAM_mem :
    return ( mesh->info.mem );
    break;
  case MMGS_IPARAM_debug :
    return ( mesh->info.ddebug );
    break;
  case MMGS_IPARAM_angle :
    if ( mesh->info.dhd <= 0. ) {
      return ( 0 );
    }
    else {
      return ( 1 );
    }
    break;
  case MMGS_IPARAM_noinsert :
    return ( mesh->info.noinsert );
    break;
  case MMGS_IPARAM_noswap :
    return ( mesh->info.noswap );
    break;
  case MMGS_IPARAM_nomove :
    return ( mesh->info.nomove );
    break;
  case MMGS_IPARAM_nreg :
    return ( mesh->info.nreg );
    break;
  case MMGS_IPARAM_numberOfLocalParam :
    return ( mesh->info.npar );
    break;
#ifdef USE_SCOTCH
  case MMGS_IPARAM_renum :
    return ( mesh->info.renum );
    break;
#endif
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    exit(EXIT_FAILURE);
  }
}

int MMGS_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val){

  switch ( dparam ) {
    /* double parameters */
  case MMGS_DPARAM_angleDetection :
    mesh->info.dhd = val;
    mesh->info.dhd = MG_MAX(0.0, MG_MIN(180.0,mesh->info.dhd));
    mesh->info.dhd = cos(mesh->info.dhd*M_PI/180.0);
    break;
  case MMGS_DPARAM_hmin :
    mesh->info.hmin     = val;
    break;
  case MMGS_DPARAM_hmax :
    mesh->info.hmax     = val;
    break;
  case MMGS_DPARAM_hgrad :
    mesh->info.hgrad    = val;
    if ( mesh->info.hgrad < 0.0 )
      mesh->info.hgrad = -1.0;
    else
      mesh->info.hgrad = log(mesh->info.hgrad);
    break;
  case MMGS_DPARAM_hausd :
    if ( val <=0 ) {
      fprintf(stderr,"  ## Error: hausdorff number must be strictly positive.\n");
      return(0);
    }
    else
      mesh->info.hausd    = val;
    break;
  case MMGS_DPARAM_ls :
    mesh->info.ls       = val;
    break;
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    return(0);
  }
  return(1);
}

int MMGS_Set_localParameter(MMG5_pMesh mesh,MMG5_pSol sol, int typ, int ref,
                            double hmin,double hmax,double hausd){
  MMG5_pPar par;
  int k;

  if ( !mesh->info.npar ) {
    fprintf(stderr,"  ## Error: You must set the number of local parameters");
    fprintf(stderr," with the MMGS_Set_iparameters function before setting");
    fprintf(stderr," values in local parameters structure. \n");
    return(0);
  }
  if ( mesh->info.npari > mesh->info.npar ) {
    fprintf(stderr,"  ## Error: unable to set a new local parameter.\n");
    fprintf(stderr,"    max number of local parameters: %d\n",mesh->info.npar);
    return(0);
  }
  if ( typ != MMG5_Triangle ) {
    fprintf(stdout,"  ## Warning: you must apply your local parameters");
    fprintf(stdout," on triangles (MMG5_Triangle or %d).\n",MMG5_Triangle);
    fprintf(stdout,"  ## Unknown type of entity: ignored.\n");
    return(0);
  }
  if ( ref < 0 ) {
    fprintf(stderr,"  ## Error: negative references are not allowed.\n");
    return(0);
  }

  for (k=0; k<mesh->info.npari; k++) {
    par = &mesh->info.par[k];

    if ( par->elt == typ && par->ref == ref ) {
      par->hausd = hausd;
      par->hmin  = hmin;
      par->hmax  = hmax;
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
          fprintf(stdout,"  ## Warning: new parameters (hausd, hmin and hmax)");
          fprintf(stdout," for entities of type %d and of ref %d\n",typ,ref);
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
  case ( MMG5_Vertex ):
    mesh->info.parTyp |= MG_Vert;
    break;
  case ( MMG5_Triangle ):
    mesh->info.parTyp |= MG_Tria;
    break;
  }

  mesh->info.npari++;

  return(1);
}

void MMGS_Free_all(const int starter,...)
{

  va_list argptr;

  va_start(argptr, starter);

  _MMGS_Free_all_var(argptr);

  va_end(argptr);

  return;
}

void MMGS_Free_structures(const int starter,...)
{

  va_list argptr;

  va_start(argptr, starter);

  _MMGS_Free_structures_var(argptr);

  va_end(argptr);

  return;
}

void MMGS_Free_names(const int starter,...)
{

  va_list argptr;

  va_start(argptr, starter);

  _MMGS_Free_names_var(argptr);

  va_end(argptr);

  return;
}
