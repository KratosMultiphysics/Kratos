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
 * \file mmg2d/API_functions_2d.c
 * \brief C API functions definitions for MMG2D library.
 * \author Algiane Froehly (Bx INP/Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 *
 * \note This file contains some internal functions for the API, see
 * the \ref mmg2d/libmmg2d.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * C API for MMG2D library.
 *
 */


#include "mmg2d.h"

void MMG2D_Init_mesh(const int starter,...) {
  va_list argptr;

  va_start(argptr, starter);

  _MMG2D_Init_mesh_var(argptr);

  va_end(argptr);

  return;
}

void MMG2D_Init_fileNames(MMG5_pMesh mesh,MMG5_pSol sol
  ) {

  MMG5_Init_fileNames(mesh,sol);
  return;
}

int MMG2D_Set_inputMeshName(MMG5_pMesh mesh, const char* meshin) {

  return(MMG5_Set_inputMeshName(mesh,meshin));
}

int MMG2D_Set_inputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solin) {
  return(MMG5_Set_inputSolName(mesh,sol,solin));
}

int MMG2D_Set_outputMeshName(MMG5_pMesh mesh, const char* meshout) {

  return(MMG5_Set_outputMeshName(mesh,meshout));
}

int MMG2D_Set_outputSolName(MMG5_pMesh mesh,MMG5_pSol sol, const char* solout) {
  return(MMG5_Set_outputSolName(mesh,sol,solout));
}
void MMG2D_Init_parameters(MMG5_pMesh mesh) {

  /* Init common parameters for mmg2d, mmgs and mmg2d. */
  _MMG5_Init_parameters(mesh);

  /* default values for integers */
  /* MMG2D_IPARAM_lag = -1 */
  mesh->info.lag      = -1;
  /* MMG2D_IPARAM_optim = 0 */
  mesh->info.optim    =  0;
  /* MMG2D_IPARAM_nosurf = 0 */
  mesh->info.nosurf   =  0;  /* [0/1]    ,avoid/allow surface modifications */

  mesh->info.renum    = 0;   /* [0]    , Turn on/off the renumbering using SCOTCH; */
  mesh->info.nreg     = 0;
  /* default values for doubles */
  mesh->info.ls       = 0.0;      /* level set value */
  mesh->info.hgrad    = 1.3;      /* control gradation; */

  mesh->info.dhd  = 135.;

  //mesh->info.imprim = -7;

  /* MMG2D_IPARAM_bucket = 64 */
  mesh->info.octree = 64;
}

int MMG2D_Set_iparameter(MMG5_pMesh mesh, MMG5_pSol sol, int iparam, int val){

  switch ( iparam ) {
    /* Integer parameters */
  case MMG2D_IPARAM_verbose :
    mesh->info.imprim   = val;
    break;
  case MMG2D_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf(stdout,"  ## Warning: maximal memory authorized must be strictly positive.\n");
      fprintf(stdout,"  Reset to default value.\n");
    }
    else
      mesh->info.mem      = val;
    _MMG2D_memOption(mesh);
    if(mesh->np && (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt )) {
      return(0);
    } else if(mesh->info.mem < 39)
      return(0);
    break;
  case MMG2D_IPARAM_bucket :
    mesh->info.octree   = val;
    break;
  case MMG2D_IPARAM_debug :
    mesh->info.ddebug   = val;
    break;
  case MMG2D_IPARAM_angle :
    /* free table that may contains old ridges */
    if ( mesh->htab.geom )
      _MMG5_DEL_MEM(mesh,mesh->htab.geom,(mesh->htab.max+1)*sizeof(MMG5_hgeom));
    if ( mesh->xpoint )
      _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
    if ( mesh->xtetra )
      _MMG5_DEL_MEM(mesh,mesh->xtetra,(mesh->xtmax+1)*sizeof(MMG5_xTetra));
    if ( !val )
      mesh->info.dhd    = -1.;
    else {
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug )
        fprintf(stdout,"  ## Warning: angle detection parameter set to default value\n");
      mesh->info.dhd    = _MMG5_ANGEDG;
    }
    break;
  case MMG2D_IPARAM_iso :
    mesh->info.iso      = val;
    break;
  case MMG2D_IPARAM_lag :
    if ( val < 0 || val > 2 )
      exit(EXIT_FAILURE);
    mesh->info.lag = val;
    break;
  case MMG2D_IPARAM_msh :
    mesh->info.nreg = val;
    break;
  case MMG2D_IPARAM_numsubdomain :
    mesh->info.renum = val;
    break;
  case MMG2D_IPARAM_noinsert :
    mesh->info.noinsert = val;
    break;
  case MMG2D_IPARAM_noswap :
    mesh->info.noswap   = val;
    break;
  case MMG2D_IPARAM_nomove :
    mesh->info.nomove   = val;
    break;
  case MMG2D_IPARAM_nosurf :
    mesh->info.nosurf   = val;
    break;
  default :
    fprintf(stdout,"  ## Error: unknown type of parameter\n");
    return(0);
  }
  /* other options */
  mesh->info.fem      = 0;
  return(1);
}

int MMG2D_Set_dparameter(MMG5_pMesh mesh, MMG5_pSol sol, int dparam, double val){

  switch ( dparam ) {
    /* double parameters */
  case MMG2D_DPARAM_angleDetection :
    mesh->info.dhd = val;
    mesh->info.dhd = MG_MAX(0.0, MG_MIN(180.0,mesh->info.dhd));
    mesh->info.dhd = 180. - mesh->info.dhd;
    break;
  case MMG2D_DPARAM_hmin :
    mesh->info.hmin     = val;
    break;
  case MMG2D_DPARAM_hmax :
    mesh->info.hmax     = val;
    break;
  case MMG2D_DPARAM_hgrad :
    mesh->info.hgrad    = val;
    if ( mesh->info.hgrad < 0.0 )
      mesh->info.hgrad = -1.0;
    break;
  case MMG2D_DPARAM_hausd :
    if ( val <=0 ) {
      fprintf(stdout,"  ## Error: hausdorff number must be strictly positive.\n");
      return(0);
    }
    else
      mesh->info.hausd    = val;
    break;
  case MMG2D_DPARAM_ls :
    mesh->info.ls       = val;
    fprintf(stdout,"  ## Warning: unstable feature.\n");
    break;
  default :
    fprintf(stdout,"  ## Error: unknown type of parameter\n");
    return(0);
  }
  return(1);
}

int MMG2D_Set_meshSize(MMG5_pMesh mesh, int np, int nt, int na) {
  int k;

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->tria || mesh->edge) )
    fprintf(stdout,"  ## Warning: new mesh\n");

  if ( mesh->point )
    _MMG5_DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(MMG5_Point));
  if ( mesh->tria )
    _MMG5_DEL_MEM(mesh,mesh->tria,(mesh->ntmax+1)*sizeof(MMG5_Tria));
  if ( mesh->edge )
    _MMG5_DEL_MEM(mesh,mesh->edge,(mesh->namax+1)*sizeof(MMG5_Edge));

  mesh->np  = np;
  mesh->nt  = nt;
  mesh->na  = na;
  mesh->npi = mesh->np;
  mesh->nti = mesh->nt;
  mesh->nai = mesh->na;

  /*tester si -m definie : renvoie 0 si pas ok et met la taille min dans info.mem */
  if( mesh->info.mem > 0) {
    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->namax < mesh->na) ) {
      _MMG2D_memOption(mesh);
      //     printf("pas de pbs ? %d %d %d %d %d %d -- %d\n",mesh->npmax,mesh->np,
      //     mesh->ntmax,mesh->nt,mesh->nemax,mesh->ne,mesh->info.mem);
      if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt)) {
        fprintf(stdout,"mem insuffisante np : %d %d nt : %d %d \n"
                ,mesh->npmax,mesh->np,
                mesh->ntmax,mesh->nt);
        return(0);
      }
      else
        return(1);
    } else if(mesh->info.mem < 39) {
      printf("mem insuffisante %d\n",mesh->info.mem);
      return(0);
    }
  } else {
    mesh->npmax = MG_MAX(1.5*mesh->np,_MMG2D_NPMAX);
    mesh->ntmax = MG_MAX(1.5*mesh->nt,_MMG2D_NEMAX);

  }
  _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                printf("  Exit program.\n");
                exit(EXIT_FAILURE));
  _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point);

  _MMG5_ADD_MEM(mesh,(mesh->ntmax+1)*sizeof(MMG5_Tria),"initial triangles",return(0));
  _MMG5_SAFE_CALLOC(mesh->tria,mesh->ntmax+1,MMG5_Tria);

  mesh->namax =  MG_MAX(mesh->na,_MMG2D_NEDMAX);
  _MMG5_ADD_MEM(mesh,(mesh->namax+1)*sizeof(MMG5_Edge),"initial edges",return(0));
  _MMG5_SAFE_CALLOC(mesh->edge,(mesh->namax+1),MMG5_Edge);

  /* keep track of empty links */
  mesh->npnil = mesh->np + 1;
  mesh->nenil = mesh->nt + 1;
  mesh->nanil = mesh->na + 1;

  for (k=mesh->npnil; k<mesh->npmax-1; k++) {
    mesh->point[k].tmp  = k+1;
  }
  for (k=mesh->nenil; k<mesh->ntmax-1; k++) {
    mesh->tria[k].v[2] = k+1;
  }
  for (k=mesh->nanil; k<mesh->namax-1; k++) {
    mesh->edge[k].b = k+1;
  }

  if ( !mesh->nt ) {
    fprintf(stdout,"  **WARNING NO GIVEN TRIANGLE\n");
  }

  /* stats */
  if ( abs(mesh->info.imprim) > 6 ) {
    fprintf(stdout,"     NUMBER OF VERTICES     %8d\n",mesh->np);
    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES        %8d\n",mesh->na);
    }
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES    %8d\n",mesh->nt);
  }
  return(1);
}

int MMG2D_Set_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int typEntity, int np, int typSol) {

  if ( ( (mesh->info.imprim > 5) || mesh->info.ddebug ) && sol->m )
    fprintf(stdout,"  ## Warning: new solution\n");

  if ( typEntity != MMG5_Vertex ) {
    fprintf(stdout,"  ## Error: MMG2D need a solution imposed on vertices\n");
    return(0);
  }
  if ( typSol == MMG5_Scalar ) {
    sol->size = 1;
  }
  else if ( typSol == MMG5_Tensor ) {
    sol->size = 3;
  }
  else {
    fprintf(stdout,"  ## Error: type of solution not yet implemented\n");
    return(0);
  }

  sol->dim = 2;
  if ( np ) {
    sol->np  = np;
    sol->npi = np;
    if ( sol->m )
      _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));

    sol->npmax = mesh->npmax;
    _MMG5_ADD_MEM(mesh,(sol->size*(sol->npmax+1))*sizeof(double),"initial solution",
                  printf("  Exit program.\n");
                  exit(EXIT_FAILURE));
    _MMG5_SAFE_CALLOC(sol->m,(sol->size*(sol->npmax+1)),double);
  }
  return(1);
}

int MMG2D_Get_solSize(MMG5_pMesh mesh, MMG5_pSol sol, int* typEntity, int* np,
                      int* typSol) {

  if ( typEntity != NULL )
    *typEntity = MMG5_Vertex;

  if ( typSol != NULL ) {
    if ( sol->size == 1 )
      *typSol    = MMG5_Scalar;
    else if ( sol->size == 2 )
      *typSol    = MMG5_Vector;
    else if ( sol->size == 3 )
      *typSol    = MMG5_Tensor;
    else
      *typSol    = MMG5_Notype;
  }

  assert( (!sol->np) || (sol->np == mesh->np));

  if ( np != NULL )
    *np = sol->np;

  return(1);
}

int MMG2D_Get_meshSize(MMG5_pMesh mesh, int* np, int* nt, int* na) {
  int k;

  if ( np != NULL )
    *np = mesh->np;
  if ( nt != NULL )
    *nt = mesh->nt;

  if ( na != NULL ) {
    // Edges are not packed, thus we must count it.
    *na = 0;
    if ( mesh->na ) {
      for (k=1; k<=mesh->na; k++) {
        if ( mesh->edge[k].a ) ++(*na);
      }
    }
  }

  return(1);
}

int MMG2D_Set_vertex(MMG5_pMesh mesh, double c0, double c1, int ref, int pos) {

  if ( !mesh->np ) {
    fprintf(stdout,"  ## Error: You must set the number of points with the");
    fprintf(stdout," MMG2D_Set_meshSize function before setting vertices in mesh\n");
    return(0);
  }

  if ( pos > mesh->npmax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new point.\n");
    fprintf(stdout,"    max number of points: %d\n",mesh->npmax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->np ) {
    fprintf(stdout,"  ## Error: attempt to set new vertex at position %d.",pos);
    fprintf(stdout," Overflow of the given number of vertices: %d\n",mesh->np);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the vertex.\n");
    return(0);
  }

  mesh->point[pos].c[0] = c0;
  mesh->point[pos].c[1] = c1;
  mesh->point[pos].ref  = ref;
  if ( mesh->nt )
    mesh->point[pos].tag  = MG_NUL;
  else
    mesh->point[pos].tag  &= ~MG_NUL;

  mesh->point[pos].flag = 0;
  mesh->point[pos].tmp = 0;

  return(1);
}

/* int MMG2D_Set_corner(MMG5_pMesh mesh, int k) { */
/*   assert ( k <= mesh->np ); */
/*   mesh->point[k].tag |= M_CORNER; */
/*   return(1); */
/* } */

/* int MMG2D_Set_requiredVertex(MMG5_pMesh mesh, int k) { */
/*   assert ( k <= mesh->np ); */
/*   mesh->point[k].tag |= M_REQUIRED; */
/*   return(1); */
/* } */

int MMG2D_Get_vertex(MMG5_pMesh mesh, double* c0, double* c1, int* ref,
                    int* isCorner, int* isRequired) {

 if ( mesh->npi == mesh->np ) {
   mesh->npi = 0;
   if ( mesh->info.ddebug ) {
    fprintf(stdout,"  ## Warning: reset the internal counter of points.\n");
    fprintf(stdout,"     You must pass here exactly one time (the first time ");
    fprintf(stdout,"you call the MMG2D_Get_vertex function).\n");
    fprintf(stdout,"     If not, the number of call of this function");
    fprintf(stdout," exceed the number of points: %d\n ",mesh->np);
   }
 }

  mesh->npi++;

  if ( mesh->npi > mesh->np ) {
    fprintf(stdout,"  ## Error: unable to get point.\n");
    fprintf(stdout,"     The number of call of MMG2D_Get_vertex function");
    fprintf(stdout," exceed the number of points: %d\n ",mesh->np);
    return(0);
  }

  *c0  = mesh->point[mesh->npi].c[0];
  *c1  = mesh->point[mesh->npi].c[1];
  if ( ref != NULL )
    *ref = mesh->point[mesh->npi].ref;

  if ( isCorner != NULL ) {
    if ( mesh->point[mesh->npi].tag & M_CORNER )
      *isCorner = 1;
    else
      *isCorner = 0;
  }

  if ( isRequired != NULL ) {
    if ( mesh->point[mesh->npi].tag & M_REQUIRED )
      *isRequired = 1;
    else
      *isRequired = 0;
  }

  return(1);
}

int  MMG2D_Set_vertices(MMG5_pMesh mesh, double *vertices,int *refs) {
  MMG5_pPoint ppt;
  int i,j;

  /*coordinates vertices*/
  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*2;
    ppt->c[0]  = vertices[j];
    ppt->c[1]  = vertices[j+1];

    ppt->flag = 0;
    ppt->tmp = 0;

    if ( refs != NULL )
      ppt->ref   = refs[i-1];

    if ( mesh->nt )
      ppt->tag  = MG_NUL;
    else
      ppt->tag  &= ~MG_NUL;
  }

  return 1;
}

int  MMG2D_Get_vertices(MMG5_pMesh mesh, double* vertices, int* refs,
                        int* areCorners, int* areRequired) {
  MMG5_pPoint ppt;
  int i,j;

  for (i=1;i<=mesh->np;i++)
  {
    ppt = &mesh->point[i];

    j = (i-1)*2;
    vertices[j] = ppt->c[0];
    vertices[j+1] = ppt->c[1];

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

int MMG2D_Set_triangle(MMG5_pMesh mesh, int v0, int v1, int v2, int ref, int pos) {
  MMG5_pPoint ppt;
  MMG5_pTria  pt;
  double      vol;
  int         i,j,ip,tmp;

  if ( !mesh->nt ) {
    fprintf(stdout,"  ## Error: You must set the number of elements with the");
    fprintf(stdout," MMG2D_Set_meshSize function before setting elements in mesh\n");
    return(0);
  }

  if ( pos > mesh->ntmax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new element.\n");
    fprintf(stdout,"    max number of element: %d\n",mesh->ntmax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->nt ) {
    fprintf(stdout,"  ## Error: attempt to set new triangle at position %d.",pos);
    fprintf(stdout," Overflow of the given number of triangle: %d\n",mesh->nt);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the triangle.\n");
    return(0);
  }

  pt = &mesh->tria[pos];
  pt->v[0] = v0;
  pt->v[1] = v1;
  pt->v[2] = v2;
  pt->ref  = ref;

  mesh->point[pt->v[0]].tag &= ~MG_NUL;
  mesh->point[pt->v[1]].tag &= ~MG_NUL;
  mesh->point[pt->v[2]].tag &= ~MG_NUL;

  for(i=0 ; i<3 ; i++)
    pt->edg[i] = 0;

  vol = MMG2_quickarea(mesh->point[pt->v[0]].c,mesh->point[pt->v[1]].c,
                           mesh->point[pt->v[2]].c);

  if ( vol == 0.0 ) {
    fprintf(stderr,"  ## Error: triangle %d has null area.\n",pos);
    for ( ip=0; ip<3; ip++ ) {
      ppt = &mesh->point[pt->v[ip]];
      for ( j=0; j<3; j++ ) {
        if ( fabs(ppt->c[j])>0. ) {
          fprintf(stderr," Check that you don't have a sliver triangle.\n");
          return(0);
        }
      }
    }
  }
  else if(vol < 0) {
    printf("Tr %d bad oriented\n",pos);
    tmp = pt->v[2];
    pt->v[2] = pt->v[1];
    pt->v[1] = tmp;
    /* mesh->xt temporary used to count reoriented tetra */
    mesh->xt++;
  }
  if ( mesh->info.ddebug && (mesh->nt == pos) && mesh->xt > 0 ) {
    fprintf(stdout,"  ## %d triangles reoriented\n",mesh->xt);
    mesh->xt = 0;
  }

  return(1);
}

/* int MMG2D_Set_requiredTriangle(MMG5_pMesh mesh, int k) { */
/*   assert ( k <= mesh->nt ); */
/*   mesh->tria[k].tag[0] |= M_REQUIRED; */
/*   mesh->tria[k].tag[1] |= M_REQUIRED; */
/*   mesh->tria[k].tag[2] |= M_REQUIRED; */
/*   return(1); */
/* } */

int MMG2D_Get_triangle(MMG5_pMesh mesh, int* v0, int* v1, int* v2, int* ref
                       ,int* isRequired) {
  MMG5_pTria  ptt;

  if ( mesh->nti == mesh->nt ) {
    mesh->nti = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of triangles.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMG2D_Get_triangle function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of triangles: %d\n ",mesh->nt);
    }
  }

  mesh->nti++;

  if ( mesh->nti > mesh->nt ) {
    fprintf(stdout,"  ## Error: unable to get triangle.\n");
    fprintf(stdout,"    The number of call of MMG2D_Get_triangle function");
    fprintf(stdout," can not exceed the number of triangles: %d\n ",mesh->nt);
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

int  MMG2D_Set_triangles(MMG5_pMesh mesh, int *tria, int *refs) {
  MMG5_pPoint ppt;
  MMG5_pTria ptt;
  double vol;
  int i, j,ip,tmp;

  mesh->xt = 0;
  for (i=1;i<=mesh->nt;i++)
  {
      j = (i-1)*3;
      ptt = &mesh->tria[i];
      ptt->v[0] = tria[j]  ;
      ptt->v[1] = tria[j+2];
      ptt->v[2] = tria[j+1];
      if ( refs != NULL )
        ptt->ref  = refs[i-1];

      mesh->point[ptt->v[0]].tag &= ~MG_NUL;
      mesh->point[ptt->v[1]].tag &= ~MG_NUL;
      mesh->point[ptt->v[2]].tag &= ~MG_NUL;

      for(i=0 ; i<3 ; i++)
        ptt->edg[i] = 0;

      vol = MMG2_quickarea(mesh->point[ptt->v[0]].c,mesh->point[ptt->v[1]].c,
                           mesh->point[ptt->v[2]].c);

      if ( vol == 0.0 ) {
        fprintf(stderr,"  ## Error: triangle %d has null area.\n",i);
        for ( ip=0; ip<3; ip++ ) {
          ppt = &mesh->point[ptt->v[ip]];
          for ( j=0; j<3; j++ ) {
            if ( fabs(ppt->c[j])>0. ) {
              fprintf(stderr," Check that you don't have a sliver triangle.\n");
              return(0);
            }
          }
        }
      }
      else if(vol < 0) {
        printf("Tr %d bad oriented\n",i);
        tmp = ptt->v[2];
        ptt->v[2] = ptt->v[1];
        ptt->v[1] = tmp;
        /* mesh->xt temporary used to count reoriented tetra */
        mesh->xt++;
      }
      if ( mesh->info.ddebug && mesh->xt > 0 ) {
        fprintf(stdout,"  ## %d triangles reoriented\n",mesh->xt);
      }
  }
  return 1;
}

int  MMG2D_Get_triangles(MMG5_pMesh mesh, int* tria, int* refs,
                         int* areRequired) {
  MMG5_pTria ptt;
  int        i, j;

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

int MMG2D_Set_edge(MMG5_pMesh mesh, int v0, int v1, int ref, int pos) {
  MMG5_pEdge pt;

  if ( !mesh->na ) {
    fprintf(stdout,"  ## Error: You must set the number of elements with the");
    fprintf(stdout," MMG2D_Set_meshSize function before setting elements in mesh\n");
    return(0);
  }

  if ( pos > mesh->namax ) {
    fprintf(stdout,"  ## Error: unable to allocate a new element.\n");
    fprintf(stdout,"    max number of element: %d\n",mesh->namax);
    _MMG5_INCREASE_MEM_MESSAGE();
    return(0);
  }

  if ( pos > mesh->na ) {
    fprintf(stdout,"  ## Error: attempt to set new edge at position %d.",pos);
    fprintf(stdout," Overflow of the given number of edge: %d\n",mesh->na);
    fprintf(stdout,"  ## Check the mesh size, its compactness or the position");
    fprintf(stdout," of the edge.\n");
    return(0);
  }

  pt = &mesh->edge[pos];
  pt->a = v0;
  pt->b = v1;
  pt->ref  = ref;

  mesh->point[pt->a].tag &= ~MG_NUL;
  mesh->point[pt->b].tag &= ~MG_NUL;

  return(1);
}

int MMG2D_Set_requiredEdge(MMG5_pMesh mesh, int k) {
  MMG5_pPoint ppt;
  MMG5_pEdge  ped;

  assert ( k <= mesh->na );

  ped = &mesh->edge[k];

  ped->tag |= M_REQUIRED;

  ppt = &mesh->point[ped->a];
  ppt->tag |= M_REQUIRED;
  ppt = &mesh->point[ped->b];
  ppt->tag |= M_REQUIRED;

  return(1);
}

int MMG2D_Get_edge(MMG5_pMesh mesh, int* e0, int* e1, int* ref
                   ,int* isRidge, int* isRequired) {
  MMG5_pEdge        ped;

  if ( mesh->nai == mesh->na ) {
    mesh->nai = 0;
    if ( mesh->info.ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of edges.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMG2D_Get_edge function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of edges.\n ");
      fprintf(stdout,"     Please, call the MMG2D_Get_meshSize function to get"
              " this number.\n ");
    }
  }

  mesh->nai++;

  if ( mesh->nai > mesh->na ) {
    fprintf(stdout,"  ## Error: unable to get edge.\n");
    fprintf(stdout,"    The number of call of MMG2D_Get_edge function");
    fprintf(stdout," can not exceed the number of edges: %d\n ",mesh->na);
    return(0);
  }

  ped = &mesh->edge[mesh->nai];

  while ( !ped->a && ++mesh->nai <= mesh->na ) {
    ped = &mesh->edge[mesh->nai];
  }


  *e0  = ped->a;
  *e1  = ped->b;

  if ( ref != NULL )
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

int MMG2D_Set_scalarSol(MMG5_pSol met, double s, int pos) {

  if ( !met->np ) {
    fprintf(stdout,"  ## Error: You must set the number of solution with the");
    fprintf(stdout," MMG2D_Set_solSize function before setting values");
    fprintf(stdout," in solution structure \n");
    return(0);
  }

  if ( pos >= met->npmax ) {
    fprintf(stdout,"  ## Error: unable to set a new solution.\n");
    fprintf(stdout,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stdout,"  ## Error: attempt to set new solution at position %d.",pos);
    fprintf(stdout," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stdout,"  ## Check the solution size, its compactness or the position");
    fprintf(stdout," of the solution.\n");
    return(0);
  }

  met->m[pos] = s;
  return(1);
}

int  MMG2D_Get_scalarSol(MMG5_pSol met, double* s)
{
  int ddebug = 0;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of points.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMG2D_Get_scalarSol function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of points: %d\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stdout,"  ## Error: unable to get solution.\n");
    fprintf(stdout,"     The number of call of MMG2D_Get_scalarSol function");
    fprintf(stdout," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  *s  = met->m[met->npi];

  return(1);
}

int  MMG2D_Set_scalarSols(MMG5_pSol met, double *s) {
  int k;

  if ( !met->np ) {
    fprintf(stderr,"  ## Error: You must set the number of solution with the");
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }

  for ( k=0; k<met->np; ++k )
    met->m[k+1] = s[k];

  return 1;
}

int  MMG2D_Get_scalarSols(MMG5_pSol met, double* s) {
  int k;

  for ( k=0; k<met->np; ++k )
    s[k]  = met->m[k+1];

  return(1);
}

int MMG2D_Set_tensorSol(MMG5_pSol met, double m11, double m12, double m22,
                        int pos) {
  int isol;

  if ( !met->np ) {
    fprintf(stdout,"  ## Error: You must set the number of solution with the");
    fprintf(stdout," MMG2D_Set_solSize function before setting values");
    fprintf(stdout," in solution structure \n");
    return(0);
  }

  if ( pos >= met->npmax ) {
    fprintf(stdout,"  ## Error: unable to set a new solution.\n");
    fprintf(stdout,"    max number of solutions: %d\n",met->npmax);
    return(0);
  }

  if ( pos > met->np ) {
    fprintf(stdout,"  ## Error: attempt to set new solution at position %d.",pos);
    fprintf(stdout," Overflow of the given number of solutions: %d\n",met->np);
    fprintf(stdout,"  ## Check the solution size, its compactness or the position");
    fprintf(stdout," of the solution.\n");
    return(0);
  }
  isol = pos * met->size;
  met->m[isol    ] = m11;
  met->m[isol + 1] = m12;
  met->m[isol + 2] = m22;
  return(1);
}

int MMG2D_Get_tensorSol(MMG5_pSol met, double *m11,double *m12,double *m22)
{
  int ddebug = 0;
  int isol;

  if ( met->npi == met->np ) {
    met->npi = 0;
    if ( ddebug ) {
      fprintf(stdout,"  ## Warning: reset the internal counter of points.\n");
      fprintf(stdout,"     You must pass here exactly one time (the first time ");
      fprintf(stdout,"you call the MMG2D_Get_tensorSol function).\n");
      fprintf(stdout,"     If not, the number of call of this function");
      fprintf(stdout," exceed the number of points: %d\n ",met->np);
    }
  }

  met->npi++;

  if ( met->npi > met->np ) {
    fprintf(stdout,"  ## Error: unable to get solution.\n");
    fprintf(stdout,"     The number of call of MMG2D_Get_tensorSol function");
    fprintf(stdout," can not exceed the number of points: %d\n ",met->np);
    return(0);
  }

  isol = met->size*(met->npi);
  *m11 = met->m[isol  ];
  *m12 = met->m[isol+1];
  *m22 = met->m[isol+2];

  return(1);
}

int MMG2D_Set_tensorSols(MMG5_pSol met, double *sols) {
  double *m;
  int k,j;

  if ( !met->np ) {
    fprintf(stderr,"  ## Error: You must set the number of solution with the");
    fprintf(stderr," MMG3D_Set_solSize function before setting values");
    fprintf(stderr," in solution structure \n");
    return(0);
  }

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j];

    m[1] = sols[j];
    m[2] = sols[j+1];
    m[3] = sols[j+2];
  }
  return(1);
}

int MMG2D_Get_tensorSols(MMG5_pSol met, double *sols) {
  double *m;
  int k,j;

  for ( k=0; k<met->np; ++k ) {
    j = 3*k;
    m = &met->m[j];

    sols[j]   = m[1];
    sols[j+1] = m[2];
    sols[j+2] = m[3];
  }

  return(1);
}

int MMG2D_Chk_meshData(MMG5_pMesh mesh,MMG5_pSol met) {

  if ( (mesh->npi != mesh->np) || (mesh->nti != mesh->nt) ) {
    fprintf(stdout,"  ## Error: if you don't use the MMG2D_loadMesh function,");
    fprintf(stdout," you must call the MMG2D_Set_meshSize function to have a");
    fprintf(stdout," valid mesh.\n");
    fprintf(stdout," Missing datas.\n");
    return(0);
  }

  if ( met->npi != met->np ) {
    fprintf(stdout,"  ## Error: if you don't use the MMG2D_loadMet function,");
    fprintf(stdout," you must call the MMG2D_Set_solSize function to have a");
    fprintf(stdout," valid solution.\n");
    fprintf(stdout," Missing datas.\n");
    return(0);
  }

  /*  Check mesh data */
  if ( mesh->info.ddebug ) {
    if ( (!mesh->np) || (!mesh->point) ||
         (!mesh->nt)  ) {
      fprintf(stdout,"  ** MISSING DATA.\n");
      fprintf(stdout," Check that your mesh contains points.\n");
      fprintf(stdout," Exit program.\n");
      return(0);
    }
  }

  if ( mesh->dim != 2 ) {
    fprintf(stdout,"  ** 2 DIMENSIONAL MESH NEEDED. Exit program.\n");
    return(0);
  }
  if ( met->dim != 2 ) {
    fprintf(stdout,"  ** WRONG DIMENSION FOR METRIC. Exit program.\n");
    return(0);
  }
  if ( !mesh->ver )  mesh->ver = 2;
  if ( !met ->ver )  met ->ver = 2;

  return(1);
}

void MMG2D_Free_all(const int starter,...)
{

  va_list argptr;

  va_start(argptr, starter);

  _MMG2D_Free_all_var(argptr);

  va_end(argptr);

  return;
}

void MMG2D_Free_structures(const int starter,...)
{

  va_list argptr;

  va_start(argptr, starter);

  _MMG2D_Free_structures_var(argptr);

  va_end(argptr);

  return;
}

void MMG2D_Free_names(const int starter,...)
{

  va_list argptr;

  va_start(argptr, starter);

  _MMG2D_Free_names_var(argptr);

  va_end(argptr);

  return;
}
