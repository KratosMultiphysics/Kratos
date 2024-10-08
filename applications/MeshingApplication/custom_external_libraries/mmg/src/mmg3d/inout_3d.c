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
 * \file mmg3d/inout_3d.c
 * \brief Input / Output Functions.
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
 * \param imprim verbosity level (muted for stdout if -1)
 * \param filename file to open
 * \param inm pointer toward the file unit
 * \param bin 1 if file will be at binary format
 * \param modeASCII mode in which to open an ascii file ("r","r+","w","w+",...)
 * \param modeASCII mode in which to open an ascii file ("r","r+","w","w+",...)
 *
 * \return 0 if fail to open file, -1 for other errors, 1 if success.
 *
 * Try to open a Medit file at asked mode (read only, write, etc) and store if
 * file is binary (depending on the extension).
 *
 */
int MMG3D_openMesh(int imprim,const char *filename,FILE **inm,int *bin,
                   char *modeASCII, char* modeBIN) {
  char *ptr,*data;
  int   out,ier;

  out = (strchr(modeASCII,'w') != NULL) ? 1 : 0;
  ier = out ? 0 : -1;

  *bin = 0;
  MMG5_SAFE_CALLOC(data,strlen(filename)+7,char,return ier);

  strcpy(data,filename);
  ptr = strstr(data,".mesh");

  if ( !ptr ) {
    /* data contains the filename without extension */
    strcat(data,".meshb");
    if( !(*inm = fopen(data,modeBIN)) ) {
      /* our file is not a .meshb file, try with .mesh ext */
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if( !(*inm = fopen(data,modeASCII)) ) {
        MMG5_SAFE_FREE(data);
        return 0;
      }
    }
    else  *bin = 1;
  }
  else {
    ptr = strstr(data,".meshb");
    if ( ptr ) {
      *bin = 1;
      if( !(*inm = fopen(data,modeBIN)) ) {
        if ( out ) {
          fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        }
        MMG5_SAFE_FREE(data);
        return 0;
      }
    }
    else {
      if( !(*inm = fopen(data,modeASCII)) ) {
        if ( out ) {
          fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        }
        MMG5_SAFE_FREE(data);
        return 0;
      }
    }
  }

  if ( imprim >= 0 ) {
    fprintf(stdout,"  %%%% %s OPENED\n",data);
  }
  MMG5_SAFE_FREE(data);

 return 1;
}

int MMG3D_loadMesh_opened(MMG5_pMesh mesh,FILE *inm,int bin) {
  MMG5_pTetra pt;
  MMG5_pPrism pp;
  MMG5_pTria  pt1;
  MMG5_pQuad  pq1;
  MMG5_pEdge  pa;
  MMG5_pPoint ppt;
  double      *norm,*n,dd;
  float       fc;
  long        posnp,posnt,posne,posned,posncor,posnpreq,posntreq,posnereq,posnedreq;
  long        posnr,posnprism,posnormal,posnc1,posnq,posnqreq;
  long        posnppar,posntpar,posnepar,posnedpar,posnqpar;
  MMG5_int    npreq,ntreq,nereq,nedreq,nqreq,ncor,ned,ng;
  MMG5_int    nppar,nedpar,ntpar,nqpar,nepar;
  MMG5_int    k,ip,idn;
  int         i,bdim,binch,iswp,bpos;
  MMG5_int    na,nr,ia,aux,nref,ref;
  char        chaine[MMG5_FILESTR_LGTH],strskip[MMG5_FILESTR_LGTH];

  posnp = posnt = posne = posncor = 0;
  posnpreq = posntreq = posnereq = posnqreq = posned = posnedreq = posnr = 0;
  posnprism = posnormal= posnc1 = posnq = 0;
  posnppar = posnedpar = posntpar = posnqpar = posnepar = 0;
  ncor = ned = npreq = ntreq = nqreq = nereq = nedreq = nr = ng = 0;
  nppar = nedpar = ntpar = nqpar = nepar = 0;
  iswp = 0;
  mesh->np = mesh->nt = mesh->ne = 0;
  nref = 0;


  if (!bin) {
    strcpy(chaine,"D");
    while(fscanf(inm,"%127s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if ( chaine[0] == '#' ) {
        fgets(strskip,MMG5_FILESTR_LGTH,inm);
        continue;
      }

      if(!strncmp(chaine,"MeshVersionFormatted",strlen("MeshVersionFormatted"))) {
        MMG_FSCANF(inm,"%d",&mesh->ver);
        continue;
      } else if(!strncmp(chaine,"Dimension",strlen("Dimension"))) {
        MMG_FSCANF(inm,"%d",&mesh->dim);
        if(mesh->dim!=3) {
          fprintf(stderr,"BAD DIMENSION : %d\n",mesh->dim);
          return -1;
        }
        continue;
      } else if(!strncmp(chaine,"Vertices",strlen("Vertices"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->npi);
        posnp = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredVertices",strlen("RequiredVertices"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&npreq);
        posnpreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"ParallelVertices",strlen("ParallelVertices"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nppar);
        posnppar = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Triangles",strlen("Triangles"))) {
        if ( !strncmp(chaine,"TrianglesP",strlen("TrianglesP")) ) continue;
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->nti);
        posnt = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredTriangles",strlen("RequiredTriangles"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&ntreq);
        posntreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"ParallelTriangles",strlen("ParallelTriangles"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&ntpar);
        posntpar = ftell(inm);
        continue;
      }
      else if(!strncmp(chaine,"Quadrilaterals",strlen("Quadrilaterals"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->nquad);
        posnq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredQuadrilaterals",strlen("RequiredQuadrilaterals"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nqreq);
        posnqreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"ParallelQuadrilaterals",strlen("ParallelQuadrilaterals"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nqpar);
        posnqpar = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Tetrahedra",strlen("Tetrahedra"))) {
        if ( !strncmp(chaine,"TetrahedraP",strlen("TetrahedraP")) ) continue;
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->nei);
        posne = ftell(inm);
        continue;
      } else if((!strncmp(chaine,"Prisms",strlen("Prisms")))||
                (!strncmp(chaine,"Pentahedra",strlen("Pentahedra")))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->nprism);
        posnprism = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredTetrahedra",strlen("RequiredTetrahedra"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nereq);
        posnereq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"ParallelTetrahedra",strlen("ParallelTetrahedra"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nepar);
        posnepar = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Corners",strlen("Corners"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&ncor);
        posncor = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Edges",strlen("Edges"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->nai);
        posned = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredEdges",strlen("RequiredEdges"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nedreq);
        posnedreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"ParallelEdges",strlen("ParallelEdges"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nedpar);
        posnedpar = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Ridges",strlen("Ridges"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nr);
        posnr = ftell(inm);
        continue;
      } else if(!ng && !strncmp(chaine,"Normals",strlen("Normals"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&ng);
        posnormal = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"NormalAtVertices",strlen("NormalAtVertices"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->nc1);
        posnc1 = ftell(inm);
        continue;
      }
    }
  } else { //binary file
    bdim = 0;
    MMG_FREAD(&mesh->ver,MMG5_SW,1,inm);
    iswp=0;
    if(mesh->ver==16777216)
      iswp=1;
    else if(mesh->ver!=1) {
      fprintf(stderr,"BAD FILE ENCODING\n");
    }
    MMG_FREAD(&mesh->ver,MMG5_SW,1,inm);
    if(iswp) mesh->ver = MMG5_swapbin(mesh->ver);
    while(fread(&binch,MMG5_SW,1,inm)!=0 && binch!=54 ) {
      if(iswp) binch=MMG5_swapbin(binch);
      if(binch==54) break;
      if(!bdim && binch==3) {  //Dimension
        MMG_FREAD(&bdim,MMG5_SW,1,inm);  //NulPos=>20
        if(iswp) bdim=MMG5_swapbin(bdim);
        MMG_FREAD(&bdim,MMG5_SW,1,inm);
        if(iswp) bdim=MMG5_swapbin(bdim);
        mesh->dim = bdim;
        if(bdim!=3) {
          fprintf(stderr,"BAD MESH DIMENSION : %d\n",mesh->dim);
          fprintf(stderr," Exit program.\n");
          return -1;
        }
        continue;
      } else if(!mesh->npi && binch==4) {  //Vertices
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&mesh->npi,MMG5_SW,1,inm);
        if(iswp) mesh->npi=MMG5_SWAPBIN(mesh->npi);
        posnp = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==15) {  //RequiredVertices
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&npreq,MMG5_SW,1,inm);
        if(iswp) npreq=MMG5_SWAPBIN(npreq);
        posnpreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->nti && binch==6) {//Triangles
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&mesh->nti,MMG5_SW,1,inm);
        if(iswp) mesh->nti=MMG5_SWAPBIN(mesh->nti);
        posnt = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==17) {  //RequiredTriangles
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&ntreq,MMG5_SW,1,inm);
        if(iswp) ntreq=MMG5_SWAPBIN(ntreq);
        posntreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      }
      else if(!mesh->nquad && binch==7) {//Quadrilaterals
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&mesh->nquad,MMG5_SW,1,inm);
        if(iswp) mesh->nquad=MMG5_SWAPBIN(mesh->nquad);
        posnq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==18) {  //RequiredQuadrilaterals
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&nqreq,MMG5_SW,1,inm);
        if(iswp) nqreq=MMG5_SWAPBIN(nqreq);
        posnqreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->nei && binch==8) {//Tetra
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&mesh->nei,MMG5_SW,1,inm);
        if(iswp) mesh->nei=MMG5_SWAPBIN(mesh->nei);
        posne = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->nprism && binch==9) {//Prism
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&mesh->nprism,MMG5_SW,1,inm);
        if(iswp) mesh->nprism=MMG5_SWAPBIN(mesh->nprism);
        posnprism = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==12) {  //RequiredTetra
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&nereq,MMG5_SW,1,inm);
        if(iswp) nereq=MMG5_SWAPBIN(nereq);
        posnereq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!ncor && binch==13) { //Corners
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&ncor,MMG5_SW,1,inm);
        if(iswp) ncor=MMG5_SWAPBIN(ncor);
        posncor = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->nai && binch==5) { //Edges
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&mesh->nai,MMG5_SW,1,inm);
        if(iswp) mesh->nai=MMG5_SWAPBIN(mesh->nai);
        posned = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==16) {  //RequiredEdges
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&nedreq,MMG5_SW,1,inm);
        if(iswp) nedreq=MMG5_SWAPBIN(nedreq);
        posnedreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      }  else if(binch==14) {  //Ridges
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&nr,MMG5_SW,1,inm);
        if(iswp) nr=MMG5_SWAPBIN(nr);
        posnr = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!ng && binch==60) {  //Normals
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&ng,MMG5_SW,1,inm);
        if(iswp) ng=MMG5_SWAPBIN(ng);
        posnormal = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(binch==20) {  //NormalAtVertices
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&mesh->nc1,MMG5_SW,1,inm);
        if(iswp) mesh->nc1=MMG5_SWAPBIN(mesh->nc1);
        posnc1 = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else {
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);

        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }
  }

  if ( !mesh->npi || !mesh->nei ) {
    fprintf(stderr,"  ** MISSING DATA.\n");
    fprintf(stderr," Check that your mesh contains points and tetrahedra.\n");
    fprintf(stderr," Exit program.\n");
    return -1;
  }
  /* memory allocation */
  mesh->np = mesh->npi;
  mesh->nt = mesh->nti;
  mesh->ne = mesh->nei;
  mesh->na = mesh->nai;
  if ( !MMG3D_zaldy(mesh) )  return 0;
  if (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne) {
    return -1;
  }

  rewind(inm);
  fseek(inm,posnp,SEEK_SET);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if (mesh->ver < 2) { /*float*/
      if (!bin) {
        for (i=0 ; i<3 ; i++) {
          MMG_FSCANF(inm,"%f",&fc);
          ppt->c[i] = (double) fc;
        }
        MMG_FSCANF(inm,"%" MMG5_PRId "",&ppt->ref);
      } else {
        for (i=0 ; i<3 ; i++) {
          MMG_FREAD(&fc,MMG5_SW,1,inm);
          if(iswp) fc=MMG5_swapf(fc);
          ppt->c[i] = (double) fc;
        }
        MMG_FREAD(&ppt->ref,MMG5_SW,1,inm);
        if(iswp) ppt->ref=MMG5_SWAPBIN(ppt->ref);
      }
    } else {
      if (!bin) {
        MMG_FSCANF(inm,"%lf %lf %lf %" MMG5_PRId "",&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
      }
      else {
        for (i=0 ; i<3 ; i++) {
          MMG_FREAD(&ppt->c[i],MMG5_SD,1,inm);
          if(iswp) ppt->c[i]=MMG5_swapd(ppt->c[i]);
        }
        MMG_FREAD(&ppt->ref,MMG5_SW,1,inm);
        if(iswp) ppt->ref=MMG5_SWAPBIN(ppt->ref);
      }
    }

    if ( ppt->ref < 0 ) {
      ppt->ref = -ppt->ref;
      ++nref;
    }

    ppt->tag  = MG_NUL;
    ppt->tmp  = 0;
  }
  /* get required vertices */
  if(npreq) {
    rewind(inm);
    fseek(inm,posnpreq,SEEK_SET);
    for (k=1; k<=npreq; k++) {
      if(!bin) {
        MMG_FSCANF(inm,"%d",&i);
      }
      else {
        MMG_FREAD(&i,MMG5_SW,1,inm);
        if(iswp) i=MMG5_swapbin(i);
      }
      if(i>mesh->np) {
        fprintf(stderr,"\n  ## Warning: %s: required Vertices number %8d"
                " ignored.\n",__func__,i);
      } else {
        ppt = &mesh->point[i];
        ppt->tag |= MG_REQ;
        ppt->tag &= ~MG_NUL;
      }
    }
  }
  /* get parallel vertices */
  if(nppar) {
    rewind(inm);
    fseek(inm,posnppar,SEEK_SET);
    for (k=1; k<=nppar; k++) {
      if(!bin) {
        MMG_FSCANF(inm,"%d",&i);
      }
      else {
        MMG_FREAD(&i,MMG5_SW,1,inm);
        if(iswp) i=MMG5_swapbin(i);
      }
      if(i>mesh->np) {
        fprintf(stderr,"\n  ## Warning: %s: parallel Vertices number %8d"
                " ignored.\n",__func__,i);
      } else {
        ppt = &mesh->point[i];
        ppt->tag |= MG_PARBDY;
      }
    }
  }


  /* get corners */
  if(ncor) {
    rewind(inm);
    fseek(inm,posncor,SEEK_SET);
    for (k=1; k<=ncor; k++) {
      if(!bin) {
        MMG_FSCANF(inm,"%d",&i);
      }
      else {
        MMG_FREAD(&i,MMG5_SW,1,inm);
        if(iswp) i=MMG5_swapbin(i);
      }
      if(i>mesh->np) {
        fprintf(stderr,"\n  ## Warning: %s: corner number %8d ignored.\n",
                __func__,i);
      } else {
        ppt = &mesh->point[i];
        ppt->tag |= MG_CRN;
      }
    }
  }

  /* read mesh triangles */
  if ( mesh->nt ) {
    rewind(inm);
    fseek(inm,posnt,SEEK_SET);
    /* Skip triangles with mesh->info.isoref refs */
    for (k=1; k<=mesh->nt; k++) {
      pt1 = &mesh->tria[k];
      if (!bin) {
        MMG_FSCANF(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "",&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
      }
      else {
        for (i=0 ; i<3 ; i++) {
          MMG_FREAD(&pt1->v[i],MMG5_SW,1,inm);
          if(iswp) pt1->v[i]=MMG5_SWAPBIN(pt1->v[i]);
        }
        MMG_FREAD(&pt1->ref,MMG5_SW,1,inm);
        if(iswp) pt1->ref=MMG5_SWAPBIN(pt1->ref);
      }
      if ( pt1->ref < 0 ) {
        pt1->ref = -pt1->ref;
        ++nref;
      }
    }
    /* get required triangles */
    if(ntreq) {
      rewind(inm);
      fseek(inm,posntreq,SEEK_SET);
      for (k=1; k<=ntreq; k++) {
        if(!bin) {
          MMG_FSCANF(inm,"%d",&i);
        }
        else {
          MMG_FREAD(&i,MMG5_SW,1,inm);
          if(iswp) i=MMG5_swapbin(i);
        }
        if ( i>mesh->nt ) {
          fprintf(stderr,"\n  ## Warning: %s: required triangle number %8d"
                  " ignored.\n",__func__,i);
        } else {
          pt1 = &mesh->tria[i];
          pt1->tag[0] |= MG_REQ;
          pt1->tag[1] |= MG_REQ;
          pt1->tag[2] |= MG_REQ;
        }
      }
    }
    /* get parallel triangles */
    if(ntpar) {
      rewind(inm);
      fseek(inm,posntpar,SEEK_SET);
      for (k=1; k<=ntpar; k++) {
        if(!bin) {
          MMG_FSCANF(inm,"%d",&i);
        }
        else {
          MMG_FREAD(&i,MMG5_SW,1,inm);
          if(iswp) i=MMG5_swapbin(i);
        }
        if ( i>mesh->nt ) {
          fprintf(stderr,"\n  ## Warning: %s: parallel triangle number %8d"
                  " ignored.\n",__func__,i);
        } else {
          pt1 = &mesh->tria[i];
          pt1->tag[0] |= MG_PARBDY;
          pt1->tag[1] |= MG_PARBDY;
          pt1->tag[2] |= MG_PARBDY;
        }
      }
    }
  } //end if mesh->nt

  /* read mesh quadrilaterals */
  if ( mesh->nquad ) {
    rewind(inm);
    fseek(inm,posnq,SEEK_SET);

    for (k=1; k<=mesh->nquad; k++) {
      pq1 = &mesh->quadra[k];
      if (!bin) {
        MMG_FSCANF(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "",&pq1->v[0],&pq1->v[1],&pq1->v[2],
                   &pq1->v[3],&pq1->ref);
      }
      else {
        for (i=0 ; i<4 ; i++) {
          MMG_FREAD(&pq1->v[i],MMG5_SW,1,inm);
          if(iswp) pq1->v[i]=MMG5_SWAPBIN(pq1->v[i]);
        }
        MMG_FREAD(&pq1->ref,MMG5_SW,1,inm);
        if(iswp) pq1->ref=MMG5_SWAPBIN(pq1->ref);
      }
      if ( pq1->ref < 0 ) {
        pq1->ref = -pq1->ref;
        ++nref;
      }
    }

    /* get required quadrilaterals */
    if(nqreq) {
      rewind(inm);
      fseek(inm,posnqreq,SEEK_SET);
      for (k=1; k<=nqreq; k++) {
        if(!bin) {
          MMG_FSCANF(inm,"%d",&i);
        }
        else {
          MMG_FREAD(&i,MMG5_SW,1,inm);
          if(iswp) i=MMG5_swapbin(i);
        }
        if ( i>mesh->nquad ) {
          fprintf(stderr,"\n  ## Warning: %s: required quadrilaterals number"
                  " %8d ignored.\n",__func__,i);
        } else {
          pq1 = &mesh->quadra[i];
          pq1->tag[0] |= MG_REQ;
          pq1->tag[1] |= MG_REQ;
          pq1->tag[2] |= MG_REQ;
          pq1->tag[3] |= MG_REQ;
        }
      }
    }

    /* get parallel quadrilaterals */
    if(nqpar) {
      rewind(inm);
      fseek(inm,posnqpar,SEEK_SET);
      for (k=1; k<=nqpar; k++) {
        if(!bin) {
          MMG_FSCANF(inm,"%d",&i);
        }
        else {
          MMG_FREAD(&i,MMG5_SW,1,inm);
          if(iswp) i=MMG5_swapbin(i);
        }
        if ( i>mesh->nquad ) {
          fprintf(stderr,"\n  ## Warning: %s: parallel quadrilaterals number"
                  " %8d ignored.\n",__func__,i);
        } else {
          pq1 = &mesh->quadra[i];
          pq1->tag[0] |= MG_PARBDY;
          pq1->tag[1] |= MG_PARBDY;
          pq1->tag[2] |= MG_PARBDY;
          pq1->tag[3] |= MG_PARBDY;
        }
      }
    }
  } //end if mesh->nquad

  /* read mesh edges */
  if ( mesh->na ) {
    na = mesh->na;

    rewind(inm);
    fseek(inm,posned,SEEK_SET);

    for (k=1; k<=na; k++) {
      pa = &mesh->edge[k];
      if (!bin) {
        MMG_FSCANF(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "",&pa->a,&pa->b,&pa->ref);
      }
      else {
        MMG_FREAD(&pa->a,MMG5_SW,1,inm);
        if(iswp) pa->a=MMG5_SWAPBIN(pa->a);
        MMG_FREAD(&pa->b,MMG5_SW,1,inm);
        if(iswp) pa->b=MMG5_SWAPBIN(pa->b);
        MMG_FREAD(&pa->ref,MMG5_SW,1,inm);
        if(iswp) pa->ref=MMG5_SWAPBIN(pa->ref);
      }
      pa->tag |= MG_REF;
      if ( pa->ref < 0 ) {
        pa->ref = -pa->ref;
        ++nref;
      }
    }

    /* get ridges */
    if ( nr ) {
      rewind(inm);
      fseek(inm,posnr,SEEK_SET);
      for (k=1; k<=nr; k++) {
        if(!bin) {
          MMG_FSCANF(inm,"%" MMG5_PRId "",&ia);
        }
        else {
          MMG_FREAD(&ia,MMG5_SW,1,inm);
          if(iswp) ia=MMG5_SWAPBIN(ia);
        }
        if(ia>na) {
          fprintf(stderr,"\n  ## Warning: %s: ridge number %8" MMG5_PRId " ignored.\n",
                  __func__,ia);
          continue;
        }
        pa = &mesh->edge[ia];
        pa->tag |= MG_GEO;
      }
    }
    /* get required edges */
    if ( nedreq ) {
      rewind(inm);
      fseek(inm,posnedreq,SEEK_SET);
      for (k=1; k<=nedreq; k++) {
        if(!bin) {
          MMG_FSCANF(inm,"%" MMG5_PRId "",&ia);
        }
        else {
          MMG_FREAD(&ia,MMG5_SW,1,inm);
          if(iswp) ia=MMG5_SWAPBIN(ia);
        }
        if(ia>na) {
          fprintf(stderr,"\n  ## Warning: %s: required Edges number %8" MMG5_PRId "/%8" MMG5_PRId ""
                  " ignored.\n",__func__,ia,na);
          continue;
        }
        pa = &mesh->edge[ia];
        pa->tag |= MG_REQ;
      }
    }
    /* get parallel edges */
    if ( nedpar ) {
      rewind(inm);
      fseek(inm,posnedpar,SEEK_SET);
      for (k=1; k<=nedpar; k++) {
        if(!bin) {
          MMG_FSCANF(inm,"%" MMG5_PRId "",&ia);
        }
        else {
          MMG_FREAD(&ia,MMG5_SW,1,inm);
          if(iswp) ia=MMG5_SWAPBIN(ia);
        }
        if(ia>na) {
          fprintf(stderr,"\n  ## Warning: %s: parallel Edges number %8" MMG5_PRId "/%8" MMG5_PRId ""
                  " ignored.\n",__func__,ia,na);
          continue;
        }
        pa = &mesh->edge[ia];
        pa->tag |= MG_PARBDY;
      }
    }
  }

  /* read mesh tetrahedra */
  rewind(inm);
  fseek(inm,posne,SEEK_SET);
  mesh->xt = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if (!bin) {
      MMG_FSCANF(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "",&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&ref);
    }
    else {
      for (i=0 ; i<4 ; i++) {
        MMG_FREAD(&pt->v[i],MMG5_SW,1,inm);
        if(iswp) pt->v[i]=MMG5_SWAPBIN(pt->v[i]);
      }
      MMG_FREAD(&ref,MMG5_SW,1,inm);
      if(iswp) ref=MMG5_SWAPBIN(ref);
    }
    if(ref < 0) {
      nref++;
    }
    pt->ref  = MMG5_abs(ref);//0;//ref ;
    for (i=0; i<4; i++) {
      ppt = &mesh->point[pt->v[i]];
      ppt->tag &= ~MG_NUL;
    }

    /* Possibly switch 2 vertices number so that each tet is positively oriented */
    if ( MMG5_orvol(mesh->point,pt->v) < 0.0 ) {
      /* mesh->xt temporary used to count reoriented tetra */
      mesh->xt++;
      aux = pt->v[2];
      pt->v[2] = pt->v[3];
      pt->v[3] = aux;
    }
  }
  if(mesh->xt) {
    fprintf(stderr,"\n     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n");
    fprintf(stderr,"         BAD ORIENTATION : vol < 0 -- %8" MMG5_PRId " tetra reoriented\n",mesh->xt);
    fprintf(stderr,"     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n");
  }
  mesh->xt = 0;
  /* get required tetrahedra */
  if(nereq) {
    rewind(inm);
    fseek(inm,posnereq,SEEK_SET);
    for (k=1; k<=nereq; k++) {
      if(!bin) {
        MMG_FSCANF(inm,"%d",&i);
      }
      else {
        MMG_FREAD(&i,MMG5_SW,1,inm);
        if(iswp) i=MMG5_swapbin(i);
      }
      if(i>mesh->ne) {
        fprintf(stderr,"\n  ## Warning: %s: required Tetra number %8d"
                " ignored.\n",__func__,i);
        continue;
      }
      pt = &mesh->tetra[i];
      pt->tag |= MG_REQ;
    }
  }
  /* get parallel tetrahedra */
  if(nepar) {
    rewind(inm);
    fseek(inm,posnepar,SEEK_SET);
    for (k=1; k<=nepar; k++) {
      if(!bin) {
        MMG_FSCANF(inm,"%d",&i);
      }
      else {
        MMG_FREAD(&i,MMG5_SW,1,inm);
        if(iswp) i=MMG5_swapbin(i);
      }
      if(i>mesh->ne) {
        fprintf(stderr,"\n  ## Warning: %s: parallel Tetra number %8d"
                " ignored.\n",__func__,i);
        continue;
      }
      pt = &mesh->tetra[i];
      pt->tag |= MG_PARBDY;
    }
  }
  /* read mesh prisms */
  rewind(inm);
  fseek(inm,posnprism,SEEK_SET);
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    if (!bin) {
      MMG_FSCANF(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "",&pp->v[0],&pp->v[1],&pp->v[2],
             &pp->v[3],&pp->v[4],&pp->v[5],&ref);
    }
    else {
      for (i=0 ; i<6 ; i++) {
        MMG_FREAD(&pp->v[i],MMG5_SW,1,inm);
        if(iswp) pp->v[i]=MMG5_SWAPBIN(pp->v[i]);
      }
      MMG_FREAD(&ref,MMG5_SW,1,inm);
      if(iswp) ref=MMG5_SWAPBIN(ref);
    }
    pp->ref  = ref;
    if ( pp-> ref < 0 ) {
      pp->ref = -pp->ref;
      ++nref;
    }
    for (i=0; i<6; i++) {
      ppt = &mesh->point[pp->v[i]];
      ppt->tag &= ~MG_NUL;
    }
  }

  if ( nref ) {
    fprintf(stdout,"\n     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n");
    fprintf(stdout,"         WARNING : %" MMG5_PRId " entities with unexpected refs (ref< 0).\n",nref);
    fprintf(stdout,"                   We take their absolute values.\n");
    fprintf(stdout,"     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n");
  }

  if ( !mesh->info.iso ) {
    /* read geometric entities */
    if ( mesh->nc1 && !ng ) {
      fprintf(stderr,"\n  ## Warning: %s: your mesh don't contains Normals but contains"
              " NormalAtVertices. The NormalAtVertices are deleted. \n",__func__);
      mesh->nc1 = 0;
    }

    if ( ng > 0 ) {
      MMG5_SAFE_CALLOC(norm,3*ng+1,double,return -1);

      rewind(inm);
      fseek(inm,posnormal,SEEK_SET);
      for (k=1; k<=ng; k++) {
        n = &norm[3*(k-1)+1];
        if ( mesh->ver == 1 ) {
          if (!bin) {
            for (i=0 ; i<3 ; i++) {
              MMG_FSCANF(inm,"%f",&fc);
              n[i] = (double) fc;
            }
          } else {
            for (i=0 ; i<3 ; i++) {
              MMG_FREAD(&fc,MMG5_SW,1,inm);
              if(iswp) fc=MMG5_swapf(fc);
              n[i] = (double) fc;
            }
          }
        }
        else {
          if (!bin) {
            MMG_FSCANF(inm,"%lf %lf %lf",&n[0],&n[1],&n[2]);
          }
          else {
            for (i=0 ; i<3 ; i++) {
              MMG_FREAD(&n[i],MMG5_SD,1,inm);
              if(iswp) n[i]=MMG5_swapd(n[i]);
            }
          }
        }
        dd = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if ( dd > MMG5_EPSD2 ) {
          dd = 1.0 / sqrt(dd);
          n[0] *= dd;
          n[1] *= dd;
          n[2] *= dd;
        }
      }

      rewind(inm);
      fseek(inm,posnc1,SEEK_SET);

      for (k=1; k<=mesh->nc1; k++) {
        if (!bin) {
          MMG_FSCANF(inm,"%" MMG5_PRId " %" MMG5_PRId "",&ip,&idn);
        }
        else {
          MMG_FREAD(&ip,MMG5_SW,1,inm);
          if(iswp) ip=MMG5_SWAPBIN(ip);
          MMG_FREAD(&idn,MMG5_SW,1,inm);
          if(iswp) idn=MMG5_SWAPBIN(idn);
        }
        if ( idn > 0 && ip < mesh->np+1 ) {
          if ( (mesh->info.iso ) &&  mesh->point[ip].xp == -1 ) {
            /* Do not store the normals at mesh->info.isoref points (ls mode) */
            continue;
          }
          memcpy(&mesh->point[ip].n,&norm[3*(idn-1)+1],3*sizeof(double));
        }
      }
      MMG5_SAFE_FREE(norm);
    }

    /* Delete the mark added iso mode */
    if ( (mesh->info.iso ) && mesh->nc1 ) {
      for (k=1; k<=mesh->np; k++) {
        mesh->point[k].xp = 0;
      }
    }
  }

  /* stats */
  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"     NUMBER OF VERTICES       %8" MMG5_PRId "\n",mesh->np);
    fprintf(stdout,"     NUMBER OF TETRAHEDRA     %8" MMG5_PRId "\n",mesh->ne);
    if ( mesh->nprism )
      fprintf(stdout,"     NUMBER OF PRISMS         %8" MMG5_PRId "\n",mesh->nprism);

    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES          %8" MMG5_PRId "\n",mesh->na);
      if ( nr )
        fprintf(stdout,"     NUMBER OF RIDGES         %8" MMG5_PRId "\n",nr);
    }
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES      %8" MMG5_PRId "\n",mesh->nt);
    if ( mesh->nquad )
      fprintf(stdout,"     NUMBER OF QUADRILATERALS %8" MMG5_PRId "\n",mesh->nquad);


    if ( npreq || nedreq || ntreq || nereq || nqreq ) {
      fprintf(stdout,"     NUMBER OF REQUIRED ENTITIES: \n");
      if ( npreq )
        fprintf(stdout,"               VERTICES       %8" MMG5_PRId " \n",npreq);
      if ( nedreq )
        fprintf(stdout,"               EDGES          %8" MMG5_PRId " \n",nedreq);
      if ( ntreq )
        fprintf(stdout,"               TRIANGLES      %8" MMG5_PRId " \n",ntreq);
      if ( nqreq )
        fprintf(stdout,"               QUADRILATERALS %8" MMG5_PRId " \n",nqreq);
      if ( nereq )
        fprintf(stdout,"               TETRAHEDRA     %8" MMG5_PRId " \n",nereq);
    }
    if(ncor)
      fprintf(stdout,"     NUMBER OF CORNERS        %8" MMG5_PRId " \n",ncor);

    if ( nppar || nedpar || ntpar || nepar || nqpar ) {
      fprintf(stdout,"     NUMBER OF PARALLEL ENTITIES: \n");
      if ( nppar )
        fprintf(stdout,"               VERTICES       %8" MMG5_PRId " \n",nppar);
      if ( nedpar )
        fprintf(stdout,"               EDGES          %8" MMG5_PRId " \n",nedpar);
      if ( ntpar )
        fprintf(stdout,"               TRIANGLES      %8" MMG5_PRId " \n",ntpar);
      if ( nqpar )
        fprintf(stdout,"               QUADRILATERALS %8" MMG5_PRId " \n",nqpar);
      if ( nepar )
        fprintf(stdout,"               TETRAHEDRA     %8" MMG5_PRId " \n",nepar);
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 *
 * \return 0 if the file is not found, -1 if we detect mismatch parameters or we
 * fail, 1 otherwise.
 *
 * Read mesh data.
 *
 */
int MMG3D_loadMesh(MMG5_pMesh mesh,const char *filename) {
  FILE*       inm;
  int         bin,ier;

  ier = MMG3D_openMesh(mesh->info.imprim,filename,&inm,&bin,"rb","rb");
  if( ier < 1 ) return ier;
  ier = MMG3D_loadMesh_opened(mesh,inm,bin);
  if( ier < 1 ) return ier;

  fclose(inm);
  return 1;
}


int MMG3D_loadMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename) {
  FILE*       inm;
  int         ier;
  long        posNodes,posElts,*posNodeData;
  int         bin,iswp;
  MMG5_int    nelts;
  int         nsols;

  mesh->dim = 3;

  ier = MMG5_loadMshMesh_part1(mesh,filename,&inm,
                               &posNodes,&posElts,&posNodeData,
                               &bin,&iswp,&nelts,&nsols);
  if ( ier < 1 ) return (ier);

  if ( nsols>1 ) {
    fprintf(stderr,"Error: SEVERAL SOLUTIONS FOUND (%d)\n",nsols);
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  if ( !MMG3D_zaldy(mesh) ) {
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  if (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne) {
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  if ( !mesh->ne ) {
    fprintf(stderr,"  ** MISSING DATA.\n");
    fprintf(stderr," Check that your mesh contains tetrahedra.\n");
    fprintf(stderr," Exit program.\n");
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  ier =  MMG5_loadMshMesh_part2( mesh, &sol,&inm,
                                 posNodes,posElts,posNodeData,
                                 bin,iswp,nelts,nsols);
  MMG5_SAFE_FREE(posNodeData);
  if ( ier < 1 ) {
    fprintf(stderr,"  ** ERROR WHEN PARSING THE INPUT FILE\n");
    return  ier;
  }

  /* Check the metric type */
  if ( sol ) {
    ier = MMG5_chkMetricType(mesh,&sol->type,&sol->entities,inm);
  }

  return ier;
}

int MMG3D_loadMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename) {
  FILE*       inm;
  int         ier;
  long        posNodes,posElts,*posNodeData;
  int         bin,iswp;
  MMG5_int    nelts;
  int         nsols;

  mesh->dim = 3;

  ier = MMG5_loadMshMesh_part1(mesh,filename,&inm,
                               &posNodes,&posElts,&posNodeData,
                               &bin,&iswp,&nelts,&nsols);
  if ( ier < 1 ) return (ier);

  mesh->nsols = nsols;
  if ( *sol )  MMG5_DEL_MEM(mesh,*sol);

  MMG5_ADD_MEM(mesh,nsols*sizeof(MMG5_Sol),"solutions array",
                printf("  Exit program.\n"); fclose(inm);
                MMG5_SAFE_FREE(posNodeData);
                return -1);
  MMG5_SAFE_CALLOC(*sol,nsols,MMG5_Sol,return -1);

  if ( !MMG3D_zaldy(mesh) ) {
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  if (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne) {
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  if ( !mesh->ne ) {
    fprintf(stderr,"  ** MISSING DATA.\n");
    fprintf(stderr," Check that your mesh contains tetrahedra.\n");
    fprintf(stderr," Exit program.\n");
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  ier =  MMG5_loadMshMesh_part2( mesh, sol,&inm,
                                 posNodes,posElts,posNodeData,
                                 bin,iswp,nelts,nsols);

  if ( ier < 1 ) {
    fprintf(stderr,"  ** ERROR WHEN PARSING THE INPUT FILE\n");
  }

  MMG5_SAFE_FREE(posNodeData);

  return ier;
}

int MMG3D_loadGenericMesh(MMG5_pMesh mesh, MMG5_pSol sol, const char *filename) {
  int ier=0;
  const char *filenameptr,*solnameptr;
  char *tmp,*soltmp;

  if ( filename && strlen(filename) ) {
    filenameptr = filename;
    solnameptr = filename;
  }
  else if (mesh->namein && strlen(mesh->namein) ) {
    filenameptr = mesh->namein;
    if ( sol && strlen(sol->namein) ) {
      solnameptr  = sol->namein;
    }
    else {
      solnameptr = mesh->namein;
    }
  }
  else {
    fprintf(stderr,"  ## Error: %s: please provide input file name"
            " (either in the mesh structure or as function argument).\n",
            __func__);
    return -1;
  }

  MMG5_SAFE_MALLOC(tmp,strlen(filenameptr)+1,char,return -1);
  strcpy(tmp,filenameptr);

  /* read mesh/sol files */
  char *ptr   = MMG5_Get_filenameExt(tmp);
  int  fmt = MMG5_Get_format(ptr,MMG5_FMT_MeditASCII);

  switch ( fmt ) {

  case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
    ier = MMG3D_loadMshMesh(mesh,sol,tmp);
    break;

  case ( MMG5_FMT_VtkVtu ):
    ier = MMG3D_loadVtuMesh(mesh,sol,tmp);
    break;

  case ( MMG5_FMT_VtkVtk ):
    ier = MMG3D_loadVtkMesh(mesh,sol,tmp);
    break;

  case ( MMG5_FMT_MeditASCII ): case ( MMG5_FMT_MeditBinary ):
    ier = MMG3D_loadMesh(mesh,tmp);
    if ( ier <  1 ) { break; }

    /* Facultative metric */
    if ( sol ) {
      MMG5_SAFE_MALLOC(soltmp,strlen(solnameptr)+1,char,return -1);
      strcpy(soltmp,solnameptr);

      if ( MMG3D_loadSol(mesh,sol,soltmp) == -1) {
        fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        ier = 0;
      }
      MMG5_SAFE_FREE(soltmp);
    }
    break;

  default:
    fprintf(stderr,"  ** I/O AT FORMAT %s NOT IMPLEMENTED.\n",MMG5_Get_formatName(fmt) );
    ier= -1;
  }

  MMG5_SAFE_FREE(tmp);

  return ier;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename pointer toward the name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 * \warning you must call the \a MMG3D_packMesh function before saving your
 * mesh.
 */
int MMG3D_saveMesh(MMG5_pMesh mesh, const char *filename) {
  FILE*         inm;
  MMG5_pPoint   ppt;
  MMG5_pTetra   pt;
  MMG5_pPrism   pp;
  MMG5_pTria    ptt;
  MMG5_pQuad    pq;
  MMG5_xPoint   *pxp;
  MMG5_int      k,na,nc,np,ne,nn,nr,nre,npar,nedreq,nedpar,ntreq,ntpar,nt,nereq,nepar;
  MMG5_int      npr,nq,nqreq,nqpar,bpos;
  int           bin,binch;
  char          chaine[MMG5_FILESTR_LGTH];
  static int8_t parWarn = 0;

  mesh->ver = 2;

  if ( !MMG3D_openMesh(mesh->info.imprim,filename,&inm,&bin,"w","wb") ) {
    return 0;
  }

  /*entete fichier*/
  binch=0; bpos=10;
  if(!bin) {
    strcpy(&chaine[0],"MeshVersionFormatted 2\n");
    fprintf(inm,"%s",chaine);
    strcpy(&chaine[0],"\n\nDimension 3\n");
    fprintf(inm,"%s ",chaine);
  } else {
    binch = 1; //MeshVersionFormatted
    fwrite(&binch,MMG5_SW,1,inm);
    binch = 2; //version
    fwrite(&binch,MMG5_SW,1,inm);
    binch = 3; //Dimension
    fwrite(&binch,MMG5_SW,1,inm);
    bpos = 20; //Pos
    fwrite(&bpos,MMG5_SW,1,inm);
    binch = 3;
    fwrite(&binch,MMG5_SW,1,inm);

  }
  /* vertices */
  np = nc = na = nr = nre = npar = 0;

  if ( !mesh->point ) {
    fprintf(stderr, "\n  ## Error: %s: points array not allocated.\n",
            __func__);
    fclose(inm);
    return 0;
  }


  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      ppt->tmp  = ++np;
      ppt->flag = 0;
      if ( ppt->tag & MG_CRN )  nc++;
      if ( ppt->tag & MG_REQ )  nre++;
      if ( ppt->tag & MG_PARBDY )  npar++;
    }
  }

  if(!bin) {
    strcpy(&chaine[0],"\n\nVertices\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%" MMG5_PRId "\n",np);
  } else {
    binch = 4; //Vertices
    fwrite(&binch,MMG5_SW,1,inm);
    bpos += (3+(1+3*mesh->ver)*np)*MMG5_SW; //NullPos
    fwrite(&bpos,MMG5_SW,1,inm);
    fwrite(&np,MMG5_SW,1,inm);
  }
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      if(!bin) {
        fprintf(inm,"%.15lg %.15lg %.15lg %" MMG5_PRId "\n",
                ppt->c[0],ppt->c[1],ppt->c[2],MMG5_abs(ppt->ref));
      } else {
        fwrite((unsigned char*)&ppt->c[0],MMG5_SD,1,inm);
        fwrite((unsigned char*)&ppt->c[1],MMG5_SD,1,inm);
        fwrite((unsigned char*)&ppt->c[2],MMG5_SD,1,inm);
        ppt->ref = MMG5_abs(ppt->ref);
        fwrite((unsigned char*)&ppt->ref,MMG5_SW,1,inm);
      }
    }
  }

  /* corners+required */
  if ( nc ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nCorners\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nc);
    } else {
      binch = 13; //
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3+nc)*MMG5_SW; //NullPos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&nc,MMG5_SW,1,inm);
    }

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && ppt->tag & MG_CRN ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId "\n",ppt->tmp);
        } else {
          fwrite(&ppt->tmp,MMG5_SW,1,inm);
        }
      }
    }
  }
  if ( nre ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nRequiredVertices\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nre);
    } else {
      binch = 15; //
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3+nre)*MMG5_SW; //NullPos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&nre,MMG5_SW,1,inm);
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && ppt->tag & MG_REQ ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId "\n",ppt->tmp);
        } else {
          fwrite(&ppt->tmp,MMG5_SW,1,inm);
        }
      }
    }
  }
  if ( npar ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nParallelVertices\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",npar);
    }
    else {
      if ( !parWarn ) {
        parWarn = 1;
        fprintf(stderr, "\n  ## Warning: %s: parallel entities can't be"
                " saved at binary format. Ignored.\n",
                __func__);

      }
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && ppt->tag & MG_PARBDY ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId "\n",ppt->tmp);
        }
      }
    }
  }

  /* tetrahedra */
  ne = nereq = nepar = 0;
  if ( mesh->tetra ) {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) {
        continue;
      }
      ne++;
      if ( pt->tag & MG_REQ ){
        nereq++;
      }
      if ( pt->tag & MG_PARBDY ){
        nepar++;
      }
    }
  }
  else {
    fprintf(stderr, "\n  ## Warning: %s: tetra array not allocated.\n",
            __func__);
  }

  if(!bin) {
    strcpy(&chaine[0],"\n\nTetrahedra\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%" MMG5_PRId "\n",ne);
  } else {
    binch = 8; //Tetra
    fwrite(&binch,MMG5_SW,1,inm);
    bpos += (3 + 5*ne)*MMG5_SW;//Pos
    fwrite(&bpos,MMG5_SW,1,inm);
    fwrite((unsigned char*)&ne,MMG5_SW,1,inm);
  }
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if ( !MG_EOK(pt) ) { continue; }

    /* Tag the tetra vertices to detect points belonging to prisms only (because
     * we don't know the normals/tangents at this points, thus we don't want to
     * save it). */
    mesh->point[pt->v[0]].flag = 1;
    mesh->point[pt->v[1]].flag = 1;
    mesh->point[pt->v[2]].flag = 1;
    mesh->point[pt->v[3]].flag = 1;

    if(!bin) {
      fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp
              ,mesh->point[pt->v[2]].tmp,mesh->point[pt->v[3]].tmp,pt->ref);
    } else {
      fwrite(&mesh->point[pt->v[0]].tmp,MMG5_SW,1,inm);
      fwrite(&mesh->point[pt->v[1]].tmp,MMG5_SW,1,inm);
      fwrite(&mesh->point[pt->v[2]].tmp,MMG5_SW,1,inm);
      fwrite(&mesh->point[pt->v[3]].tmp,MMG5_SW,1,inm);
      fwrite(&pt->ref,MMG5_SW,1,inm);
    }
  }

  if ( nereq ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nRequiredTetrahedra\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nereq);
    } else {
      binch = 12; //RequiredTetra
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3 + nereq)*MMG5_SW;//Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&nereq,MMG5_SW,1,inm);
    }
    ne = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;
      ne++;
      if ( pt->tag & MG_REQ ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId "\n",ne);
        } else {
          fwrite(&ne,MMG5_SW,1,inm);
        }
      }
    }
  }

  if ( nepar ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nParallelTetrahedra\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nepar);
    } else {
      if ( !parWarn ) {
        parWarn = 1;
        fprintf(stderr, "\n  ## Warning: %s: parallel entities can't be"
                " saved at binary format. Ignored.\n",
                __func__);

      }
    }
    ne = 0;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;
      ne++;
      if ( pt->tag & MG_PARBDY ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId "\n",ne);
        }
      }
    }
  }

  /* prisms */
  npr = 0;
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) ) continue;
    npr++;
  }

  if ( npr ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nPrisms\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",npr);
    } else {
      binch = 9; //Prism
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3 + 7*npr)*MMG5_SW;//Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite((unsigned char*)&npr,MMG5_SW,1,inm);
    }
    for (k=1; k<=mesh->nprism; k++) {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) ) continue;

      if(!bin) {
        fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n"
                ,mesh->point[pp->v[0]].tmp,mesh->point[pp->v[1]].tmp
                ,mesh->point[pp->v[2]].tmp,mesh->point[pp->v[3]].tmp
                ,mesh->point[pp->v[4]].tmp,mesh->point[pp->v[5]].tmp,pp->ref);
      } else {
        fwrite(&mesh->point[pp->v[0]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[pp->v[1]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[pp->v[2]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[pp->v[3]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[pp->v[4]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[pp->v[5]].tmp,MMG5_SW,1,inm);
        fwrite(&pp->ref,MMG5_SW,1,inm);
      }
    }
  }


  nn = nt = 0;
  if ( mesh->xp && mesh->xpoint ) {
    /* Count tangents and normals */
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || (!ppt->flag) || MG_SIN(ppt->tag) ) {
        /* We compute normals at some required points but not all of them in analysis
         * thus it don't works to save it at the end of the process */
        continue;
      }
      else if ( (ppt->tag & MG_BDY)
                && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM)) ) {
         /* Indeed non-manifold point along connected mesh have a unique normal if
          * pxp->nnor==0 but saving it gives a weird visualization and keeping the
          * normal info is useless. */
        nn++;
      }
      if ( MG_EDG_OR_NOM(ppt->tag) ) {
        nt++;
      }
    }

    /* write normals */
    if(!bin) {
      strcpy(&chaine[0],"\n\nNormals\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nn);
    } else {      binch = 60; //normals
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3+(3*mesh->ver)*nn)*MMG5_SW; //Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&nn,MMG5_SW,1,inm);
    }

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || (!ppt->flag) || MG_SIN(ppt->tag) ) {
        continue;
      }
      else if ( (ppt->tag & MG_BDY)
                && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM)) ) {
        pxp = &mesh->xpoint[ppt->xp];
        if(!bin) {
          fprintf(inm,"%.15lg %.15lg %.15lg \n",pxp->n1[0],pxp->n1[1],pxp->n1[2]);
        } else {
          fwrite((unsigned char*)&pxp->n1[0],MMG5_SD,1,inm);
          fwrite((unsigned char*)&pxp->n1[1],MMG5_SD,1,inm);
          fwrite((unsigned char*)&pxp->n1[2],MMG5_SD,1,inm);
        }
      }
    }

    if(!bin) {
      strcpy(&chaine[0],"\n\nNormalAtVertices\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nn);
    } else {
      binch = 20; //normalatvertices
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3 + 2*nn)*MMG5_SW;//Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&nn,MMG5_SW,1,inm);
    }
    nn = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) || (!ppt->flag) || MG_SIN(ppt->tag) ) {
        continue;
      }
      else if ( (ppt->tag & MG_BDY)
                && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM) ) ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId "\n",ppt->tmp,++nn);
        } else {
          fwrite(&ppt->tmp,MMG5_SW,1,inm);
          ++nn;
          fwrite(&nn,MMG5_SW,1,inm);
        }
      }
    }

    if ( nt ) {
      /* Write tangents */
      if(!bin) {
        strcpy(&chaine[0],"\n\nTangents\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",nt);
      } else {
        binch = 59; //tangent
        fwrite(&binch,MMG5_SW,1,inm);
        bpos += (3+(3*mesh->ver)*nt)*MMG5_SW; //Pos
        fwrite(&bpos,MMG5_SW,1,inm);
        fwrite(&nt,MMG5_SW,1,inm);
      }

      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || (!ppt->flag) || MG_SIN(ppt->tag) ) {
          continue;
        }
        else if ( MG_EDG_OR_NOM(ppt->tag) ) {
          if(!bin) {
            fprintf(inm,"%.15lg %.15lg %.15lg \n",ppt->n[0],ppt->n[1],ppt->n[2]);
          } else {
            fwrite((unsigned char*)&ppt->n[0],MMG5_SD,1,inm);
            fwrite((unsigned char*)&ppt->n[1],MMG5_SD,1,inm);
            fwrite((unsigned char*)&ppt->n[2],MMG5_SD,1,inm);
          }
        }
      }


      if(!bin) {
        strcpy(&chaine[0],"\n\nTangentAtVertices\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",nt);
      } else {
        binch = 61; //tangentatvertices
        fwrite(&binch,MMG5_SW,1,inm);
        bpos += (3 + 2*nt)*MMG5_SW;//Pos
        fwrite(&bpos,MMG5_SW,1,inm);
        fwrite(&nt,MMG5_SW,1,inm);
      }
      nt = 0;
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) || (!ppt->flag) || MG_SIN(ppt->tag) ) {
          continue;
        }
        else if ( MG_EDG_OR_NOM(ppt->tag) ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId "\n",ppt->tmp,++nt);
          } else {
            fwrite(&ppt->tmp,MMG5_SW,1,inm);
            ++nt;
            fwrite(&(nn),MMG5_SW,1,inm);
          }
        }
      }
    }
  }

  /* boundary mesh */
  /* tria + required tria */
  ntreq = ntpar = 0;

  if ( mesh->nt ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nTriangles\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",mesh->nt);
    } else {
      binch = 6; //Triangles
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3+4*mesh->nt)*MMG5_SW; //Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&mesh->nt,MMG5_SW,1,inm);
    }
    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];
      if ( ptt->tag[0] & MG_REQ && ptt->tag[1] & MG_REQ && ptt->tag[2] & MG_REQ ) {
        ntreq++;
      }
      if ( ptt->tag[0] & MG_PARBDY && ptt->tag[1] & MG_PARBDY && ptt->tag[2] & MG_PARBDY ) {
        ntpar++;
      }
      if(!bin) {
        fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",mesh->point[ptt->v[0]].tmp,mesh->point[ptt->v[1]].tmp
                ,mesh->point[ptt->v[2]].tmp,ptt->ref);
      } else {
        fwrite(&mesh->point[ptt->v[0]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[ptt->v[1]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[ptt->v[2]].tmp,MMG5_SW,1,inm);
        fwrite(&ptt->ref,MMG5_SW,1,inm);
      }
    }
    if ( ntreq ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nRequiredTriangles\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",ntreq);
      } else {
        binch = 17; //ReqTriangles
        fwrite(&binch,MMG5_SW,1,inm);
        bpos += (3+ntreq)*MMG5_SW; //Pos
        fwrite(&bpos,MMG5_SW,1,inm);
        fwrite(&ntreq,MMG5_SW,1,inm);
      }
      for (k=0; k<=mesh->nt; k++) {
        ptt = &mesh->tria[k];
        if ( (ptt->tag[0] & MG_REQ) && (ptt->tag[1] & MG_REQ)
             && ptt->tag[2] & MG_REQ ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",k);
          } else {
            fwrite(&k,MMG5_SW,1,inm);
          }
        }
      }
    }
    if ( ntpar ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nParallelTriangles\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",ntpar);
      } else {
        if ( !parWarn ) {
          parWarn = 1;
          fprintf(stderr, "\n  ## Warning: %s: parallel entities can't be"
                  " saved at binary format. Ignored.\n",
                  __func__);
        }
      }
      for (k=0; k<=mesh->nt; k++) {
        ptt = &mesh->tria[k];
        if ( (ptt->tag[0] & MG_PARBDY) && (ptt->tag[1] & MG_PARBDY)
             && ptt->tag[2] & MG_PARBDY ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",k);
          }
        }
      }
    }
  }

  /* quad + required quad */
  nq = nqreq = nqpar = 0;

  if ( mesh->nquad ) {

    for (k=1; k<=mesh->nquad; k++) {
      pq = &mesh->quadra[k];
      if ( !MG_EOK(pq) ) {
        continue;
      }
      nq++;
      if ( pq->tag[0] & MG_REQ && pq->tag[1] & MG_REQ &&
           pq->tag[2] & MG_REQ && pq->tag[3] & MG_REQ ) {
        nqreq++;
      }
      if ( pq->tag[0] & MG_PARBDY && pq->tag[1] & MG_PARBDY &&
           pq->tag[2] & MG_PARBDY && pq->tag[3] & MG_PARBDY ) {
        nqpar++;
      }
    }
  }

  if ( nq ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nQuadrilaterals\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nq);
    } else {
      binch = 7; //Quadrilaterals
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3+5*nq)*MMG5_SW; //Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&nq,MMG5_SW,1,inm);
    }
    for (k=1; k<=mesh->nquad; k++) {
      pq = &mesh->quadra[k];
      if ( !MG_EOK(pq) ) continue;

      if(!bin) {
        fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",mesh->point[pq->v[0]].tmp,
                mesh->point[pq->v[1]].tmp,mesh->point[pq->v[2]].tmp,
                mesh->point[pq->v[3]].tmp, pq->ref);
      } else {
        fwrite(&mesh->point[pq->v[0]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[pq->v[1]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[pq->v[2]].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[pq->v[3]].tmp,MMG5_SW,1,inm);
        fwrite(&pq->ref,MMG5_SW,1,inm);
      }
    }
    if ( nqreq ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nRequiredQuadrilaterals\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",nqreq);
      } else {
        binch = 18; //ReqQuad
        fwrite(&binch,MMG5_SW,1,inm);
        bpos += (3+nqreq)*MMG5_SW; //Pos
        fwrite(&bpos,MMG5_SW,1,inm);
        fwrite(&nqreq,MMG5_SW,1,inm);
      }
      for (k=0; k<=mesh->nquad; k++) {
        pq = &mesh->quadra[k];
        if ( (pq->tag[0] & MG_REQ) && (pq->tag[1] & MG_REQ)
             && pq->tag[2] & MG_REQ && pq->tag[3] & MG_REQ ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",k);
          } else {
            fwrite(&k,MMG5_SW,1,inm);
          }
        }
      }
    }
    if ( nqpar ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nParallelQuadrilaterals\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",nqpar);
      } else {
        if ( !parWarn ) {
          parWarn = 1;
          fprintf(stderr, "\n  ## Warning: %s: parallel entities can't be"
                  " saved at binary format. Ignored.\n",
                  __func__);
        }
      }
      for (k=0; k<=mesh->nquad; k++) {
        pq = &mesh->quadra[k];
        if ( (pq->tag[0] & MG_PARBDY) && (pq->tag[1] & MG_PARBDY) &&
             pq->tag[2] & MG_PARBDY && pq->tag[3] & MG_PARBDY ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",k);
          }
        }
      }
    }
  }

  nr = nedreq = nedpar = 0;
  if ( mesh->na ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nEdges\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",mesh->na);
    } else {
      binch = 5; //Edges
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3 + 3*mesh->na)*MMG5_SW;//Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&mesh->na,MMG5_SW,1,inm);
    }
    for (k=1; k<=mesh->na; k++) {
      if(!bin) {
        fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",mesh->point[mesh->edge[k].a].tmp,
                mesh->point[mesh->edge[k].b].tmp,mesh->edge[k].ref);
      } else {
        fwrite(&mesh->point[mesh->edge[k].a].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->point[mesh->edge[k].b].tmp,MMG5_SW,1,inm);
        fwrite(&mesh->edge[k].ref,MMG5_SW,1,inm);
      }
      if ( mesh->edge[k].tag & MG_GEO ) nr++;
      if ( mesh->edge[k].tag & MG_REQ ) nedreq++;
      if ( mesh->edge[k].tag & MG_PARBDY ) nedpar++;
    }

    if ( nr ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nRidges\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",nr);
      } else {
        binch = 14; //Ridges
        fwrite(&binch,MMG5_SW,1,inm);
        bpos += (3 + nr)*MMG5_SW;//Pos
        fwrite(&bpos,MMG5_SW,1,inm);
        fwrite(&nr,MMG5_SW,1,inm);
      }
      na = 0;
      for (k=1; k<=mesh->na; k++) {
        na++;
        if ( mesh->edge[k].tag & MG_GEO ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",na);
          } else {
            fwrite(&na,MMG5_SW,1,inm);
          }
        }
      }
    }

    if ( nedreq ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nRequiredEdges\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",nedreq);
      } else {
        binch = 16; //RequiredEdges
        fwrite(&binch,MMG5_SW,1,inm);
        bpos += (3 + nedreq)*MMG5_SW;//Pos
        fwrite(&bpos,MMG5_SW,1,inm);
        fwrite(&nedreq,MMG5_SW,1,inm);
      }
      na = 0;
      for (k=1; k<=mesh->na; k++) {
        na++;
        if (  mesh->edge[k].tag & MG_REQ ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",na);
          } else {
            fwrite(&na,MMG5_SW,1,inm);
          }
        }
      }
    }
    if ( nedpar ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nParallelEdges\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",nedpar);
      } else {
        if ( !parWarn ) {
          parWarn = 1;
          fprintf(stderr, "\n  ## Warning: %s: parallel entities can't be"
                  " saved at binary format. Ignored.\n",
                  __func__);
        }
      }
      na = 0;
      for (k=1; k<=mesh->na; k++) {
        na++;
        if (  mesh->edge[k].tag & MG_PARBDY ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",na);
          }
        }
      }
    }
  }


  if ( mesh->info.imprim > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES       %8" MMG5_PRId "   CORNERS %8" MMG5_PRId ""
            "   REQUIRED %8" MMG5_PRId "\n",np,nc,nre);
    fprintf(stdout,"     NUMBER OF TETRAHEDRA     %8" MMG5_PRId "   REQUIRED  %8" MMG5_PRId "\n",
            ne,nereq);
    if ( npr )
      fprintf(stdout,"     NUMBER OF PRISMS         %8" MMG5_PRId "\n",npr);

    if ( na )
      fprintf(stdout,"     NUMBER OF EDGES          %8" MMG5_PRId "   RIDGES  %8" MMG5_PRId ""
              "   REQUIRED  %8" MMG5_PRId "\n",na,nr,nedreq);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES      %8" MMG5_PRId "   REQUIRED  %8" MMG5_PRId "\n",
              mesh->nt, ntreq);
    if ( nq )
      fprintf(stdout,"     NUMBER OF QUADRILATERALS %8" MMG5_PRId "   REQUIRED  %8" MMG5_PRId "\n",
              nq,nqreq);

    if ( npar || nedpar || ntpar || nepar || nqpar ) {
      fprintf(stdout,"     NUMBER OF PARALLEL ENTITIES: \n");
      if ( npar )
        fprintf(stdout,"                  VERTICES       %8" MMG5_PRId " \n",npar);
      if ( nedpar )
        fprintf(stdout,"                  EDGES          %8" MMG5_PRId " \n",nedpar);
      if ( ntpar )
        fprintf(stdout,"                  TRIANGLES      %8" MMG5_PRId " \n",ntpar);
      if ( nqpar )
        fprintf(stdout,"                  QUADRILATERALS %8" MMG5_PRId " \n",nqpar);
      if ( nepar )
        fprintf(stdout,"                  TETRAHEDRA    %8" MMG5_PRId " \n",nepar);
    }
  }

  /*fin fichier*/
  if(!bin) {
    strcpy(&chaine[0],"\n\nEnd\n");
    fprintf(inm,"%s",chaine);
  } else {
    binch = 54; //End
    fwrite(&binch,MMG5_SW,1,inm);
    bpos += 2*MMG5_SW; //bpos + End key
    fwrite(&bpos,MMG5_SW,1,inm);
  }
  fclose(inm);
  return 1;
}

int MMG3D_saveGenericMesh(MMG5_pMesh mesh, MMG5_pSol sol, const char *filename) {
  int ier=0;
  const char *filenameptr,*solnameptr;
  char *tmp,*soltmp;

  if ( filename && strlen(filename) ) {
    filenameptr = filename;
    solnameptr = filename;
  }
  else if (mesh->namein && strlen(mesh->namein) ) {
    filenameptr = mesh->namein;
    if ( sol && strlen(sol->namein) ) {
      solnameptr  = sol->namein;
    }
    else {
      solnameptr = mesh->namein;
    }
  }
  else {
    fprintf(stderr,"  ## Error: %s: please provide input file name"
            " (either in the mesh structure or as function argument).\n",
            __func__);
    return 0;
  }

  MMG5_SAFE_MALLOC(tmp,strlen(filenameptr)+1,char,return 0);
  strcpy(tmp,filenameptr);

  /* read mesh/sol files */
  char *ptr   = MMG5_Get_filenameExt(tmp);
  int  fmt = MMG5_Get_format(ptr,MMG5_FMT_MeditASCII);

  int8_t savesolFile = 0;

  switch ( fmt ) {
  case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
    ier = MMG3D_saveMshMesh(mesh,sol,tmp);
    break;
  case ( MMG5_FMT_VtkVtu ):
    ier = MMG3D_saveVtuMesh(mesh,sol,tmp);
    break;
  case ( MMG5_FMT_VtkVtk ):
    ier = MMG3D_saveVtkMesh(mesh,sol,tmp);
    break;
  case ( MMG5_FMT_Tetgen ):
    ier = MMG3D_saveTetgenMesh(mesh,tmp);
    savesolFile = 1;
    break;
  default:
    ier = MMG3D_saveMesh(mesh,tmp);
    savesolFile = 1;
    break;
  }

  if ( ier && savesolFile ) {
    /* Medit or tetgen output: save the solution in a .sol file */
    if ( sol && sol->np ) {
      MMG5_SAFE_MALLOC(soltmp,strlen(solnameptr)+1,char,return 0);
      strcpy(soltmp,solnameptr);

      if ( MMG3D_saveSol(mesh,sol,soltmp) == -1) {
        fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        ier = 0;
      }
      MMG5_SAFE_FREE(soltmp);
    }
  }

  return ier;
}

int MMG3D_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename) {

  return MMG5_saveMshMesh(mesh,&sol,filename,1);
}

int MMG3D_saveMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename) {

  return MMG5_saveMshMesh(mesh,sol,filename,0);
}

int MMG3D_loadSol(MMG5_pMesh mesh,MMG5_pSol met, const char *filename) {
  FILE       *inm;
  long       posnp;
  int        iswp,ier,ver,bin,*type,nsols,dim;
  MMG5_int   k,np;

  /** Read the file header */
  ier =  MMG5_loadSolHeader(filename,3,&inm,&ver,&bin,&iswp,&np,&dim,&nsols,
                             &type,&posnp,mesh->info.imprim);
  if ( ier < 1 ) return ier;

  if ( nsols!=1 ) {
    fprintf(stderr,"Error: SEVERAL SOLUTIONS FOUND (%d)\n",nsols);
    fclose(inm);
    MMG5_SAFE_FREE(type);
    return -1;
  }

  if ( mesh->np != np ) {
    fprintf(stderr,"  ** MISMATCHES DATA: THE NUMBER OF VERTICES IN "
            "THE MESH (%" MMG5_PRId ") DIFFERS FROM THE NUMBER OF VERTICES IN "
            "THE SOLUTION (%" MMG5_PRId ") \n",mesh->np,np);
    fclose(inm);
    MMG5_SAFE_FREE(type);
    return -1;
  }

  /* #MMG5_loadSolHeader function reads only solution at vertices so we don't
      have to check the entites on which the metric applies */
  int entities = MMG5_Vertex;
  ier = MMG5_chkMetricType(mesh,type,&entities,inm);

  if ( ier < 1 ) {
    MMG5_SAFE_FREE(type);
    return ier;
  }

  /* Allocate and store the header informations for each solution */
  if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type[0]) ) {
    fclose(inm);
    MMG5_SAFE_FREE(type);
    return -1;
  }
  /* For binary file, we read the verson inside the file */
  if ( ver ) met->ver = ver;

  MMG5_SAFE_FREE(type);

  /* Read mesh solutions */
  rewind(inm);
  fseek(inm,posnp,SEEK_SET);

  if ( met->ver == 1 ) {
    /* Simple precision */
    for (k=1; k<=mesh->np; k++) {
      if ( MMG5_readFloatSol3D(met,inm,bin,iswp,k) < 0 ) return -1;
    }
  }
  else {
    /* Double precision */
    for (k=1; k<=mesh->np; k++) {
      if ( MMG5_readDoubleSol3D(met,inm,bin,iswp,k) < 0 ) return -1;
    }
  }

  fclose(inm);

  /* stats */
  MMG5_printMetStats(mesh,met);

  return 1;
}

int MMG3D_loadAllSols(MMG5_pMesh mesh,MMG5_pSol *sol, const char *filename) {
  MMG5_pSol   psl;
  FILE        *inm;
  long        posnp;
  int         iswp,ier,ver,bin,*type,nsols,dim;
  MMG5_int    j,k,np;
  char        data[16];
  static char mmgWarn = 0;

  /** Read the file header */
  ier =  MMG5_loadSolHeader(filename,3,&inm,&ver,&bin,&iswp,&np,&dim,&nsols,
                            &type,&posnp,mesh->info.imprim);
  if ( ier < 1 ) return ier;

  if ( mesh->np != np ) {
    fprintf(stderr,"  ** MISMATCHES DATA: THE NUMBER OF VERTICES IN "
            "THE MESH (%" MMG5_PRId ") DIFFERS FROM THE NUMBER OF VERTICES IN "
            "THE SOLUTION (%" MMG5_PRId ") \n",mesh->np,np);
    fclose(inm);
    MMG5_SAFE_FREE(type);
    return -1;
  }

  /** Sol tab allocation */
  mesh->nsols = nsols;

  if ( nsols > MMG5_NSOLS_MAX ) {
    fprintf(stderr,"\n  ## Error: %s: unexpected number of data (%d).\n",
            __func__,nsols);
    MMG5_SAFE_FREE(type);
    fclose(inm);
    return -1;
  }

  if ( *sol )  MMG5_DEL_MEM(mesh,*sol);

  MMG5_ADD_MEM(mesh,nsols*sizeof(MMG5_Sol),"solutions array",
                printf("  Exit program.\n"); fclose(inm);
                MMG5_SAFE_FREE(type);
                return -1);
  MMG5_SAFE_CALLOC(*sol,nsols,MMG5_Sol,return -1);

  for ( j=0; j<nsols; ++j ) {
    psl = *sol + j;

    /* Give an arbitrary name to the solution because the Medit format has non
     * name field */
    sprintf(data,"sol_%" MMG5_PRId "",j);
    if ( !MMG3D_Set_inputSolName(mesh,psl,data) ) {
      if ( !mmgWarn ) {
        mmgWarn = 1;
        fprintf(stderr,"\n  ## Warning: %s: unable to set solution name for"
                " at least 1 solution.\n",__func__);
      }
    }

    /* Allocate and store the header informations for each solution */
    if ( !MMG3D_Set_solSize(mesh,psl,MMG5_Vertex,mesh->np,type[j]) ) {
      MMG5_SAFE_FREE(type);
      fclose(inm);
      return -1;
    }

    /* For binary file, we read the verson inside the file */
    if ( ver ) psl->ver = ver;
  }
  MMG5_SAFE_FREE(type);

  /* read mesh solutions */
  rewind(inm);
  fseek(inm,posnp,SEEK_SET);

  if ( (*sol)[0].ver == 1 ) {
    /* Simple precision */
    for (k=1; k<=mesh->np; k++) {
      for ( j=0; j<nsols; ++j ) {
        psl = *sol + j;
        if ( MMG5_readFloatSol3D(psl,inm,bin,iswp,k) < 0 ) return -1;
      }
    }
  }
  else {
    /* Double precision */
    for (k=1; k<=mesh->np; k++) {
      for ( j=0; j<nsols; ++j ) {
        psl = *sol + j;
        if ( MMG5_readDoubleSol3D(psl,inm,bin,iswp,k) < 0 ) return -1;
      }
    }
  }
  fclose(inm);

  /* stats */
  MMG5_printSolStats(mesh,sol);

  return 1;
}

int MMG3D_saveSol(MMG5_pMesh mesh,MMG5_pSol met, const char *filename) {
  FILE*        inm;
  MMG5_pPoint  ppt;
  int          binch,bin,ier;
  MMG5_int     k,bpos;

  if ( !met->m ) {
    fprintf(stderr,"\n  ## Warning: %s: no metric data to save.\n",__func__);
    return 1;
  }

  met->ver = 2;

  bpos = 0;
  ier = MMG5_saveSolHeader( mesh,filename,&inm,met->ver,&bin,&bpos,mesh->np,
                            met->dim,1,&met->entities,&met->type,&met->size);

  if ( ier < 1 )  return ier;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;

    MMG5_writeDoubleSol3D(mesh,met,inm,bin,k,1);
    fprintf(inm,"\n");
  }

  /* End file */
  if(!bin) {
    fprintf(inm,"\n\nEnd\n");
  } else {
    binch = 54; //End
    fwrite(&binch,MMG5_SW,1,inm);
  }
  fclose(inm);
  return 1;
}

int MMG3D_saveAllSols(MMG5_pMesh mesh,MMG5_pSol *sol, const char *filename) {
  MMG5_pSol    psl;
  FILE*        inm;
  MMG5_pPoint  ppt;
  MMG5_pTetra  pt;
  int          binch,bin,ier,npointSols,ncellSols;
  MMG5_int     j,k,bpos;
  int          *type,*entities,*size;

  if ( !(*sol)[0].m )  return -1;

  (*sol)[0].ver = 2;

  MMG5_SAFE_CALLOC(type,mesh->nsols,int,return 0);
  MMG5_SAFE_CALLOC(size,mesh->nsols,int,MMG5_SAFE_FREE(type);return 0);
  MMG5_SAFE_CALLOC(entities,mesh->nsols,int,
                   MMG5_SAFE_FREE(type);MMG5_SAFE_FREE(size);return 0);

  npointSols = 0;
  ncellSols = 0;
  for (k=0; k<mesh->nsols; ++k ) {
    if ( ((*sol)[k].entities==MMG5_Noentity) || ((*sol)[k].entities==MMG5_Vertex) ) {
      ++npointSols;
    }
    else if ( (*sol)[k].entities == MMG5_Tetrahedron ) {
      ++ncellSols;
    }
    else {
      printf("\n  ## Warning: %s: unexpected entity type for solution %" MMG5_PRId ": %s."
             "\n Ignored.\n",
             __func__,k,MMG5_Get_entitiesName((*sol)[k].entities));
    }
    type[k]     = (*sol)[k].type;
    size[k]     = (*sol)[k].size;
    entities[k] = (*sol)[k].entities;
  }

  bpos = 0;
  ier = MMG5_saveSolHeader( mesh,filename,&inm,(*sol)[0].ver,&bin,&bpos,mesh->np,
                            (*sol)[0].dim,mesh->nsols,entities,type,size);

  if ( ier < 1 )  return ier;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;

    for ( j=0; j<mesh->nsols; ++j ) {
      psl = *sol+j;

      if ( (psl->entities==MMG5_Noentity) || (psl->entities==MMG5_Vertex) ) {
        MMG5_writeDoubleSol3D(mesh,psl,inm,bin,k,0);
      }
    }
    fprintf(inm,"\n");
  }

  MMG5_saveSolAtTetrahedraHeader( mesh,inm,(*sol)[0].ver,bin,&bpos,mesh->nsols,
                                  ncellSols,entities,type,size );

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for ( j=0; j<mesh->nsols; ++j ) {
      psl = *sol+j;
      if ( psl->entities==MMG5_Tetrahedron ) {
        MMG5_writeDoubleSol3D(mesh,psl,inm,bin,k,0);
      }
    }
    fprintf(inm,"\n");
  }


  MMG5_SAFE_FREE(type);
  MMG5_SAFE_FREE(size);
  MMG5_SAFE_FREE(entities);

  /* End file */
  if(!bin) {
    fprintf(inm,"\n\nEnd\n");
  } else {
    binch = 54; //End
    fwrite(&binch,MMG5_SW,1,inm);
  }
  fclose(inm);
  return 1;
}

static inline
int MMG3D_saveEle(MMG5_pMesh mesh,const char *filename) {
  FILE*             inm;
  MMG5_pTetra       pt;
  MMG5_int          k,i,ne;
  char              *ptr,*data;

  if ( !mesh->ne ) {
    return 1;
  }

  if ( (!filename) || !(*filename) ) {
    filename = mesh->nameout;
  }
  if ( (!filename) || !(*filename) ) {
    printf("\n  ## Error: %s: unable to save a file without a valid filename\n.",
           __func__);
    return 0;
  }

  /* Name of file */
  MMG5_SAFE_CALLOC(data,strlen(filename)+5,char,return 0);
  strcpy(data,filename);
  ptr = strstr(data,".node");
  if ( ptr ) {
    *ptr = '\0';
  }

  /* Add .node ext  */
  strcat(data,".ele");
  if( !(inm = fopen(data,"wb")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
    MMG5_SAFE_FREE(data);
    return 0;
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

  ne    = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    ne++;
  }

  /* Save elt number, node number per elt, 1 bdy marker */
  fprintf(inm, "%" MMG5_PRId " %d %d\n\n",ne,mesh->dim+1,1);

  ne = 0;
  for ( k=1; k<=mesh->ne; ++k ) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) {
      /* Save elt idx */
      fprintf(inm, "%" MMG5_PRId " ",++ne);

      /* Save connectivity */
      for ( i=0; i<=mesh->dim; ++i ) {
        fprintf(inm, "%" MMG5_PRId " ",mesh->point[pt->v[i]].tmp);
      }

      /* Save bdy marker */
      fprintf(inm, "%" MMG5_PRId "\n",pt->ref);
    }
  }
  fprintf(stdout,"     NUMBER OF ELEMENT       %8" MMG5_PRId "\n",ne);

  fclose(inm);

  return 1;
}

static inline
int MMG3D_saveNeigh(MMG5_pMesh mesh,const char *filename) {
  FILE*             inm;
  MMG5_pTetra       pt;
  MMG5_int          k,i,ne;
  MMG5_int          idx;
  char              *ptr,*data;

  if ( !mesh->na ) {
    return 1;
  }

  if ( (!filename) || !(*filename) ) {
    filename = mesh->nameout;
  }
  if ( (!filename) || !(*filename) ) {
    printf("\n  ## Error: %s: unable to save a file without a valid filename\n.",
           __func__);
    return 0;
  }

  /* Name of file */
  MMG5_SAFE_CALLOC(data,strlen(filename)+7,char,return 0);
  strcpy(data,filename);
  ptr = strstr(data,".node");
  if ( ptr ) {
    *ptr = '\0';
  }

  /* Add .node ext  */
  strcat(data,".neigh");
  if( !(inm = fopen(data,"wb")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
    MMG5_SAFE_FREE(data);
    return 0;
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

  if ( ! mesh->adja ) {
    if ( !MMG3D_hashTetra(mesh,1) ) {
      printf("\n  ## Error: %s: unable to compute triangle adjacencies\n.",__func__);
      return 0;
    }
  }

  ne    = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    ne++;
  }

  /* Save elt number, number of neighbors per elt */
  fprintf(inm, "%" MMG5_PRId " %d\n\n",ne,mesh->dim+1);

  ne = 0;
  for ( k=1; k<=mesh->ne; ++k ) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) ) {
      /* Save elt idx */
      fprintf(inm, "%" MMG5_PRId " ",++ne);

      /* Save neighbors */
      for ( i=1; i<=mesh->dim+1; ++i ) {
        /* The triangle conventions is that no neighbors <=> -1 */
        idx = ( mesh->adja[4*(k-1)+i] > 0 ) ? mesh->adja[4*(k-1)+i]/4 : -1;
        fprintf(inm, "%"MMG5_PRId" ",idx);
      }
      /* Save bdy marker */
      fprintf(inm, "\n");
    }
  }

  fclose(inm);

  return 1;
}

static inline
int MMG3D_saveFace(MMG5_pMesh mesh,const char *filename) {
  FILE*             inm;
  MMG5_pTria        pt;
  MMG5_int          k,i;
  char              *ptr,*data;

  if ( !mesh->nt ) {
    return 1;
  }

  if ( (!filename) || !(*filename) ) {
    filename = mesh->nameout;
  }
  if ( (!filename) || !(*filename) ) {
    printf("\n  ## Error: %s: unable to save a file without a valid filename\n.",
           __func__);
    return 0;
  }

  /* Name of file */
  MMG5_SAFE_CALLOC(data,strlen(filename)+6,char,return 0);
  strcpy(data,filename);
  ptr = strstr(data,".node");
  if ( ptr ) {
    *ptr = '\0';
  }

  /* Add .node ext  */
  strcat(data,".face");
  if( !(inm = fopen(data,"wb")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
    MMG5_SAFE_FREE(data);
    return 0;
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

  /* Save tria number, 1 bdy marker */
  fprintf(inm, "%" MMG5_PRId " %d\n\n",mesh->nt,1);

  for ( k=1; k<=mesh->nt; ++k ) {
    pt = &mesh->tria[k];
    if ( MG_EOK(pt) ) {
      /* Save elt idx */
      fprintf(inm, "%" MMG5_PRId " ",k);

      /* Save connectivity */
      for ( i=0; i<mesh->dim; ++i ) {
        fprintf(inm, "%" MMG5_PRId " ",mesh->point[pt->v[i]].tmp);
      }

      /* Save bdy marker */
      fprintf(inm, "%" MMG5_PRId "\n",pt->ref);
    }
  }
  fprintf(stdout,"     NUMBER OF TRIANGLES       %8" MMG5_PRId "\n",mesh->nt);

  fclose(inm);

  return 1;
}


int MMG3D_saveTetgenMesh(MMG5_pMesh mesh,const char *filename) {

  if ( !MMG5_saveNode(mesh,filename) ) {
    return 0;
  }

  if ( !MMG3D_saveEle(mesh,filename) ) {
    return 0;
  }

  if ( !MMG3D_saveFace(mesh,filename) ) {
    return 0;
  }

  if ( !MMG5_saveEdge(mesh,filename,".edge") ) {
    return 0;
  }

  if ( !MMG3D_saveNeigh(mesh,filename) ) {
    return 0;
  }

  return 1;
}
