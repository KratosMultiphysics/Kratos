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
 * \file mmgs/inout_s.c
 * \brief Input / Output Functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */


#include "libmmgs.h"
#include "libmmgs_private.h"
#include <math.h>

int MMGS_loadMesh(MMG5_pMesh mesh, const char *filename) {
  FILE        *inm;
  MMG5_pTria  pt1,pt2;
  MMG5_pPoint ppt;
  double      *norm,*n,dd;
  float       fc;
  long        posnp,posnt,posne,posncor,posnq,posned,posnr;
  long        posntreq,posnpreq,posnormal,posnc1;
  MMG5_int    k,ia,nq,nri,ip,idn,ng,npreq,nref;
  int         i;
  MMG5_int    ncor,nedreq,ntreq,posnedreq,bpos;
  int         binch,bin,iswp,bdim;
  MMG5_int    na;
  char        *ptr,*data;
  char        chaine[MMG5_FILESTR_LGTH],strskip[MMG5_FILESTR_LGTH];

  posnp = posnt = posne = posncor = posnq = 0;
  posned = posnr = posnpreq = posntreq = posnc1 = npreq = 0;
  posnedreq = posnormal = 0;
  ncor = nri = ng = nedreq = nq = ntreq = 0;
  bin = 0;
  iswp = 0;
  mesh->np = mesh->nt = mesh->nti = mesh->npi = 0;

  nref = 0;

  MMG5_SAFE_CALLOC(data,strlen(filename)+7,char,return -1);

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    /* data contains the filename without extension */
    strcat(data,".meshb");
    if( !(inm = fopen(data,"rb")) ) {
      /* our file is not a .meshb file, try with .mesh ext */
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if( !(inm = fopen(data,"rb")) ) {
        MMG5_SAFE_FREE(data);
        return 0;
      }
    }
    else  bin = 1;
  }
  else {
    ptr = strstr(data,".meshb");
    if ( ptr )  bin = 1;
    if( !(inm = fopen(data,"rb")) ) {
      MMG5_SAFE_FREE(data);
      return 0;
    }
  }

  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

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
      } else if(!strncmp(chaine,"Triangles",strlen("Triangles"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->nti);
        posnt = ftell(inm);
        continue;
      }
      else if(!strncmp(chaine,"RequiredTriangles",strlen("RequiredTriangles"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&ntreq);
        posntreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Quadrilaterals",strlen("Quadrilaterals"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nq);
        posnq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Corners",strlen("Corners"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&ncor);
        posncor = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Edges",strlen("Edges"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&mesh->na);
        posned = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredEdges",strlen("RequiredEdges"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nedreq);
        posnedreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Ridges",strlen("Ridges"))) {
        MMG_FSCANF(inm,"%" MMG5_PRId "",&nri);
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
      fprintf(stdout,"BAD FILE ENCODING\n");
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
      } else if(binch==7) {//Quadrilaterals
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&nq,MMG5_SW,1,inm);
        if(iswp) nq=MMG5_SWAPBIN(nq);
        posnq = ftell(inm);
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
      } else if(!mesh->na && binch==5) { //Edges
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&mesh->na,MMG5_SW,1,inm);
        if(iswp) mesh->na=MMG5_SWAPBIN(mesh->na);
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
      } else if(binch==14) {  //Ridges
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&nri,MMG5_SW,1,inm);
        if(iswp) nri=MMG5_SWAPBIN(nri);
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

  if ( !mesh->npi || !mesh->nti ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return -1;
  }
  mesh->np = mesh->npi;
  mesh->nt = mesh->nti + 2*nq;

  /* mem alloc */
  if ( !MMGS_zaldy(mesh) )  return -1;

  /* read vertices */

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
    ppt->tag = MG_NUL;
  }

  /* read triangles and set seed */
  if ( mesh->nti ) {
    rewind(inm);
    fseek(inm,posnt,SEEK_SET);
    for (k=1; k<=mesh->nti; k++) {
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
      for (i=0; i<3; i++) {
        ppt = &mesh->point[pt1->v[i]];
        ppt->tag &= ~MG_NUL;
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
        if ( i>mesh->nti ) {
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
  }

  /* read quads: automatic conversion into tria */
  if ( nq > 0 ) {
    rewind(inm);
    fseek(inm,posnq,SEEK_SET);

    printf("  ## Warning: %s: quadrangles automatically converted into"
           " triangles\n.",__func__);

    for (k=1; k<=nq; k++) {
      mesh->nti++;
      pt1 = &mesh->tria[mesh->nti];
      mesh->nti++;
      pt2 = &mesh->tria[mesh->nti];

      if (!bin) {
        MMG_FSCANF(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "",&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt2->v[2],&pt1->ref);
      }
      else {
        for (i=0 ; i<3 ; i++) {
          MMG_FREAD(&pt1->v[i],MMG5_SW,1,inm);
          if(iswp) pt1->v[i]=MMG5_SWAPBIN(pt1->v[i]);
        }
        MMG_FREAD(&pt2->v[2],MMG5_SW,1,inm);
        if(iswp) pt2->v[2]=MMG5_SWAPBIN(pt2->v[2]);
        MMG_FREAD(&pt1->ref,MMG5_SW,1,inm);
        if(iswp) pt1->ref=MMG5_SWAPBIN(pt1->ref);
      }
      if ( pt1->ref < 0 ) {
        pt1->ref = -pt1->ref;
        nref += 2;
      }

      pt2->v[0] = pt1->v[0];
      pt2->v[1] = pt1->v[2];
      pt2->ref  = pt1->ref;
      for (i=0; i<3; i++) {
        ppt = &mesh->point[pt1->v[i]];
        ppt->tag &= ~MG_NUL;
        ppt = &mesh->point[pt2->v[i]];
        ppt->tag &= ~MG_NUL;
      }
    }
    mesh->nt = mesh->nti;
  }

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

  /* read required vertices */
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
        fprintf(stderr,"\n  ## Warning: %s: required Vertices number %8d ignored\n",
                __func__,i);
      } else {
        ppt = &mesh->point[i];
        ppt->tag |= MG_REQ;
        ppt->tag &= ~MG_NUL;
      }
    }
  }

  /* Read mesh edges */
  na = mesh->na;
  if ( mesh->na ) {
    rewind(inm);
    fseek(inm,posned,SEEK_SET);

    for (k=1; k<=mesh->na; k++) {
      if (!bin) {
        MMG_FSCANF(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "",&mesh->edge[k].a,&mesh->edge[k].b,&mesh->edge[k].ref);
      }
      else {
        MMG_FREAD(&mesh->edge[k].a,MMG5_SW,1,inm);
        if(iswp) mesh->edge[k].a=MMG5_SWAPBIN(mesh->edge[k].a);
        MMG_FREAD(&mesh->edge[k].b,MMG5_SW,1,inm);
        if(iswp) mesh->edge[k].b=MMG5_SWAPBIN(mesh->edge[k].b);
        MMG_FREAD(&mesh->edge[k].ref,MMG5_SW,1,inm);
        if(iswp) mesh->edge[k].ref=MMG5_SWAPBIN(mesh->edge[k].ref);
      }
      if ( mesh->edge[k].ref < 0 ) {
        mesh->edge[k].ref = -mesh->edge[k].ref;
        ++nref;
      }
      mesh->edge[k].tag |= MG_REF;
      mesh->point[mesh->edge[k].a].tag |= MG_REF;
      mesh->point[mesh->edge[k].b].tag |= MG_REF;
    }

    /* get ridges */
    if ( nri ) {
      rewind(inm);
      fseek(inm,posnr,SEEK_SET);
      for (k=1; k<=nri; k++) {
        if(!bin) {
          MMG_FSCANF(inm,"%" MMG5_PRId "",&ia);
        }
        else {
          MMG_FREAD(&ia,MMG5_SW,1,inm);
          if(iswp) ia=MMG5_SWAPBIN(ia);
        }
        if ( (ia>na) || (ia<0) ) {
          fprintf(stderr,"\n  ## Warning: %s: ridge number %8" MMG5_PRId " ignored.\n",
                  __func__,ia);
        }
        else {
          mesh->edge[ia].tag |= MG_GEO;
        }
      }
    }

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
        if ( (ia>na) || (ia<0) ) {
          fprintf(stderr,"\n  ## Warning: %s: required edge number %8" MMG5_PRId " ignored\n",
                  __func__,ia);
        }
        else {
          mesh->edge[ia].tag |= MG_REQ;
        }
      }
    }
  }

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
      if ( idn > 0 && ip < mesh->np+1 )
        memcpy(&mesh->point[ip].n,&norm[3*(idn-1)+1],3*sizeof(double));
    }
    MMG5_SAFE_FREE(norm);
  }

  if ( nref ) {
    fprintf(stdout,"\n     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n");
    fprintf(stdout,"         WARNING : %" MMG5_PRId " entities with unexpected refs (ref< 0).\n",nref);
    fprintf(stdout,"                   We take their absolute values.\n");
    fprintf(stdout,"     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n");
  }

  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8" MMG5_PRId " / %8" MMG5_PRId "   CORNERS/REQ. %" MMG5_PRId " / %" MMG5_PRId "\n",mesh->npi,mesh->npmax,ncor,npreq);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8" MMG5_PRId " / %8" MMG5_PRId "\n",mesh->nti,mesh->ntmax);

    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8" MMG5_PRId "  RIDGES %6" MMG5_PRId "\n",mesh->na,nri);
  }
  fclose(inm);
  return 1;
}

int MMGS_loadMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename) {
  FILE*       inm;
  long        posNodes,posElts,*posNodeData;
  int         ier,bin,iswp,nsols;
  MMG5_int    nelts;

  mesh->dim = 3;

  ier = MMG5_loadMshMesh_part1(mesh,filename,&inm,
                               &posNodes,&posElts,&posNodeData,
                               &bin,&iswp,&nelts,&nsols);
  if ( ier < 1 )  return (ier);

  if ( nsols > 1 ) {
    fprintf(stderr,"Error: SEVERAL SOLUTIONS FOUND (%d)\n",nsols);
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  if ( !MMGS_zaldy(mesh) ) {
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  mesh->ne = mesh->nprism = 0;

  if ( !mesh->nt ) {
    fprintf(stderr,"  ** MISSING DATA.\n");
    fprintf(stderr," Check that your mesh contains triangles.\n");
    fprintf(stderr," Exit program.\n");
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }


  if (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt ) {
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  ier = MMG5_loadMshMesh_part2( mesh, &sol,&inm,
                                posNodes,posElts,posNodeData,
                                bin,iswp,nelts,nsols);
  MMG5_SAFE_FREE(posNodeData);
  if ( ier < 1 ) {
    fprintf(stderr,"  ** ERROR WHEN PARSING THE INPUT FILE\n");
    return  ier;
  }

  if ( sol ) {
    /* Check the metric type */
    ier = MMG5_chkMetricType(mesh,&sol->type,&sol->entities,inm);
  }

  return ier;
}

int MMGS_loadMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename) {
  FILE*       inm;
  long        posNodes,posElts,*posNodeData;
  int         ier,bin,iswp,nsols;
  MMG5_int    nelts;

  mesh->dim = 3;

  ier = MMG5_loadMshMesh_part1(mesh,filename,&inm,
                               &posNodes,&posElts,&posNodeData,
                               &bin,&iswp,&nelts,&nsols);
  if ( ier < 1 )  return (ier);

  mesh->nsols = nsols;
  if ( *sol )  MMG5_DEL_MEM(mesh,*sol);
  MMG5_ADD_MEM(mesh,nsols*sizeof(MMG5_Sol),"solutions array",
               printf("  Exit program.\n"); fclose(inm);
               MMG5_SAFE_FREE(posNodeData);
               return -1);
  MMG5_SAFE_CALLOC(*sol,nsols,MMG5_Sol,return -1);

  if ( !MMGS_zaldy(mesh) ) {
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  mesh->ne = mesh->nprism = 0;

  if ( !mesh->nt ) {
    fprintf(stderr,"  ** MISSING DATA.\n");
    fprintf(stderr," Check that your mesh contains triangles.\n");
    fprintf(stderr," Exit program.\n");
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }


  if (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt ) {
    fclose(inm);
    MMG5_SAFE_FREE(posNodeData);
    return -1;
  }

  ier = MMG5_loadMshMesh_part2( mesh, sol,&inm,
                                posNodes,posElts,posNodeData,
                                bin,iswp,nelts,nsols);

  if ( ier < 1 ) {
    fprintf(stderr,"  ** ERROR WHEN PARSING THE INPUT FILE\n");
  }

  MMG5_SAFE_FREE(posNodeData);

  return ier;
}

int MMGS_loadGenericMesh(MMG5_pMesh mesh, MMG5_pSol sol, const char *filename) {
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
  int  fmtin = MMG5_Get_format(ptr,MMG5_FMT_MeditASCII);

  switch ( fmtin ) {

  case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
    ier = MMGS_loadMshMesh(mesh,sol,tmp);
    break;

  case ( MMG5_FMT_VtkVtp ):
    ier = MMGS_loadVtpMesh(mesh,sol,tmp);
    break;

  case ( MMG5_FMT_VtkVtu ):
    ier = MMGS_loadVtuMesh(mesh,sol,tmp);
    break;

  case ( MMG5_FMT_VtkVtk ):
    ier = MMGS_loadVtkMesh(mesh,sol,tmp);
    break;

  case ( MMG5_FMT_MeditASCII ): case ( MMG5_FMT_MeditBinary ):
    ier = MMGS_loadMesh(mesh,tmp);
    if ( ier <  1 ) { break; }

    /* Facultative metric */
    if ( sol ) {
      MMG5_SAFE_MALLOC(soltmp,strlen(solnameptr)+1,char,return -1);
      strcpy(soltmp,solnameptr);

      if ( MMGS_loadSol(mesh,sol,tmp) == -1) {
        fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        ier = 0;
      }
      MMG5_SAFE_FREE(soltmp);
    }
    break;

  default:
    fprintf(stderr,"  ** I/O AT FORMAT %s NOT IMPLEMENTED.\n",MMG5_Get_formatName(fmtin) );
    ier= 0;
  }

  MMG5_SAFE_FREE(tmp);

  return ier;
}


int MMGS_saveMesh(MMG5_pMesh mesh, const char* filename) {
  FILE         *inm;
  MMG5_pPoint  ppt;
  MMG5_pTria   pt;
  MMG5_pxPoint go;
  MMG5_int     k,np,nt,nc,ng,nn,nr,nre,ntreq,bpos;
  int          bin,binch;
  // int          outm;
  char         *data,*ptr,chaine[MMG5_FILESTR_LGTH];

  mesh->ver = 2;

  bin = 0;

  MMG5_SAFE_CALLOC(data,strlen(filename)+7,char,return 0);

  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if( !(inm = fopen(data,"wb")) ) {
      ptr  = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if( !(inm = fopen(data,"w")) ) {
        MMG5_SAFE_FREE(data);
        return 0;
      }
    } else {
      bin = 1;
    }
  }
  else {
    ptr = strstr(data,".meshb");
    if( ptr ) {
      bin = 1;
      if( !(inm = fopen(data,"wb")) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        MMG5_SAFE_FREE(data);
        return 0;
      }
    } else {
      if( !(inm = fopen(data,"w")) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        MMG5_SAFE_FREE(data);
        return 0;
      }
    }
  }
  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

  /*entete fichier*/
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
  np = nc = ng = nn = nre = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    ppt->tmp = 0;
    if ( MG_VOK(ppt) ) {
      np++;
      ppt->tmp = np;
      if ( ppt->tag & MG_CRN )  nc++;
      if ( ppt->tag & MG_REQ )  nre++;
      if ( MG_EDG(ppt->tag) )   ng++;
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
        fprintf(inm,"%.15lg %.15lg %.15lg %" MMG5_PRId "\n",ppt->c[0],ppt->c[1],ppt->c[2],MMG5_abs(ppt->ref));
      } else {
        fwrite((unsigned char*)&ppt->c[0],MMG5_SD,1,inm);
        fwrite((unsigned char*)&ppt->c[1],MMG5_SD,1,inm);
        fwrite((unsigned char*)&ppt->c[2],MMG5_SD,1,inm);
        ppt->ref = MMG5_abs(ppt->ref);
        fwrite((unsigned char*)&ppt->ref,MMG5_SW,1,inm);
      }
      if ( !((ppt->tag & MG_GEO) || ppt->tag & MG_CRN) )  nn++;
    }
  }

  nt = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    nt++;
  }

  /* write corners */
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

  /* write edges, ridges */
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
    nre = nr = 0;
    for (k=1; k<=mesh->na; k++) {
      if(!bin) {
        fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",
                mesh->edge[k].a,mesh->edge[k].b,mesh->edge[k].ref);
      } else {
        fwrite(&mesh->edge[k].a,MMG5_SW,1,inm);
        fwrite(&mesh->edge[k].b,MMG5_SW,1,inm);
        fwrite(&mesh->edge[k].ref,MMG5_SW,1,inm);
      }
      if ( mesh->edge[k].tag & MG_REQ )  nre++;
      if ( mesh->edge[k].tag & MG_GEO )  nr++;
    }
    /* ridges */
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
      for (k=1; k<=mesh->na; k++) {
        if ( mesh->edge[k].tag & MG_GEO ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",k);
          } else {
            fwrite(&k,MMG5_SW,1,inm);
          }
        }
      }
    }
    if ( nre ) {
      if(!bin) {
        strcpy(&chaine[0],"\n\nRequiredEdges\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%" MMG5_PRId "\n",nre);
      } else {
        binch = 16; //RequiredEdges
        fwrite(&binch,MMG5_SW,1,inm);
        bpos += (3 + nre)*MMG5_SW;//Pos
        fwrite(&bpos,MMG5_SW,1,inm);
        fwrite(&nre,MMG5_SW,1,inm);
      }
      for (k=1; k<=mesh->na; k++)
        if ( mesh->edge[k].tag & MG_REQ )  {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",k);
          } else {
            fwrite(&k,MMG5_SW,1,inm);
          }
        }
    }
  }

  /* write triangles */
  if ( mesh->nt ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nTriangles\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nt);
    } else {
      binch = 6; //Triangles
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3+4*nt)*MMG5_SW; //Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&nt,MMG5_SW,1,inm);
    }

    ntreq = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( MG_EOK(pt) ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId " %" MMG5_PRId "\n",mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp
                  ,mesh->point[pt->v[2]].tmp,MMG5_abs(pt->ref));
        } else {
          fwrite(&mesh->point[pt->v[0]].tmp,MMG5_SW,1,inm);
          fwrite(&mesh->point[pt->v[1]].tmp,MMG5_SW,1,inm);
          fwrite(&mesh->point[pt->v[2]].tmp,MMG5_SW,1,inm);
          pt->ref = MMG5_abs(pt->ref);
          fwrite(&pt->ref,MMG5_SW,1,inm);
        }

        if ( (pt->tag[0] & MG_REQ) && (pt->tag[1] & MG_REQ) && (pt->tag[2] & MG_REQ) ) {
          ntreq++;
        }
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
        pt = &mesh->tria[k];
        if ( (pt->tag[0] & MG_REQ) && (pt->tag[1] & MG_REQ) && pt->tag[2] & MG_REQ ) {
          if(!bin) {
            fprintf(inm,"%" MMG5_PRId "\n",k);
          } else {
            fwrite(&k,MMG5_SW,1,inm);
          }
        }
      }
    }
  }

  if ( (!mesh->xp) || (!mesh->xpoint) ) nn = ng = 0;

  /* write normals */
  if ( nn ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nNormals\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",nn);
    } else {
      binch = 60; //normals
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3+(3*mesh->ver)*nn)*MMG5_SW; //Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&nn,MMG5_SW,1,inm);
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )  continue;
      else if ( !(ppt->tag & MG_GEO) && !(ppt->tag & MG_CRN) ) {
        if ( ppt->tag & MG_REF ) {
          go = &mesh->xpoint[ppt->xp];
          if(!bin) {
            fprintf(inm,"%.15lg %.15lg %.15lg \n",go->n1[0],go->n1[1],go->n1[2]);
          } else {
            fwrite((unsigned char*)&go->n1[0],MMG5_SD,1,inm);
            fwrite((unsigned char*)&go->n1[1],MMG5_SD,1,inm);
            fwrite((unsigned char*)&go->n1[2],MMG5_SD,1,inm);
          }
        }
        else {
          if(!bin) {
            fprintf(inm,"%.15lg %.15lg %.15lg \n",ppt->n[0],ppt->n[1],ppt->n[2]);
          } else {
            fwrite((unsigned char*)&ppt->n[0],MMG5_SD,1,inm);
            fwrite((unsigned char*)&ppt->n[1],MMG5_SD,1,inm);
            fwrite((unsigned char*)&ppt->n[2],MMG5_SD,1,inm);
          }
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
      if ( MG_VOK(ppt) && !((ppt->tag & MG_GEO) || (ppt->tag & MG_CRN) ) ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId "\n",ppt->tmp,++nn);
        } else {
          fwrite(&ppt->tmp,MMG5_SW,1,inm);
          ++nn;
          fwrite(&nn,MMG5_SW,1,inm);
        }
      }
    }
/*
  nn = 0;
  for (k=1; k<=mesh->nt; k++) {
  pt = &mesh->tria[k];
  if ( !MG_EOK(pt) )  continue;
  for (i=0; i<3; i++) {
  ppt = &mesh->point[pt->v[i]];
  if ( ppt->tag & MG_GEO )  nn++;
  }
  }
  GmfSetKwd(outm,GmfNormalAtTriangleVertices,nn);
  for (k=1; k<=mesh->nt; k++) {
  pt = &mesh->tria[k];
  if ( !MG_EOK(pt) )  continue;
  for (i=0; i<3; i++) {
  ppt = &mesh->point[pt->v[i]];
  if ( ppt->tag & MG_GEO )  nn++;
  }
*/
  }

  /* write tangents (only if we have already analyzed the mesh =>
   * mesh->xpoint is allocated ) */
  if ( ng && mesh->xpoint ) {
    /* Write tangents */
    if(!bin) {
      strcpy(&chaine[0],"\n\nTangents\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%" MMG5_PRId "\n",ng);
    } else {
      binch = 59; //tangent
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3+(3*mesh->ver)*ng)*MMG5_SW; //Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&ng,MMG5_SW,1,inm);
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && MG_EDG(ppt->tag) ) {
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
      fprintf(inm,"%" MMG5_PRId "\n",ng);
    } else {
      binch = 61; //tangentatvertices
      fwrite(&binch,MMG5_SW,1,inm);
      bpos += (3 + 2*ng)*MMG5_SW;//Pos
      fwrite(&bpos,MMG5_SW,1,inm);
      fwrite(&ng,MMG5_SW,1,inm);
    }
    ng = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( MG_VOK(ppt) && MG_EDG(ppt->tag) ) {
        if(!bin) {
          fprintf(inm,"%" MMG5_PRId " %" MMG5_PRId "\n",ppt->tmp,++ng);
        } else {
          fwrite(&ppt->tmp,MMG5_SW,1,inm);
          ++ng;
          fwrite(&(ng),MMG5_SW,1,inm);
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 4 ) {
    fprintf(stdout,"     NUMBER OF VERTICES   %8" MMG5_PRId "  CORNERS    %6" MMG5_PRId "\n",np,nc);
    fprintf(stdout,"     NUMBER OF TRIANGLES  %8" MMG5_PRId "\n",nt);

    if ( mesh->na )
      fprintf(stdout,"     NUMBER OF EDGES      %8" MMG5_PRId "  RIDGES     %6" MMG5_PRId "\n",mesh->na,nr);
    if ( nn+ng )
      fprintf(stdout,"     NUMBER OF NORMALS    %8" MMG5_PRId "  TANGENTS   %6" MMG5_PRId "\n",nn,ng);
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

int MMGS_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename) {

  return MMG5_saveMshMesh(mesh,&sol,filename,1);
}

int MMGS_saveMshMesh_and_allData(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename) {

  return MMG5_saveMshMesh(mesh,sol,filename,0);
}

int MMGS_loadSol(MMG5_pMesh mesh,MMG5_pSol met,const char* filename) {

  FILE       *inm;
  long        posnp;
  int         iswp,ier,*type,ver,bin,nsols,dim;
  MMG5_int    k,np;

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
  if ( !MMGS_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,type[0]) ) {
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

  /* isotropic metric */
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

int MMGS_loadAllSols(MMG5_pMesh mesh,MMG5_pSol *sol, const char *filename) {
  MMG5_pSol   psl;
  FILE        *inm;
  long        posnp;
  int         iswp,ier,*type,ver,bin,nsols,dim;
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
    if ( !MMGS_Set_inputSolName(mesh,psl,data) ) {
      if ( !mmgWarn ) {
        mmgWarn = 1;
        fprintf(stderr,"\n  ## Warning: %s: unable to set solution name for"
                " at least 1 solution.\n",__func__);
      }
    }

    /* Allocate and store the header informations for each solution */
    if ( !MMGS_Set_solSize(mesh,psl,MMG5_Vertex,mesh->np,type[j]) ) {
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

int MMGS_saveSol(MMG5_pMesh mesh,MMG5_pSol met, const char *filename) {
  FILE*        inm;
  MMG5_pPoint  ppt;
  int          binch,bin,ier;
  MMG5_int     bpos,k;

  if ( !met->m ) {
    fprintf(stderr,"\n  ## Warning: %s: no metric data to save.\n",__func__);
    return 1;
  }

  met->ver = 2;
  bpos = 0;
  ier = MMG5_saveSolHeader( mesh,filename,&inm,met->ver,&bin,&bpos,mesh->np,met->dim,
                            1,&met->entities,&met->type,&met->size);

  if ( ier < 1 )  return ier;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;

    MMG5_writeDoubleSol3D(mesh,met,inm,bin,k,1);
    fprintf(inm,"\n");
  }

  /*fin fichier*/
  if(!bin) {
    fprintf(inm,"\n\nEnd\n");
  } else {
    binch = 54; //End
    fwrite(&binch,MMG5_SW,1,inm);
  }
  fclose(inm);
  return 1;
}

int MMGS_saveAllSols(MMG5_pMesh mesh,MMG5_pSol *sol, const char *filename) {
  MMG5_pSol    psl;
  FILE*        inm;
  MMG5_pPoint  ppt;
  int          binch,bin,ier,npointSols,ncellSols;
  int          *type,*entities,j,*size;
  MMG5_int     k,bpos;

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
    else if ( (*sol)[k].entities == MMG5_Triangle ) {
      ++ncellSols;
    }
    else {
      printf("\n  ## Warning: %s: unexpected entity type for solution %" MMG5_PRId ": %s."
             "\n Ignored.\n",
             __func__,k,MMG5_Get_entitiesName((*sol)[k].entities));
    }

    type[k] = (*sol)[k].type;
    size[k] = (*sol)[k].size;
    entities[k] = (*sol)[k].entities;
  }

  bpos = 0;
  ier = MMG5_saveSolHeader( mesh,filename,&inm,(*sol)[0].ver,&bin,&bpos,mesh->np,
                            (*sol)[0].dim,mesh->nsols,entities,type,size);

  if ( ier < 1 )  return ier;

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) )  continue;

    for ( j=0; j<mesh->nsols; ++j ) {
      psl = *sol + j;

      if ( (psl->entities==MMG5_Noentity) || (psl->entities==MMG5_Vertex) ) {
        MMG5_writeDoubleSol3D(mesh,psl,inm,bin,k,0);
      }
    }
    fprintf(inm,"\n");
  }

  MMG5_saveSolAtTrianglesHeader( mesh,inm,(*sol)[0].ver,bin,&bpos,mesh->nsols,
                                 ncellSols,entities,type,size );

  for (k=1; k<=mesh->nt; k++) {
    MMG5_pTria ptt = &mesh->tria[k];
    if ( !MG_EOK(ptt) )  continue;

    for ( j=0; j<mesh->nsols; ++j ) {
      psl = *sol + j;
      if ( psl->entities==MMG5_Triangle ) {
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

int MMGS_saveGenericMesh(MMG5_pMesh mesh, MMG5_pSol sol, const char *filename) {
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
    ier = MMGS_saveMshMesh(mesh,sol,tmp);
    break;
  case ( MMG5_FMT_VtkVtp ):
    ier = MMGS_saveVtpMesh(mesh,sol,tmp);
    break;
  case ( MMG5_FMT_VtkVtu ):
    ier = MMGS_saveVtuMesh(mesh,sol,tmp);
    break;
  case ( MMG5_FMT_VtkVtk ):
    ier = MMGS_saveVtkMesh(mesh,sol,tmp);
    break;
  default:
    ier = MMGS_saveMesh(mesh,tmp);
    savesolFile = 1;
    break;
  }

  if ( ier && savesolFile ) {
    /* Medit or tetgen output: save the solution in a .sol file */
    if ( sol && sol->np ) {
      MMG5_SAFE_MALLOC(soltmp,strlen(solnameptr)+1,char,return 0);
      strcpy(soltmp,solnameptr);

      if ( MMGS_saveSol(mesh,sol,soltmp) == -1) {
        fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        ier = 0;
      }
      MMG5_SAFE_FREE(soltmp);
    }
  }

  return ier;
}
