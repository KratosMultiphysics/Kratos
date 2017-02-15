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
#include "mmg2d.h"

#define sw 4
#define sd 8

int MMG_swapbin(int sbin)
{
  int inv;
  char *p_in = (char *) &sbin;
  char *p = (char *)&inv;


  p[0] = p_in[3];
  p[1] = p_in[2];
  p[2] = p_in[1];
  p[3] = p_in[0];

  return(inv);
  /*unsigned char c1, c2, c3, c4;

    c1 = sbin & 255;
    c2 = (sbin >> 8) & 255;
    c3 = (sbin >> 16) & 255;
    c4 = (sbin >> 24) & 255;

    return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;   */

}
float MMG_swapf(float sbin)
{
  float out;
  char *p_in = (char *) &sbin;
  char *p_out = (char *) &out;
  p_out[0] = p_in[3];
  p_out[1] = p_in[2];
  p_out[2] = p_in[1];
  p_out[3] = p_in[0];

  return(out);
}
double MMG_swapd(double sbin)
{
  float out;
  char *p_in = (char *) &sbin;
  char *p_out = (char *) &out;
  int i;

  for(i=0;i<8;i++)
  {
    p_out[i] = p_in[7-i];
  }
  //printf("CONVERTION DOUBLE\n");
  return(out);
}

/* read mesh data */
int MMG2D_loadMesh(MMG5_pMesh mesh,const char *filename) {
  FILE        *inm;
  MMG5_pPoint       ppt;
  MMG5_pEdge        ped;
  MMG5_pTria        pt;
  float             fc;
  long         posnp,posnt,posncor,posned,posnq,posreq,posreqed;
  int          k,ref,tmp,ncor,norient,nreq,nreqed,bin,iswp,nq;
  char        *ptr,data[128],chaine[128];
  double       air,dtmp;
  int          i,bdim,binch,bpos;


  posnp = posnt = posncor = posned = posnq = posreq = posreqed = 0;
  ncor = nreq = nreqed = 0;
  bin = 0;
  iswp = 0;
  mesh->np = mesh->nt = mesh->na = mesh->xp = 0;
  nq = 0;


  strcpy(data,filename);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if (!(inm = fopen(data,"rb")) ) {
      ptr  = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if (!(inm = fopen(data,"rb")) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    }
    else  bin = 1;
  }
  else {
    ptr = strstr(data,".meshb");

    if ( ptr )  bin = 1;

    if( !(inm = fopen(data,"rb")) ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
      return(0);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if (!bin) {
    strcpy(chaine,"D");
    while(fscanf(inm,"%s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if(!strncmp(chaine,"MeshVersionFormatted",strlen("MeshVersionFormatted"))) {
        fscanf(inm,"%d",&mesh->ver);
        continue;
      } else if(!strncmp(chaine,"Dimension",strlen("Dimension"))) {
        fscanf(inm,"%d",&mesh->dim);
        if(mesh->info.nreg==2) {
          if(mesh->dim!=3) {
            fprintf(stdout,"WRONG USE OF -msh \n");
            return(0);
          }
          mesh->dim = 2;
        }
        if(mesh->dim!=2) {
          fprintf(stdout,"BAD DIMENSION : %d\n",mesh->dim);
          return(0);
        }
        continue;
      } else if(!strncmp(chaine,"Vertices",strlen("Vertices"))) {
        fscanf(inm,"%d",&mesh->np);
        posnp = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Triangles",strlen("Triangles"))) {
        fscanf(inm,"%d",&mesh->nt);
        posnt = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Corners",strlen("Corners"))) {
        fscanf(inm,"%d",&ncor);
        posncor = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredVertices",strlen("RequiredVertices"))) {
        fscanf(inm,"%d",&nreq);
        posreq = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Edges",strlen("Edges"))) {
        fscanf(inm,"%d",&mesh->na);
        posned = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"RequiredEdges",strlen("RequiredEdges"))) {
        fscanf(inm,"%d",&nreqed);
        posreqed = ftell(inm);
        continue;
      } else if(!strncmp(chaine,"Quadrilaterals",strlen("Quadrilaterals"))) {
        fscanf(inm,"%d",&nq);
        posnq = ftell(inm);
        continue;
      }
    }
  } else {
    bdim = 0;
    fread(&mesh->ver,sw,1,inm);
    iswp=0;
    if(mesh->ver==16777216)
      iswp=1;
    else if(mesh->ver!=1) {
      fprintf(stdout,"BAD FILE ENCODING\n");
    }
    fread(&mesh->ver,sw,1,inm);
    if(iswp) mesh->ver = MMG_swapbin(mesh->ver);
    while(fread(&binch,sw,1,inm)!=0 && binch!=54 ) {
      if(iswp) binch=MMG_swapbin(binch);
      if(binch==54) break;
      if(!bdim && binch==3) {  //Dimension
        fread(&bdim,sw,1,inm);  //NulPos=>20
        if(iswp) bdim=MMG_swapbin(bdim);
        fread(&bdim,sw,1,inm);
        if(iswp) bdim=MMG_swapbin(bdim);
        mesh->dim = bdim;
        if(bdim!=2) {
          fprintf(stdout,"BAD MESH DIMENSION : %d\n",mesh->dim);
          exit(EXIT_FAILURE);
          return(1);
        }
        continue;
      } else if(!mesh->np && binch==4) {  //Vertices
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=MMG_swapbin(bpos);
        fread(&mesh->np,sw,1,inm);
        if(iswp) mesh->np=MMG_swapbin(mesh->np);
        posnp = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      }  else if(!mesh->nt && binch==6) {//MMG5_Triangles
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=MMG_swapbin(bpos);
        fread(&mesh->nt,sw,1,inm);
        if(iswp) mesh->nt=MMG_swapbin(mesh->nt);
        posnt = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!ncor && binch==13) {
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=MMG_swapbin(bpos);
        fread(&ncor,sw,1,inm);
        if(iswp) ncor=MMG_swapbin(ncor);
        posncor = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!mesh->na && binch==5) { //Edges
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=MMG_swapbin(bpos);
        fread(&mesh->na,sw,1,inm);
        if(iswp) mesh->na=MMG_swapbin(mesh->na);
        posned = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!nreqed && binch==16) { //RequiredEdges
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=MMG_swapbin(bpos);
        fread(&nreqed,sw,1,inm);
        if(iswp) nreqed=MMG_swapbin(nreqed);
        posreqed = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if(!nreq && binch==15) { //RequiredVertices
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=MMG_swapbin(bpos);
        fread(&nreq,sw,1,inm);
        if(iswp) nreq=MMG_swapbin(nreq);
        posreq = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else {
        //printf("on traite ? %d\n",binch);
        fread(&bpos,sw,1,inm); //NulPos
        if(iswp) bpos=MMG_swapbin(bpos);
        //printf("on avance... Nulpos %d\n",bpos);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }

  }

  if ( abs(mesh->info.imprim) > 5 )
    fprintf(stdout,"  -- READING DATA FILE %s\n",data);

  if ( !mesh->np  ) {
    fprintf(stdout,"  ** MISSING DATA : no point\n");
    return(0);
  }
  if (!mesh->nt) {
    fprintf(stdout,"  **WARNING NO GIVEN TRIANGLE\n");
  }

  mesh->npi  = mesh->np;
  mesh->nai  = mesh->na;
  mesh->nti  = mesh->nt;
  if ( !mesh->np ) {
    fprintf(stdout,"  ** MISSING DATA\n");
    return(0);
  }

  /* mem alloc */
  if ( !MMG2D_zaldy(mesh) )  return(0);

  /* read vertices */
  rewind(inm);
  fseek(inm,posnp,SEEK_SET);
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if (mesh->ver < 2) { /*float*/
      if (!bin) {
        if(mesh->info.nreg==2) {
          for (i=0 ; i<3 ; i++) {
            fscanf(inm,"%f",&fc);
            if(i==2) break;
            ppt->c[i] = (double) fc;
          }
        } else {
          for (i=0 ; i<2 ; i++) {
            fscanf(inm,"%f",&fc);
            ppt->c[i] = (double) fc;
          }
        }
        fscanf(inm,"%d",&ppt->ref);
      } else {
        if(mesh->info.nreg==2) {
          fprintf(stdout,"binary not available with -msh option\n");
          return(0);
        }
        for (i=0 ; i<2 ; i++) {
          fread(&fc,sw,1,inm);
          if(iswp) fc=MMG_swapf(fc);
          ppt->c[i] = (double) fc;
        }
        fread(&ppt->ref,sw,1,inm);
        if(iswp) ppt->ref=MMG_swapbin(ppt->ref);
      }
    } else {
      if (!bin) {
        if(mesh->info.nreg==2) {
          fscanf(inm,"%lf %lf %lf %d",&ppt->c[0],&ppt->c[1],&dtmp,&ppt->ref);
        } else {
          fscanf(inm,"%lf %lf %d",&ppt->c[0],&ppt->c[1],&ppt->ref);
        }
      }
      else {
        for (i=0 ; i<2 ; i++) {
          fread(&ppt->c[i],sd,1,inm);
          if(iswp) ppt->c[i]=MMG_swapd(ppt->c[i]);
        }
        fread(&ppt->ref,sw,1,inm);
        if(iswp) ppt->ref=MMG_swapbin(ppt->ref);
      }
    }
    ppt->tag = 0;
    ppt->tag = MG_NUL;
  }

  /* read edges */
  rewind(inm);
  fseek(inm,posned,SEEK_SET);
  for (k=1; k<=mesh->na; k++) {
    ped = &mesh->edge[k];
    if (!bin)
      fscanf(inm,"%d %d %d",&ped->a,&ped->b,&ped->ref);
    else {
      fread(&ped->a,sw,1,inm);
      if(iswp) ped->a=MMG_swapbin(ped->a);
      fread(&ped->b,sw,1,inm);
      if(iswp) ped->b=MMG_swapbin(ped->b);
      fread(&ped->ref,sw,1,inm);
      if(iswp) ped->ref=MMG_swapbin(ped->ref);
    }
  }

  /* read triangles */
  if(mesh->nt) {
    rewind(inm);
    fseek(inm,posnt,SEEK_SET);
    norient = 0;
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if (!bin)
        fscanf(inm,"%d %d %d %d",&pt->v[0],&pt->v[1],&pt->v[2],&pt->ref);
      else {
        for (i=0 ; i<3 ; i++) {
          fread(&pt->v[i],sw,1,inm);
          if(iswp) pt->v[i]=MMG_swapbin(pt->v[i]);
        }
        fread(&pt->ref,sw,1,inm);
        if(iswp) pt->ref=MMG_swapbin(pt->ref);
      }
      for (i=0; i<3; i++) {
        ppt = &mesh->point[ pt->v[i] ];
        ppt->tag &= ~MG_NUL;
      }
      for(i=0 ; i<3 ; i++)
        pt->edg[i] = 0;
      air = MMG2_quickarea(mesh->point[pt->v[0]].c,mesh->point[pt->v[1]].c,
                           mesh->point[pt->v[2]].c);

      if(air < 0) {
        if ( mesh->info.ddebug || mesh->info.imprim > 6 )
          printf("Tr %d bad oriented\n",k);
        norient++;
        tmp = pt->v[2];
        pt->v[2] = pt->v[1];
        pt->v[1] = tmp;
      }
    }
    if ( norient ) {
      fprintf(stdout,"\n     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n");
      fprintf(stdout,"         BAD ORIENTATION : vol < 0 -- %8d element(s) reoriented\n",norient);
      fprintf(stdout,"     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n");
    }
  } else {
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[ k ];
      ppt->tag &= ~MG_NUL;
    }
  }

  /* read corners*/
  if (ncor) {
    rewind(inm);
    fseek(inm,posncor,SEEK_SET);
    for (k=1; k<=ncor; k++) {
      if (!bin)
        fscanf(inm,"%d",&ref);
      else {
        fread(&ref,sw,1,inm);
        if(iswp) ref=MMG_swapbin(ref);
      }
      ppt = &mesh->point[ref];
      ppt->tag |= M_CORNER;
    }
  }

  /* read required vertices*/
  if (nreq) {
    rewind(inm);
    fseek(inm,posreq,SEEK_SET);
    for (k=1; k<=nreq; k++) {
      if (!bin)
        fscanf(inm,"%d",&ref);
      else {
        fread(&ref,sw,1,inm);
        if(iswp) ref=MMG_swapbin(ref);
      }
      ppt = &mesh->point[ref];
      ppt->tag |= M_REQUIRED;
    }
  }

  /* read required edges*/
  if (nreqed) {
    rewind(inm);
    fseek(inm,posreqed,SEEK_SET);
    for (k=1; k<=nreqed; k++) {
      if (!bin)
        fscanf(inm,"%d",&ref);
      else {
        fread(&ref,sw,1,inm);
        if(iswp) ref=MMG_swapbin(ref);
      }
      ped = &mesh->edge[ref];
      ped->tag |= M_REQUIRED;
      ppt = &mesh->point[ped->a];
      ppt->tag |= M_REQUIRED;
      ppt = &mesh->point[ped->b];
      ppt->tag |= M_REQUIRED;
    }
  }

  fclose(inm);

  /*maill periodique : remettre toutes les coord entre 0 et 1*/
  if(mesh->info.renum==-10) {
    if ( mesh->info.imprim > 4 || mesh->info.ddebug )
      printf("  ## Periodic mesh: %d points %d triangles\n",mesh->np,mesh->nt);
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      while (ppt->c[0] > 1 + 5e-3) {
        ppt->c[0] -= 1;
      }
      while (ppt->c[0] < 0 - 5e-3) {
        ppt->c[0] += 1;
      }
      while (ppt->c[1] > 1 + 5e-3) {
        ppt->c[1] -= 1;
      }
      while (ppt->c[1] < 0 - 5e-3) {
        ppt->c[1] += 1;
      }
    }
  }


  return(1);
}

/* Load mesh file at gmsh format */
int MMG2D_loadMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename) {
  FILE*       inm;
  MMG5_pPoint ppt;
  double      z;
  long        posNodes,posElts,posNodeData;
  int         ier,k;
  int         bin,iswp,nelts;

  mesh->dim = 2;

  ier = MMG5_loadMshMesh_part1(mesh,sol,filename,&inm,
                               &posNodes,&posElts,&posNodeData,
                               &bin,&iswp,&nelts);
  if ( ier < 1 ) return (ier);

  if ( !MMG2D_zaldy(mesh) )  return(0);

  if ( mesh->ne || mesh->nprism ) {
    fprintf(stderr,"  ## Error: Input mesh must be a two-dimensional mesh.\n");
    return(-1);
  }
  if ( !mesh->nt )
      fprintf(stdout,"  **WARNING NO GIVEN TRIANGLE\n");

  if (mesh->npmax < mesh->np || mesh->ntmax < mesh->nt )
    return(-1);

  ier = MMG5_loadMshMesh_part2( mesh, sol,&inm,
                                posNodes,posElts,posNodeData,
                                bin,iswp,nelts);

  if ( ier < 1 ) return ( ier );

  z = 0.;
  for ( k=1; k<=mesh->np; ++k ) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;

    z += fabs(ppt->c[2]);
  }
  if ( z > _MMG5_EPSOK ) {
    fprintf(stderr,"  ## Error: Input mesh must be a two-dimensional mesh.\n");
    return(-1);
  }

  return(1);
}


/* load metric */
int MMG2D_loadSol(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename) {
  FILE       *inm;
  float       fsol;
  double      dsol;
  long        posnp;
  int         binch,bdim,iswp;
  int         k,i,isol,type,bin,dim,btyp,bpos;
  int         msh;
  char        *ptr,data[128],chaine[128];

  msh = mesh->info.nreg;

  bin = 0;

  strcpy(data,filename);

  ptr = strstr(data,".sol");
  if ( ptr ) {
    // filename contains the solution extension
    ptr = strstr(data,".solb");

    if ( ptr )  bin = 1;

    if( !(inm = fopen(data,"rb")) ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
      return(0);
    }
  }
  else
  {
    // filename don't contains the solution extension
    ptr = strstr(data,".mesh");
    if ( ptr ) *ptr = '\0';

    strcat(data,".solb");
    if (!(inm = fopen(data,"rb")) ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      if (!(inm = fopen(data,"rb")) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        return(0);
      }
    }
    else  bin = 1;
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);

  if(!bin) {
    strcpy(chaine,"DDD");
    while(fscanf(inm,"%s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if(!strncmp(chaine,"Dimension",strlen("Dimension"))) {
        fscanf(inm,"%d",&dim);
        if(dim!=2 && !msh) {
          fprintf(stdout,"  -- BAD SOL DIMENSION : %d\n",dim);
          return(-1);
        } else if(dim!=2 && msh) {
          if(dim==3)
            fprintf(stdout,"  -- READ 3D SOLUTION : %d\n",dim);
          else {
            fprintf(stdout,"  -- BAD SOL DIMENSION : %d\n",dim);
            return(-1);
          }
        }
        continue;
      } else if(!strncmp(chaine,"SolAtVertices",strlen("SolAtVertices"))) {
        fscanf(inm,"%d",&sol->np);
        fscanf(inm,"%d",&type);
        if(type!=1) {
          fprintf(stdout,"SEVERAL SOLUTION => IGNORED : %d\n",type);
          return(-1);
        }
        fscanf(inm,"%d",&btyp);
        posnp = ftell(inm);
        break;
      }
    }
  } else {
    fread(&binch,sw,1,inm);
    iswp=0;
    if(binch==16777216) iswp=1;
    else if(binch!=1) {
      fprintf(stdout,"BAD FILE ENCODING \n");
    }
    fread(&sol->ver,sw,1,inm);
    if(iswp) sol->ver = MMG_swapbin(sol->ver);
    while(fread(&binch,sw,1,inm)!=EOF && binch!=54 ) {
      if(iswp) binch=MMG_swapbin(binch);
      if(binch==54) break;
      if(binch==3) {  //Dimension
        fread(&bdim,sw,1,inm);  //NulPos=>20
        if(iswp) bdim=MMG_swapbin(bdim);
        fread(&bdim,sw,1,inm);
        if(iswp) bdim=MMG_swapbin(bdim);
        dim = bdim;
        if(bdim!=2) {
          fprintf(stdout,"BAD SOL DIMENSION : %d\n",bdim);
          exit(EXIT_FAILURE);
          return(-1);
        }
        continue;
      } else if(binch==62) {  //SolAtVertices
        fread(&binch,sw,1,inm); //NulPos
        if(iswp) binch=MMG_swapbin(binch);
        fread(&sol->np,sw,1,inm);
        if(iswp) sol->np=MMG_swapbin(sol->np);
        fread(&binch,sw,1,inm); //nb sol
        if(iswp) binch=MMG_swapbin(binch);
        if(binch!=1) {
          fprintf(stdout,"SEVERAL SOLUTION => IGNORED : %d\n",type);
          return(-1);
        }
        fread(&btyp,sw,1,inm); //typsol
        if(iswp) btyp=MMG_swapbin(btyp);
        posnp = ftell(inm);
        break;
      } else {
        fread(&bpos,sw,1,inm); //Pos
        if(iswp) bpos=MMG_swapbin(bpos);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }

  }

  if ( !sol->np ) {
    fprintf(stdout,"  ** MISSING DATA.\n");
    return(-1);
  }

  if ( sol->np != mesh->np ) {
    fprintf(stdout,"  ** WRONG DATA. IGNORED\n");
    return(-1);
  }

  if ( btyp!= 1 && btyp!=3 ) {
    fprintf(stdout,"  ** DATA IGNORED\n");
    sol->size = 1;
    sol->np = 0;
    return(-1);
  }
  sol->size = btyp;

  /* mem alloc */
  _MMG5_ADD_MEM(mesh,(sol->size*(mesh->npmax+1))*sizeof(double),
                "initial solution",return(0));
  _MMG5_SAFE_CALLOC(sol->m,(sol->size*(mesh->npmax+1)),double);

  /* read mesh solutions */
  rewind(inm);
  fseek(inm,posnp,SEEK_SET);
  for (k=1; k<=sol->np; k++) {
    isol = k * sol->size;
    if (sol->ver == 1) {
      for (i=0; i<sol->size; i++) {
        if(!bin){
          fscanf(inm,"%f",&fsol);
          sol->m[isol + i] = (double) fsol;
        } else {
          fread(&fsol,sw,1,inm);
          if(iswp) fsol=MMG_swapf(fsol);
          sol->m[isol + i] = (double) fsol;
        }
      }
      if(dim==3 && (sol->size>1)) {
        if(!bin){
          fscanf(inm,"%f",&fsol);
          fscanf(inm,"%f",&fsol);
          fscanf(inm,"%f",&fsol);
        } else {
          fread(&fsol,sw,1,inm);
          fread(&fsol,sw,1,inm);
          fread(&fsol,sw,1,inm);
        }
      }
    } else {
      for (i=0; i<sol->size; i++) {
        if(!bin){
          fscanf(inm,"%lf",&sol->m[isol + i]);

        } else {
          fread(&sol->m[isol + i],sd,1,inm);
          if(iswp) sol->m[isol + i]=MMG_swapd(sol->m[isol + i]);
        }
      }
      if(dim==3 && (sol->size>1)) {
        if(!bin){
          fscanf(inm,"%lf",&dsol);
          fscanf(inm,"%lf",&dsol);
          fscanf(inm,"%lf",&dsol);
        } else {
          fread(&dsol,sw,1,inm);
          fread(&dsol,sw,1,inm);
          fread(&dsol,sw,1,inm);
        }
      }
    }
  }

  sol->npi = sol->np;

  fclose(inm);
  return(1);
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 * \warning you must call the \a MMG2_tassage function before saving your
 * mesh.
 */
int MMG2D_saveMesh(MMG5_pMesh mesh,const char *filename) {
  FILE*             inm;
  MMG5_pPoint       ppt;
  MMG5_pEdge        ped;
  MMG5_pTria        pt;
  double            dblb;
  int               k,ne,ref/*,nn*/,ntang;
  int               bin, binch, bpos;
  char              *ptr,data[128],chaine[128];

  mesh->ver = 2;
  bin = 0;

  strcpy(data,filename);

  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if( !(inm = fopen(data,"wb")) ) {
      ptr  = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if( !(inm = fopen(data,"wb")) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        return(0);
      }
    } else {
      bin = 1;
    }
  }
  else {
    ptr = strstr(data,".meshb");
    if( ptr )  bin = 1;
    if( !(inm = fopen(data,"wb")) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
      return(0);
    }
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /*entete fichier*/
  binch=0; bpos=10;
  if (!bin) {
    strcpy(&chaine[0],"MeshVersionFormatted 2\n");
    fprintf(inm,"%s",chaine);
    if(mesh->info.nreg) {
      strcpy(&chaine[0],"\n\nDimension 3\n");
    } else {
      strcpy(&chaine[0],"\n\nDimension 2\n");
    }
    fprintf(inm,"%s ",chaine);
  }
  else
  {
    binch = 1; //MeshVersionFormatted
    fwrite(&binch,sw,1,inm);
    binch = 2; //version
    fwrite(&binch,sw,1,inm);
    binch = 3; //Dimension
    fwrite(&binch,sw,1,inm);
    bpos = 20; //Pos
    fwrite(&bpos,sw,1,inm);
    if(mesh->info.nreg) binch = 3; //Dimension
    else binch = 2;
    fwrite(&binch,sw,1,inm);
  }

  // //HACK SAVE ONLY SD ref 2
  //   for (k=1; k<=mesh->np; k++) {
  //     ppt = &mesh->point[k];
  //     ppt->tag = MG_NUL;
  //   }
  //   for (k=1; k<=mesh->nt; k++) {
  //     pt = &mesh->tria[k];
  //  if (!M_EOK(pt)) continue;
  //     if (pt->ref==2) {
  //    mesh->point[pt->v[0]].tag &= ~MG_NUL;
  //    mesh->point[pt->v[1]].tag &= ~MG_NUL;
  //    mesh->point[pt->v[2]].tag &= ~MG_NUL;
  //     } else {
  //    pt->v[0] = 0;
  //     }
  //   }
  //   for (k=1; k<=mesh->na; k++) {
  //     ped = &mesh->edge[k];
  //    if(!ped->a) continue;
  //     ppt = &mesh->point[ped->a];
  //    if(!M_VOK(ppt)) {
  //      ped->a = 0;
  //      continue;
  //    }
  //     ppt = &mesh->point[ped->b];
  //    if(!M_VOK(ppt)) {
  //      ped->a = 0;
  //      continue;
  //    }
  //    ped->tag |= M_REQUIRED;
  //    }
  //
  // //END HACK

  /* vertices */
  ne = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( M_VOK(ppt) )  ne++;
    ppt->tmp = ne;
  }

  if(!bin)
  {
    strcpy(&chaine[0],"\n\nVertices\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%d\n",ne);
  }
  else
  {
    binch = 4; //Vertices
    fwrite(&binch,sw,1,inm);
    if ( mesh->info.nreg )
      bpos += 12+(1+3*mesh->ver)*4*ne; //NullPos
    else
      bpos += 12+(1+2*mesh->ver)*4*ne; //NullPos

    fwrite(&bpos,sw,1,inm);
    fwrite(&ne,sw,1,inm);
  }
  fflush(inm);

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( M_VOK(ppt) ) {
      ref = ppt->ref;
      if(mesh->info.nreg) {
        if ( !bin )
          fprintf(inm,"%.15lg %.15lg 0. %d\n",ppt->c[0],ppt->c[1],ref);
        else {
          dblb = 0.;
          fwrite((unsigned char*)&ppt->c[0],sd,1,inm);
          fwrite((unsigned char*)&ppt->c[1],sd,1,inm);
          fwrite((unsigned char*)&dblb,sd,1,inm);
          fwrite((unsigned char*)&ref,sw,1,inm);
        }
      } else {
        if ( !bin ) {
          fprintf(inm,"%.15lg %.15lg %d\n",ppt->c[0],ppt->c[1],ref);
          fflush(inm);
        }
        else {
          fwrite(&ppt->c[0],sd,1,inm);
          fwrite(&ppt->c[1],sd,1,inm);
          fwrite(&ref,sw,1,inm);
        }
      }
    }
  }

  /* corners */
  ne = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( M_VOK(ppt) && (ppt->tag & M_CORNER) /*&& (ppt->tag & M_BDRY)*/)  ne++;
  }

  if ( ne ) {
    if ( !bin ) {
      strcpy(&chaine[0],"\n\nCorners\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",ne);
    }
    else
    {
      binch = 13; //
      fwrite(&binch,sw,1,inm);
      bpos += 12+4*ne; //NullPos
      fwrite(&bpos,sw,1,inm);
      fwrite(&ne,sw,1,inm);
    }

    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( M_VOK(ppt) && (ppt->tag & M_CORNER) /*&& (ppt->tag & M_BDRY)*/) {
        if(!bin) {
          fprintf(inm,"%d\n",ppt->tmp);
        }
        else {
          fwrite(&ppt->tmp,sw,1,inm);
        }
      }
    }
  }

  /* required vertex */
  ne = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( M_VOK(ppt) ) {
      if ( mesh->info.nosurf && (ppt->tag & M_NOSURF) ) continue;
      if (ppt->tag & M_REQUIRED) ne++;
    }
  }
  if ( ne ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nRequiredVertices\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",ne);
    }
    else {
      binch = 15; //
      fwrite(&binch,sw,1,inm);
      bpos += 12+4*ne; //NullPos
      fwrite(&bpos,sw,1,inm);
      fwrite(&ne,sw,1,inm);
    }
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( M_VOK(ppt) ) {
        if ( mesh->info.nosurf && ( ppt->tag & M_NOSURF )) continue;
        if ((ppt->tag & M_REQUIRED)
            /*&& ( (ppt->tag & M_BDRY) || (ppt->tag & M_SD) ) */ ) {
          if(!bin)
            fprintf(inm,"%d\n",ppt->tmp);
          else
            fwrite(&ppt->tmp,sw,1,inm);
        }
      }
    }
  }

  /* edges */
  // Warning : edges are not packed
  ne = 0;
  if ( mesh->na ) {
    for (k=1; k<=mesh->na; k++) {
      ped = &mesh->edge[k];
      if ( ped->a && ped->b ) ++ne;
    }
  }

  if ( ne ) {
    if(!bin) {
      strcpy(&chaine[0],"\n\nEdges\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",ne);
    }
    else {
      binch = 5; //Edges
      fwrite(&binch,sw,1,inm);
      bpos += 12 + 3*4*ne;//Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&ne,sw,1,inm);
    }
    for (k=1; k<=mesh->na; k++) {
      ped = &mesh->edge[k];
      if((!ped->a)||(!ped->b)) continue;

      ref = ped->ref;
      if(!bin)
        fprintf(inm,"%d %d %d\n",mesh->point[ped->a].tmp,mesh->point[ped->b].tmp,ref);
      else
      {
        fwrite(&mesh->point[ped->a].tmp,sw,1,inm);
        fwrite(&mesh->point[ped->b].tmp,sw,1,inm);
        fwrite(&ped->ref,sw,1,inm);
      }
    }

    /* required edges */
    ne = 0;
    for (k=1; k<=mesh->na; k++) {
      ped = &mesh->edge[k];
      if ( ped->a && (ped->tag & M_REQUIRED))  ne++;
    }
    if ( ne ) {
      if ( !bin ) {
        strcpy(&chaine[0],"\n\nRequiredEdges\n");
        fprintf(inm,"%s",chaine);
        fprintf(inm,"%d\n",ne);
      }
      else {
        binch = 16; //RequiredEdges
        fwrite(&binch,sw,1,inm);
        bpos += 12 + 4*ne;//Pos
        fwrite(&bpos,sw,1,inm);
        fwrite(&ne,sw,1,inm);
      }
      ne = 0;
      for (k=1; k<=mesh->na; k++) {
        ped = &mesh->edge[k];
        if ( !ped->a ) continue;
        ++ne;
        if ( (ped->tag & M_REQUIRED) ) {
          if ( !bin )
            fprintf(inm,"%d\n",ne);
          else {
            fwrite(&ne,sw,1,inm);
          }
        }
      }
    }
  }

  /* elements */
  ne = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( M_EOK(pt) )  ne++;
  }

  if ( ne ) {
    if ( !bin ) {
      strcpy(&chaine[0],"\n\nTriangles\n");
      fprintf(inm,"%s",chaine);
      fprintf(inm,"%d\n",ne);
    }
    else {
      binch = 6; //Triangles
      fwrite(&binch,sw,1,inm);
      bpos += 12+16*ne; //Pos
      fwrite(&bpos,sw,1,inm);
      fwrite(&ne,sw,1,inm);
    }
    for (k=1; k<=mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( M_EOK(pt) ) {
        ref = pt->ref;
        if ( !bin ) {
          fprintf(inm,"%d %d %d %d\n",mesh->point[pt->v[0]].tmp,
                  mesh->point[pt->v[1]].tmp,
                  mesh->point[pt->v[2]].tmp,ref);
        }
        else {
          fwrite(&mesh->point[pt->v[0]].tmp,sw,1,inm);
          fwrite(&mesh->point[pt->v[1]].tmp,sw,1,inm);
          fwrite(&mesh->point[pt->v[2]].tmp,sw,1,inm);
          fwrite(&ref,sw,1,inm);
        }
      }
    }
  }


  /*savetangent*/
  ntang=0;
  for(k=1 ; k<=mesh->np ; k++) {
    ppt = &mesh->point[k];
    if ( M_VOK(ppt) ) {
      if(!(ppt->tag & M_BDRY)) continue;
      if(ppt->tag & M_CORNER) continue;
      ntang++;
    }
  }


  /* Remark: here we save the tangents but there is a bug in medit (it crashes
   * if it try to read tangents without normals. It is easy to patch, in
   * zaldy1.c, seek the " if ( mesh->ntg )" field and replace
   * "assert(mesh->extra->n);" by "assert(mesh->extra->t);").
   * To not have to modify medit, here we save the tangents as if it were normals. */
  /* if ( ntang ) { */
  /*   if ( !bin ) { */
  /*     strcpy(&chaine[0],"\n\nNormals\n"); //be careful it is tangent!! */
  /*     fprintf(inm,"%s",chaine); */
  /*     fprintf(inm,"%d\n",ntang); */
  /*   } */
  /*   else */
  /*   { */
  /*     binch = 60; //normals */
  /*     fwrite(&binch,sw,1,inm); */
  /*     if ( mesh->info.nreg ) */
  /*       bpos += 12+(3*mesh->ver)*4*ntang; //Pos */
  /*     else */
  /*       bpos += 12+(2*mesh->ver)*4*ntang; //Pos */
  /*     fwrite(&bpos,sw,1,inm); */
  /*     fwrite(&ntang,sw,1,inm); */
  /*   } */

  /*   for(k=1 ; k<=mesh->np ; k++) { */
  /*     ppt = &mesh->point[k]; */
  /*     if(!M_VOK(ppt)) continue; */
  /*     if(!(ppt->tag & M_BDRY)) continue; */
  /*     if(ppt->tag & M_CORNER) continue; */
  /*     if(mesh->info.nreg) { */
  /*       if ( !bin ) */
  /*         fprintf(inm,"%lf %lf %lf\n",ppt->n[0],ppt->n[1],0.e0); */
  /*       else { */
  /*         dblb = 0; */
  /*         fwrite((unsigned char*)&ppt->n[0],sd,1,inm); */
  /*         fwrite((unsigned char*)&ppt->n[1],sd,1,inm); */
  /*         fwrite(&dblb,sd,1,inm); */
  /*       } */
  /*     } */
  /*     else */
  /*     { */
  /*       if ( !bin ) */
  /*         fprintf(inm,"%lf %lf \n",ppt->n[0],ppt->n[1]); */
  /*       else { */
  /*         fwrite((unsigned char*)&ppt->n[0],sd,1,inm); */
  /*         fwrite((unsigned char*)&ppt->n[1],sd,1,inm); */
  /*       } */
  /*     } */
  /*   } */

  /*   if ( !bin ) { */
  /*     strcpy(&chaine[0],"\n\nNormalAtVertices\n"); */
  /*     fprintf(inm,"%s",chaine); */
  /*     fprintf(inm,"%d\n",ntang); */
  /*   } */
  /*   else { */
  /*     binch = 20; //normalatvertices */
  /*     fwrite(&binch,sw,1,inm); */
  /*     bpos += 12 + 2*4*ntang;//Pos */
  /*     fwrite(&bpos,sw,1,inm); */
  /*     fwrite(&ntang,sw,1,inm); */
  /*   } */
  /*   nn=1; */
  /*   for(k=1 ; k<=mesh->np ; k++) { */
  /*     ppt = &mesh->point[k]; */
  /*     if ( !M_VOK(ppt) ) continue; */
  /*     if(!(ppt->tag & M_BDRY)) continue; */
  /*     if(ppt->tag & M_CORNER) continue; */

  /*     if(!bin) { */
  /*       fprintf(inm,"%d %d \n",ppt->tmp,nn++); */
  /*     } */
  /*     else { */
  /*       fwrite(&ppt->tmp,sw,1,inm); */
  /*       ++nn; */
  /*       fwrite(&nn,sw,1,inm); */
  /*     } */
  /*   } */
  /* } */

  if(!bin) {
    strcpy(&chaine[0],"\n\nEnd\n");
    fprintf(inm,"%s",chaine);
  }
  else {
    binch = 54; //End
    fwrite(&binch,sw,1,inm);
  }

  fclose(inm);

  return(1);
}

int MMG2D_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename) {
  return(MMG5_saveMshMesh(mesh,sol,filename));
}

int MMG2_loadVect(MMG5_pMesh mesh,const char *filename) {
  printf("comment for merge of the data structure\n");
  exit(EXIT_FAILURE);
  /* FILE       *inm; */
  /* Displ       pd; */
  /* float       fsol; */
  /* int         binch,bdim,iswp; */
  /* int         k,i,type,bin,dim,btyp,bpos,iadr; */
  /* long        posnp; */
  /* char        *ptr,data[128],chaine[128]; */

  /* pd = mesh->disp; */

  /* posnp = 0; */
  /* bin   = 0; */
  /* iswp  = 0; */

  /* strcpy(data,filename); */
  /* ptr = strstr(data,".mesh"); */
  /* if ( ptr )  *ptr = '\0'; */
  /* strcat(data,".solb"); */
  /* if( !(inm = fopen(data,"rb")) ) { */
  /*   ptr  = strstr(data,".solb"); */
  /*   *ptr = '\0'; */
  /*   strcat(data,".sol"); */
  /*   if( !(inm = fopen(data,"rb")) ) { */
  /*     fprintf(stderr,"  ** %s  NOT FOUND.\n",data); */
  /*     return(0); */
  /*   } */
  /* } else { */
  /*   bin = 1; */
  /* } */
  /* fprintf(stdout,"  %%%% %s OPENED\n",data); */

  /* if(!bin) { */
  /*   strcpy(chaine,"DDD"); */
  /*   while(fscanf(inm,"%s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) { */
  /*     if(!strncmp(chaine,"MeshVersionFormatted",strlen("MeshVersionFormatted"))) { */
  /*       fscanf(inm,"%d",&pd.ver); */
  /*       continue; */
  /*     } if(!strncmp(chaine,"Dimension",strlen("Dimension"))) { */
  /*  fscanf(inm,"%d",&dim); */
  /*  if(dim!=2) { */
  /*    fprintf(stdout,"BAD SOL DIMENSION : %d\n",dim); */
  /*    return(0); */
  /*  } */
  /*  continue; */
  /*     } else if(!strncmp(chaine,"SolAtVertices",strlen("SolAtVertices"))) { */
  /*  fscanf(inm,"%d",&pd.np); */
  /*  fscanf(inm,"%d",&type); */
  /*  if(type!=1) { */
  /*    fprintf(stdout,"SEVERAL SOLUTION => IGNORED : %d\n",type); */
  /*    return(0); */
  /*  } */
  /*  fscanf(inm,"%d",&btyp); */
  /*  posnp = ftell(inm); */
  /*  break; */
  /*     } */
  /*   } */
  /* } else { */
  /*   fread(&pd.ver,sw,1,inm); */
  /*   iswp=0; */
  /*   if(pd.ver==16777216) iswp=1; */
  /*   else if(pd.ver!=1) { */
  /*     fprintf(stdout,"BAD FILE ENCODING\n"); */
  /*   } */
  /*   fread(&pd.ver,sw,1,inm); */
  /*   if(iswp) pd.ver = MMG_swapbin(pd.ver); */
  /*   while(fread(&binch,sw,1,inm)!=EOF && binch!=54 ) { */
  /*     if(iswp) binch=MMG_swapbin(binch); */
  /*     if(binch==54) break; */
  /*     if(binch==3) {  //Dimension */
  /*  fread(&bdim,sw,1,inm);  //Pos=>20 */
  /*  if(iswp) bdim=MMG_swapbin(bdim); */
  /*  fread(&bdim,sw,1,inm); */
  /*  if(iswp) bdim=MMG_swapbin(bdim); */
  /*  dim = bdim */
  /*  if(bdim!=2) { */
  /*    fprintf(stdout,"BAD SOL DIMENSION : %d\n",dim); */
  /*    exit(EXIT_FAILURE); */
  /*    return(0); */
  /*  } */
  /*  continue; */
  /*     } else if(binch==62) {  //SolAtVertices */
  /*  fread(&binch,sw,1,inm); //Pos */
  /*  if(iswp) binch=MMG_swapbin(binch); */
  /*  fread(&pd.np,sw,1,inm); */
  /*  if(iswp) pd.np=MMG_swapbin(pd.np); */
  /*  fread(&binch,sw,1,inm); //nb sol */
  /*  if(iswp) binch=MMG_swapbin(binch); */
  /*  if(binch!=1) { */
  /*    fprintf(stdout,"SEVERAL SOLUTION => IGNORED : %d\n",type); */
  /*    return(0); */
  /*  } */
  /*  fread(&btyp,sw,1,inm); //typsol */
  /*  if(iswp) btyp=MMG_swapbin(btyp); */
  /*  posnp = ftell(inm); */
  /*  break; */
  /*     } else { */
  /*  fread(&bpos,sw,1,inm); //Pos */
  /*  if(iswp) bpos=MMG_swapbin(bpos); */
  /*  rewind(inm); */
  /*  fseek(inm,bpos,SEEK_SET); */
  /*     } */
  /*   } */

  /* } */
  /* if ( !pd.np ) { */
  /*   fprintf(stdout,"  ** MISSING DATA\n"); */
  /*   return(0); */
  /* } */
  /* else if ( pd.np != mesh->np ) { */
  /*   fprintf(stdout,"  ** WRONG DATA\n"); */
  /*   return(0); */
  /* } */

  /* if ( btyp != 2 ) { */
  /*   fprintf(stdout,"  ** DISPLACEMENT DATA IGNORED\n"); */
  /*   return(0); */
  /* } */

  /* if ( abs(mesh->info.imprim) > 3 ) */
  /*   fprintf(stdout,"  -- READING DISPLACEMENT FILE %s\n",data); */

  /* /\* read mesh solutions *\/   */
  /* rewind(inm); */
  /* fseek(inm,posnp,SEEK_SET); */
  /* for (k=1; k<=pd.np; k++) { */
  /*   iadr = (k - 1) * btyp + 1; */
  /*   if (pd.ver < 2) { */
  /*     for (i=0; i<2; i++) { */
  /*  if(!bin){ */
  /*    fscanf(inm,"%f",&fsol); */
  /*    pd.mv[iadr + i] = (double) fsol; */
  /*  } else { */
  /*    fread(&fsol,sw,1,inm); */
  /*    if(iswp) fsol=MMG_swapf(fsol); */
  /*    pd.mv[iadr + i] = (double) fsol; */
  /*  } */
  /*     } */
  /*   } else { */
  /*     for (i=0; i<2; i++) { */
  /*  if(!bin){ */
  /*    fscanf(inm,"%lf",&pd.mv[iadr + i]); */
  /*  } else { */
  /*    fread(&pd.mv[iadr + i],sd,1,inm); */
  /*    if(iswp) pd.mv[iadr + i]=MMG_swapd(pd.mv[iadr + i]); */
  /*  } */
  /*     } */
  /*   } */
  /* } */

  /* if ( abs(mesh->info.imprim) > 3 ) */
  /*   fprintf(stdout,"     NUMBER OF GIVEN DATA       %8d\n",pd.np); */

  /* fclose(inm); */
  return(1);
}
int MMG2D_saveSol(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename) {
  FILE*        inm;
  MMG5_pPoint       ppt;
  float        fsol;
  double       dsol;
  int          i,k,nbl,isol,bin,bpos,typ;
  char        *ptr,data[128],chaine[128];
  int          binch,msh;

  msh = mesh->info.nreg;
  if ( !sol->np )  return(1);
  bin = 0;

  strcpy(data,filename);

  ptr = strstr(data,".sol");
  if ( ptr ) {
    // filename contains the solution extension
    ptr = strstr(data,".solb");

    if ( ptr )  bin = 1;

    if( !(inm = fopen(data,"wb")) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
      return(0);
    }
  }
  else
  {
    // filename don't contains the solution extension
    ptr = strstr(data,".mesh");
    if ( ptr ) *ptr = '\0';

    strcat(data,".sol");
    if (!(inm = fopen(data,"wb")) ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      if (!(inm = fopen(data,"wb")) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        return(0);
      }
      else bin = 1;
    }
  }

  fprintf(stdout,"  %%%% %s OPENED\n",data);

  /*entete fichier*/
  if(!bin) {
    strcpy(&chaine[0],"MeshVersionFormatted 2\n");
    fprintf(inm,"%s",chaine);
    if(msh) {
      strcpy(&chaine[0],"\n\nDimension 3\n");
    } else {
      strcpy(&chaine[0],"\n\nDimension 2\n");
    }
    fprintf(inm,"%s ",chaine);
  } else {
    binch = 1; //MeshVersionFormatted
    fwrite(&binch,sw,1,inm);
    binch = 2; //version
    fwrite(&binch,sw,1,inm);
    binch = 3; //Dimension
    fwrite(&binch,sw,1,inm);
    bpos = 20; //Pos
    fwrite(&bpos,sw,1,inm);
    if (msh )
      binch = 3;
    else
      binch = 2;
    fwrite(&binch,sw,1,inm);

  }

  switch(sol->size) {
  case 1:
    typ = 1;
    break;
  case 3:
    typ = 3;
    break;
  default:
    fprintf(stdout,"  ** DATA IGNORED\n");
    return(0);
  }

  /* write data */
  nbl = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !M_VOK(ppt) )  continue;
    nbl++;
  }
  if(!bin) {
    strcpy(&chaine[0],"\n\nSolAtVertices\n");
    fprintf(inm,"%s",chaine);
    fprintf(inm,"%d\n",nbl);
    fprintf(inm,"%d %d\n",1,typ);
  } else {
    binch = 62; //Vertices
    fwrite(&binch,sw,1,inm);
    bpos += 20+(sol->size*sol->ver)*4*nbl; //Pos
    fwrite(&bpos,sw,1,inm);
    fwrite(&nbl,sw,1,inm);
    binch = 1; //nb sol
    fwrite(&binch,sw,1,inm);
    binch = typ; //typ sol
    fwrite(&binch,sw,1,inm);
  }

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( !M_VOK(ppt) )  continue;
    isol = k * sol->size;
    if (sol->ver < 2) {
      if(!bin) {
        if(msh && sol->size > 1) {
          for (i=0; i<sol->size; i++) {
            fsol = (float) sol->m[isol + i];
            fprintf(inm,"%f ",fsol);
          }
          fprintf(inm,"%f %f %f",0.,0.,1.);
        } else {
          for (i=0; i<sol->size; i++) {
            fsol = (float) sol->m[isol + i];
            fprintf(inm,"%f ",fsol);
          }
        }
        fprintf(inm,"\n");
      } else {
        if(msh && sol->size > 1) {
          for (i=0; i<sol->size; i++) {
            fsol = (float) sol->m[isol + i];
            fwrite(&fsol,sw,1,inm);
          }
          fsol = 0;
          fwrite(&fsol,sw,1,inm);
          fsol = 0;
          fwrite(&fsol,sw,1,inm);
          fsol = 1;
          fwrite(&fsol,sw,1,inm);
        } else {
          for (i=0; i<sol->size; i++) {
            fsol = (float) sol->m[isol + i];
            fwrite(&fsol,sw,1,inm);
          }
        }
      }
    } else {
      if(!bin) {
        if(msh && sol->size > 1) {
          for (i=0; i<sol->size; i++)
            fprintf(inm,"%.15lg ",sol->m[isol + i]);
          fprintf(inm,"%.15lg %.15lg %.15lg",0.,0.,1.);
        } else {
          for (i=0; i<sol->size; i++)
            fprintf(inm,"%.15lg ",sol->m[isol + i]);
        }
        fprintf(inm,"\n");
      } else {
        if(msh && sol->size > 1) {
          for (i=0; i<sol->size; i++)
            fwrite(&sol->m[isol + i],sd,1,inm);
          dsol = 0.;
          fwrite(&dsol,sd,1,inm);
          dsol = 0.;
          fwrite(&dsol,sd,1,inm);
          dsol = 0.;
          fwrite(&dsol,sd,1,inm);
        } else {
          for (i=0; i<sol->size; i++)
            fwrite(&sol->m[isol + i],sd,1,inm);
        }
      }

    }
  }

  /*fin fichier*/
  if(!bin) {
    strcpy(&chaine[0],"\n\nEnd\n");
    fprintf(inm,"%s",chaine);
  } else {
    binch = 54; //End
    fwrite(&binch,sw,1,inm);
  }
  fclose(inm);

  return(1);
}

int MMG2D_saveVect(MMG5_pMesh mesh,MMG5_pSol sol,const char *filename,double lambda) {
  printf("comment for merge of the data structure\n");
  exit(EXIT_FAILURE);
  /* FILE        *inm,*f,*inm2; */
  /* Displ        pd; */
  /* MMG5_pPoint       ppt; */
  /* double       dsol,dd,lreal; */
  /* int          i,k,nbl,bin,bpos,typ; */
  /* char        *ptr,data[128],data2[128],chaine[128]; */
  /* unsigned char binch; */

  /* if(lambda > 1) { */
  /*   lreal = 1.; */
  /* } else { */
  /*   lreal = lambda; */
  /* } */
  /* /\*ecriture du lambda dans un fichier*\/ */
  /* f = fopen("lambda.dat","w"); */
  /* fprintf(f,"%10.12f\n",lambda); */
  /* fclose(f); */

  /* pd = mesh->disp; */
  /* pd.ver = 2; */
  /* bin = 0; */

  /* strcpy(data,filename); */
  /* strcpy(data2,filename); */
  /* ptr = strstr(data,".meshb"); */
  /* if ( ptr )  *ptr = '\0'; */
  /* else { */
  /*   ptr = strstr(data,".mesh"); */
  /*   if ( ptr ) { */
  /*     *ptr = '\0'; */
  /*     bin  = 0; */
  /*   } */
  /* } */
  /* strcat(data,"-lambda.sol"); */
  /* ptr = strstr(data2,".meshb"); */
  /* if ( ptr )  *ptr = '\0'; */
  /* else { */
  /*   ptr = strstr(data2,".mesh"); */
  /*   if ( ptr ) { */
  /*     *ptr = '\0'; */
  /*     bin  = 0; */
  /*   } */
  /* } */
  /* strcat(data2,"-1-lambda.sol"); */

  /* if( bin && !(inm = fopen(data,"wb")) ) { */
  /*   fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data); */
  /*   return(0); */
  /* } else { */
  /*   if( !(inm = fopen(data,"w")) ) { */
  /*     fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data); */
  /*     return(0); */
  /*   } */
  /* } */
  /* fprintf(stdout,"  %%%% %s OPENED\n",data); */

  /*  if( bin && !(inm2 = fopen(data2,"wb")) ) { */
  /*   fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data2); */
  /*   return(0); */
  /* } else { */
  /*   if( !(inm2 = fopen(data2,"w")) ) { */
  /*     fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data2); */
  /*     return(0); */
  /*   } */
  /* } */
  /* fprintf(stdout,"  %%%% %s OPENED\n",data2); */

  /* /\*entete fichier*\/ */
  /* if(!bin) { */
  /*   strcpy(&chaine[0],"MeshVersionFormatted 2\n"); */
  /*   fprintf(inm,"%s",chaine); */
  /*   strcpy(&chaine[0],"\n\nDimension 2\n"); */
  /*   fprintf(inm,"%s ",chaine); */
  /* } else { */
  /*   binch = 1; //MeshVersionFormatted */
  /*   fwrite(&binch,sw,1,inm); */
  /*   binch = pd.ver; //version */
  /*   fwrite(&binch,sw,1,inm); */
  /*   binch = 2; //Dimension */
  /*   fwrite(&binch,sw,1,inm); */
  /*   bpos = 20; //Pos */
  /*   fwrite(&bpos,sw,1,inm); */
  /*   binch = 3; */
  /*   fwrite(&binch,sw,1,inm); */

  /* } */
  /* /\*entete fichier 2*\/ */
  /* if(!bin) { */
  /*   strcpy(&chaine[0],"MeshVersionFormatted 2\n"); */
  /*   fprintf(inm2,"%s",chaine); */
  /*   strcpy(&chaine[0],"\n\nDimension 2\n"); */
  /*   fprintf(inm2,"%s ",chaine); */
  /* } else { */
  /*   binch = 1; //MeshVersionFormatted */
  /*   fwrite(&binch,sw,1,inm2); */
  /*   binch = pd.ver; //version */
  /*   fwrite(&binch,sw,1,inm2); */
  /*   binch = 2; //Dimension */
  /*   fwrite(&binch,sw,1,inm2); */
  /*   bpos = 20; //Pos */
  /*   fwrite(&bpos,sw,1,inm2); */
  /*   binch = 3; */
  /*   fwrite(&binch,sw,1,inm2); */

  /* } */
  /* typ = 2; */

  /* /\* write data *\/ */
  /* nbl = 0; */
  /* for (k=1; k<=mesh->np; k++) { */
  /*   ppt = &mesh->point[k]; */
  /*   if ( !M_VOK(ppt) )  continue; */
  /*   nbl++; */
  /* } */

  /* if(!bin) { */
  /*   strcpy(&chaine[0],"\n\nSolAtVertices\n"); */
  /*   fprintf(inm,"%s",chaine); */
  /*   fprintf(inm,"%d\n",nbl); */
  /*   fprintf(inm,"%d %d\n",1,typ); */
  /*   fprintf(inm2,"%s",chaine); */
  /*   fprintf(inm2,"%d\n",nbl); */
  /*   fprintf(inm2,"%d %d\n",1,typ); */
  /* } else { */
  /*   binch = 62; //SolAtVertices */
  /*   fwrite(&binch,sw,1,inm); */
  /*   bpos += 20+(2*pd.ver)*4*nbl; //Pos */
  /*   fwrite(&bpos,sw,1,inm); */
  /*   fwrite(&nbl,sw,1,inm); */
  /*   binch = 1; //nb sol */
  /*   fwrite(&binch,sw,1,inm); */
  /*   binch = typ; //typ sol */
  /*   fwrite(&binch,sw,1,inm); */
  /*   binch = 62; //SolAtVertices */
  /*   fwrite(&binch,sw,1,inm2); */
  /*   bpos += 20+(2*pd.ver)*4*nbl; //Pos */
  /*   fwrite(&bpos,sw,1,inm2); */
  /*   fwrite(&nbl,sw,1,inm2); */
  /*   binch = 1; //nb sol */
  /*   fwrite(&binch,sw,1,inm2); */
  /*   binch = typ; //typ sol */
  /*   fwrite(&binch,sw,1,inm2); */
  /* } */


  /* for (k=1; k<=mesh->np; k++) { */
  /*   ppt = &mesh->point[k]; */
  /*   if ( !M_VOK(ppt) )  continue; */
  /*   for (i=0; i<2; i++) { */
  /*     dsol = lreal*pd.mv[2*(k-1) + 1 + i]; */
  /*     if(!bin) { */
  /*  fprintf(inm,"%.15lg ",dsol); */
  /*     } else { */
  /*  fwrite(&dsol,sd,1,inm); */
  /*     } */
  /*   } */
  /*   for (i=0; i<2; i++) { */
  /*     dsol = (1.-lreal)*pd.mv[2*(k-1) + 1 + i]; */
  /*     if(!bin) { */
  /*  fprintf(inm2,"%.15lg ",dsol); */
  /*     } else { */
  /*  fwrite(&dsol,sd,1,inm2); */
  /*     } */
  /*   } */
  /*   if(!bin) { */
  /*     fprintf(inm,"\n"); */
  /*     fprintf(inm2,"\n"); */
  /*   } */
  /* } */
  /* /\*fin fichier*\/ */
  /* if(!bin) { */
  /*   strcpy(&chaine[0],"\n\nEnd\n"); */
  /*   fprintf(inm,"%s",chaine); */
  /*   fprintf(inm2,"%s",chaine); */
  /* } else { */
  /*   binch = 54; //End */
  /*   fwrite(&binch,sw,1,inm); */
  /*   fwrite(&binch,sw,1,inm2); */
  /* } */
  /* fclose(inm); */
  /* fclose(inm2); */
  return(1);
}
