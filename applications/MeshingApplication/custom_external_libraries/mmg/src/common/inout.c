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
 * \file common/inout.c
 * \brief Input / Output Functions.
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "mmgcommon.h"

#define sw 4
#define sd 8

static int MMG5_swapbin(int sbin)
{
  int inv;
  char *p_in = (char *) &sbin;
  char *p = (char *)&inv;


  p[0] = p_in[3];
  p[1] = p_in[2];
  p[2] = p_in[1];
  p[3] = p_in[0];

  return inv;
  /*unsigned char c1, c2, c3, c4;

    c1 = sbin & 255;
    c2 = (sbin >> 8) & 255;
    c3 = (sbin >> 16) & 255;
    c4 = (sbin >> 24) & 255;

    return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;   */

}
static float MMG5_swapf(float sbin)
{
  float out;
  char *p_in = (char *) &sbin;
  char *p_out = (char *) &out;
  p_out[0] = p_in[3];
  p_out[1] = p_in[2];
  p_out[2] = p_in[1];
  p_out[3] = p_in[0];

  return out;
}
static double MMG5_swapd(double sbin)
{
  float out;
  char *p_in = (char *) &sbin;
  char *p_out = (char *) &out;
  int i;

  for(i=0;i<8;i++)
  {
    p_out[i] = p_in[7-i];
  }

  return out;
}

static
int MMG5_countBinaryElts(FILE **inm, const int nelts,const int iswp,
                          int *np, int *na, int* nt,int *nq, int *ne, int *npr)
{
  int    typ,num,tagNum,i,k,l,idx;
  static char mmgWarn = 0;

  k = 0;

  while ( k<nelts ) {
    fread(&typ,sw,1,(*inm));
    if(iswp) typ = MMG5_swapbin(typ);

    switch (typ) {
    case 1:
      /* Edge */
      fread(&num,sw,1,(*inm));
      fread(&tagNum,sw,1,(*inm));
      if(iswp) {
        num = MMG5_swapbin(num);
        tagNum = MMG5_swapbin(tagNum);
      }
      for ( idx=0; idx<num; ++idx ) {
        fread(&i,sw,1,(*inm));
        for ( l=0; l<tagNum; ++l ) fread(&i,sw,1,(*inm));
        fread(&i,sw,1,(*inm)); // edge->a
        fread(&i,sw,1,(*inm)); // edge->b
      }
      (*na) += num;
      k  += num;
      break;
    case 2:
      /* Tria */
      fread(&num,sw,1,(*inm));
      fread(&tagNum,sw,1,(*inm));
      if(iswp) {
        num = MMG5_swapbin(num);
        tagNum = MMG5_swapbin(tagNum);
      }
      for ( idx=0; idx<num; ++idx ) {
        fread(&i,sw,1,(*inm));
        for ( l=0; l<tagNum; ++l ) fread(&i,sw,1,(*inm));
        fread(&i,sw,1,(*inm)); // tria->v[0]
        fread(&i,sw,1,(*inm)); // tria->v[1]
        fread(&i,sw,1,(*inm)); // tria->v[2]
      }
      (*nt) += num;
      k  += num;
      break;
    case 3:
      /* Quad */
      fread(&num,sw,1,(*inm));
      fread(&tagNum,sw,1,(*inm));
      if(iswp) {
        num = MMG5_swapbin(num);
        tagNum = MMG5_swapbin(tagNum);
      }
      for ( idx=0; idx<num; ++idx ) {
        fread(&i,sw,1,(*inm));
        for ( l=0; l<tagNum; ++l ) fread(&i,sw,1,(*inm));
        fread(&i,sw,1,(*inm)); // quadra->v[0]
        fread(&i,sw,1,(*inm)); // quadra->v[1]
        fread(&i,sw,1,(*inm)); // quadra->v[2]
        fread(&i,sw,1,(*inm)); // quadra->v[3]
      }
      (*nq) += num;
      k  += num;
      break;
    case 4:
      /* Tetra */
      fread(&num,sw,1,(*inm));
      fread(&tagNum,sw,1,(*inm));
      if(iswp) {
        num = MMG5_swapbin(num);
        tagNum = MMG5_swapbin(tagNum);
      }
      for ( idx=0; idx<num; ++idx ) {
        fread(&i,sw,1,(*inm));
        for ( l=0; l<tagNum; ++l ) fread(&i,sw,1,(*inm));
        fread(&i,sw,1,(*inm)); // tetra->v[0]
        fread(&i,sw,1,(*inm)); // tetra->v[1]
        fread(&i,sw,1,(*inm)); // tetra->v[2]
        fread(&i,sw,1,(*inm)); // tetra->v[3]
      }
      (*ne) += num;
      k  += num;
      break;
    case 6:
      /* Prism */
      fread(&num,sw,1,(*inm));
      fread(&tagNum,sw,1,(*inm));
      if(iswp) {
        num = MMG5_swapbin(num);
        tagNum = MMG5_swapbin(tagNum);
      }
      for ( idx=0; idx<num; ++idx ){
        fread(&i,sw,1,(*inm));
        for ( l=0; l<tagNum; ++l ) fread(&i,sw,1,(*inm));
        fread(&i,sw,1,(*inm)); // prism->v[0]
        fread(&i,sw,1,(*inm)); // prism->v[1]
        fread(&i,sw,1,(*inm)); // prism->v[2]
        fread(&i,sw,1,(*inm)); // prism->v[3]
        fread(&i,sw,1,(*inm)); // prism->v[4]
        fread(&i,sw,1,(*inm)); // prism->v[5]
      }
      (*npr) += num;
      k  += num;
      break;
    case 15:
      /* Node */
      fread(&num,sw,1,(*inm));
      fread(&tagNum,sw,1,(*inm));
      if(iswp) {
        num = MMG5_swapbin(num);
        tagNum = MMG5_swapbin(tagNum);
      }
      for ( idx=0; idx<num; ++idx ){
        fread(&i,sw,1,(*inm));
        for ( l=0; l<tagNum; ++l ) fread(&i,sw,1,(*inm));
        fread(&i,sw,1,(*inm)); //node idx
      }
      (*np) += num;
      k  += num;
      break;
    default:
      if ( !mmgWarn ) {
        fprintf(stderr,"\n  ## Warning: %s: unexpected type of element (%d) for at"
                " least 1 element (%d).\n",__func__,typ,k);
        mmgWarn = 1;
      }
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param filename pointer toward the name of file
 * \param inm pointer toward the file pointer
 * \param posNodes pointer toward the position of nodes data in file
 * \param posElts pointer toward the position of elts data in file
 * \param posNodeData pointer toward the list of the positions of data in file
 * \param bin 1 if binary format
 * \param nelts number of elements in file
 * \param nsol number of data in file
 * \return 1 if success, 0 if file is not found, -1 if fail.
 *
 * Begin to read mesh at MSH file format. Read the mesh size informations.
 *
 */
int MMG5_loadMshMesh_part1(MMG5_pMesh mesh,const char *filename,
                           FILE **inm,
                           long *posNodes, long *posElts,
                           long **posNodeData, int *bin, int *iswp,
                           int *nelts,int *nsols) {
  double      dbuf[9];
  float       fbuf[9];
  int         ver,oneBin,k,i;
  int         nt,na,nq,ne,npr,np;
  int         typ,tagNum,posNodeDataSize,initPosNodeDataSize;
  char        *ptr,*data,chaine[128],verNum[5];

  ver = oneBin = 0;
  *posNodes = 0;
  *posElts = 0;
  *nelts = 0;
  *nsols = 0;
  *bin = 0;
  *iswp = 0;
  mesh->np = mesh->nt = mesh->ne = 0;
  nt = na = nq = ne = npr = np = 0;

  MMG5_SAFE_CALLOC(data,strlen(filename)+7,char,return 0);

  /* Allocation of the posNodeData array: we assume that we have less than 20
   * solutions in the file (for a greater number of sol, posNoteData is
   * reallocated) */
  initPosNodeDataSize = posNodeDataSize = 20;
  MMG5_SAFE_CALLOC(*posNodeData,posNodeDataSize,long,return 0);

  strcpy(data,filename);
  ptr = strstr(data,".msh");
  if ( !ptr ) {
    /* data contains the filename without extension */
    strcat(data,".mshb");
    if (!(*inm = fopen(data,"rb")) ) {
      ptr  = strstr(data,".msh");
      *ptr = '\0';
      strcat(data,".msh");
      if( !((*inm) = fopen(data,"rb")) ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
        MMG5_SAFE_FREE(data);
        MMG5_SAFE_FREE(*posNodeData);
        return 0;
      }
    }
  }
  else {
    if( !((*inm) = fopen(data,"rb")) ) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
      MMG5_SAFE_FREE(data);
      MMG5_SAFE_FREE(*posNodeData);
      return 0;
    }
  }

  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);


  /* Detection of the different fields of the file */
  strcpy(chaine,"D");
  while(fscanf((*inm),"%127s ",&chaine[0])!=EOF ) {
    if(!strncmp(chaine,"$MeshFormat",strlen("$MeshFormat"))) {
      fscanf((*inm),"%4s %d %d ",verNum,bin,&ver);
      mesh->ver = ver/4;
      if ( strncmp(verNum,"2.2",3) ) {
        fprintf(stderr,"\n  ## Warning: %s: format version (%s) may be not supported."
                " Please, use the format version 2.2.\n",__func__,verNum);
      }
      if ( *bin ) {
        fread(&oneBin,sw,1,(*inm));
        if ( oneBin!=1 ) {
          assert(oneBin==16777216);
          *iswp=1;
          oneBin = MMG5_swapbin(oneBin);
        }
      }
      continue;
    } else if(!strncmp(chaine,"$EndMeshFormat",strlen("EndMeshFormat"))) {
      continue;
    } else if(!strncmp(chaine,"$Nodes",strlen("$Nodes"))) {
      fscanf((*inm),"%d ",&mesh->npi);
      *posNodes = ftell((*inm));
      if ( *bin ) {
        /* Skip the binary nodes data */
        if ( mesh->ver==1 ) {
          for ( k=1; k<=mesh->npi; ++k ) {
            fread(&i,sw,1,(*inm));
            fread( &fbuf[0],sw,3,(*inm) );
          }
        }
        else {
          for ( k=1; k<=mesh->npi; ++k ) {
            fread(&i,sw,1,(*inm));
            fread( &dbuf[0],sd,3,(*inm) );
          }
        }
      }
      continue;
    } else if(!strncmp(chaine,"$EndNodes",strlen("$EndNodes"))) {
      continue;
    } else if(!strncmp(chaine,"$NodeData",strlen("$NodeData"))) {

      (*posNodeData)[*nsols] = ftell((*inm));
      if ( ++(*nsols) == posNodeDataSize ) {
        MMG5_SAFE_RECALLOC(*posNodeData,*nsols,
                            posNodeDataSize+initPosNodeDataSize,
                            long,"posNodeData",return 0);
        posNodeDataSize += initPosNodeDataSize;
      }

      if ( *bin ) {
        /* Skip the binary nodes data */
        /* String tags */
        fscanf((*inm),"%d ",&tagNum);
        for ( k=0; k<tagNum; ++k ) {
          fscanf((*inm),"%*[^\n]%*c");
        }

        /* Real tags */
        fscanf((*inm),"%d ",&tagNum);
        for ( k=0; k<tagNum; ++k ) {
          fscanf((*inm),"%*[^\n]%*c");
        }

        /* Integer tags */
        fscanf((*inm),"%d ",&tagNum);
        fscanf((*inm),"%d ",&i); //time step
        fscanf((*inm),"%d ",&typ); //type of solution
        fscanf((*inm),"%d ",&np);
        for ( k=3; k<tagNum; ++k ) {
          fscanf((*inm),"%d",&i);
        }
        if ( mesh->ver==1 ) {
          for ( k=1; k<=np; ++k ) {
            fread(&i,sw,1,(*inm));
            fread( &fbuf[0],sw,typ,(*inm) );
          }
        }
        else {
          for ( k=1; k<=np; ++k ) {
            fread(&i,sw,1,(*inm));
            fread( &dbuf[0],sd,typ,(*inm) );
          }
        }
      }
      continue;
    } else if(!strncmp(chaine,"$EndNodeData",strlen("$EndNodeData"))) {
      continue;
    } else if(!strncmp(chaine,"$Elements",strlen("$Elements"))) {
      fscanf((*inm),"%d ",nelts);
      *posElts = ftell((*inm));

      /* Count the elements */
      if ( !*bin ) {
        for ( k=0; k<*nelts; ++k) {
          fscanf((*inm),"%d %d ",&i,&typ);
          switch (typ) {
          case 1:
            /* Edge */
            ++na;
            break;
          case 2:
            /* Tria */
            ++nt;
            break;
          case 3:
            /* Quad */
            ++nq;
            break;
          case 4:
            /* Tetra */
            ++ne;
            break;
          case 6:
            /* Prism */
            ++npr;
            break;
          case 15:
            /* Node */
            ++np;
            break;
          }
          fscanf((*inm),"%*[^\n]%*c");
        }
      }
      else {
        if ( !MMG5_countBinaryElts(inm,*nelts,*iswp,&np,&na,&nt,&nq,&ne,&npr) ) {
          fclose(*inm);
          MMG5_SAFE_FREE(*posNodeData);
          return -1;
        }
      }
      continue;
    }
    else if(!strncmp(chaine,"$EndElements",strlen("$EndElements"))) {
      continue;
    }
  }

  if ( !mesh->npi ) {
    fprintf(stderr,"  ** MISSING DATA.\n");
    fprintf(stderr,"     Check that your mesh contains points and elements.\n");
    fprintf(stderr,"     Exit program.\n");
    fclose(*inm);
    MMG5_SAFE_FREE(*posNodeData);
    return -1;
  }

  /* memory allocation */
  mesh->np = mesh->npi;
  mesh->nt = mesh->nti = nt;
  mesh->ne = mesh->nei = ne;
  mesh->na = mesh->nai = na;
  mesh->nprism = npr;
  mesh->nquad = nq;

  if ( !mesh->np ) {
    fprintf(stderr,"  ** MISSING DATA.\n");
    fprintf(stderr,"     Check that your mesh contains points.\n");
    fprintf(stderr,"     Exit program.\n");
    fclose(*inm);
    MMG5_SAFE_FREE(*posNodeData);
    return -1;
  }

  return 1;
}


/**
 * \param mesh pointer toward the mesh
 * \param sol pointer toward the solutions array
 * \param inm pointer toward the file pointer
 * \param posNodes position of nodes data in file
 * \param posElts position of elts data in file
 * \param posNodeData position of solution data in file
 * \param bin 1 if binary format
 * \param nelts number of elements in file
 * \param nsols number of silutions in file
 * \return 1 if success, 0 if fail.
 *
 * End to read mesh and solution array at MSH file format after the
 * mesh/solution array alloc.
 *
 */
int MMG5_loadMshMesh_part2(MMG5_pMesh mesh,MMG5_pSol *sol,FILE **inm,
                           const long posNodes,const long posElts,
                           const long *posNodeData,const int bin,const int iswp,
                           const int nelts,const int nsols) {
  MMG5_pTetra pt;
  MMG5_pPrism pp;
  MMG5_pTria  ptt;
  MMG5_pQuad  pq1;
  MMG5_pEdge  pa;
  MMG5_pPoint ppt;
  MMG5_pSol   psl;
  double      aux, dbuf[9];
  float       fbuf[9],fc;
  int         k,i,l,nref,iadr;
  int         *ina_t,*ina_a,nt,na,nq,ne,npr;
  int         nbl_t,nbl_a,typ,tagNum,ref,idx,num;
  int         v[4],isol;
  char        chaine[128],*ptr;
  static char mmgWarn=0, mmgWarn1=0;

  ina_t = ina_a = NULL;

  /** Second step: read the nodes and elements */
  rewind((*inm));
  fseek((*inm),posNodes,SEEK_SET);

  if ( mesh->ver < 2 ) {
    for ( k=0; k< mesh->np; ++k)
    {
      if ( !bin ) {
        fscanf((*inm),"%d ",&idx);
        ppt = &mesh->point[idx];
        for (i=0 ; i<mesh->dim ; i++) {
          fscanf((*inm),"%f ",&fc);
          ppt->c[i] = (double) fc;
        }
      }
      else {
        fread(&idx,sw,1,(*inm));
        if ( iswp ) idx = MMG5_swapbin(idx);
        ppt = &mesh->point[idx];
        for (i=0 ; i<mesh->dim ; i++) {
          fread(&fc,sw,1,(*inm));
          if(iswp) fc=MMG5_swapf(fc);
          ppt->c[i] = (double) fc;
        }
      }
      ppt->tag  = MG_NUL;
      ppt->tmp  = 0;
      ppt->ref = 0;
    }
  }
  else {
    for ( k=0; k< mesh->np; ++k)
    {
      if ( !bin ) {
        fscanf((*inm),"%d ",&i);
        ppt = &mesh->point[i];
        fscanf((*inm),"%lf %lf %lf ",&ppt->c[0],&ppt->c[1],&ppt->c[2]);
      }
      else {
        fread(&i,sw,1,(*inm));
        if ( iswp ) i = MMG5_swapbin(i);
        ppt = &mesh->point[i];
        for (i=0 ; i<3 ; i++) {
          fread(&ppt->c[i],sd,1,(*inm));
          if(iswp) ppt->c[i]=MMG5_swapd(ppt->c[i]);
        }
      }
      ppt->tag  = MG_NUL;
      ppt->tmp  = 0;
      ppt->ref = 0;
    }
  }

  rewind((*inm));
  fseek((*inm),posElts,SEEK_SET);

  nbl_a = nbl_t = nt = na = nq = ne = npr = 0;
  nref = 0;

  /* Skip triangles and edges with MG_ISO refs */
  if( mesh->info.iso ) {
    if ( mesh->nt ) {
      MMG5_SAFE_CALLOC(ina_t,mesh->nt+1,int,return 0);
    }
    if ( mesh->na ) {
      MMG5_SAFE_CALLOC(ina_a,mesh->na+1,int,return 0);
    }
  }

  if ( !bin ) {
    for ( k=0; k<nelts; ++k)
    {
      fscanf((*inm),"%d %d %d ",&i,&typ, &tagNum);
      if ( tagNum < 2 ) {
        fprintf(stderr,"\n  ## Error: %s: elt %d (type %d): Expected at least 2 tags (%d given).\n",
                __func__,k,typ,tagNum);
        if ( ina_t ) MMG5_SAFE_FREE(ina_t);
        if ( ina_a ) MMG5_SAFE_FREE(ina_a);
        fclose(*inm);
        return -1;
      }
      fscanf((*inm),"%d %d ",&ref,&i);
      for ( l=2; l<tagNum; ++l ) {
        fscanf((*inm),"%d ",&i);
      }

      switch (typ) {
      case 1:
        /* Edge */
        /* Skip edges with MG_ISO refs */
        if ( !mesh->info.iso ) {
          pa = &mesh->edge[++na];
          fscanf((*inm),"%d %d ",&pa->a,&pa->b);
          pa->ref = abs(ref);
          pa->tag |= MG_REF;
        }
        else {
          if ( abs(ref)!= MG_ISO ) {
            pa = &mesh->edge[++na];
            fscanf((*inm),"%d %d ",&pa->a,&pa->b);
            pa->ref = abs(ref);
            pa->tag |= MG_REF;
            ina_a[na+nbl_a]=na;
          }
          else {
            /* Skip this edge but advance the file pointer */
            pa = &mesh->edge[0];
            fscanf((*inm),"%d %d ",&pa->a,&pa->b);
            ++nbl_a;
          }
        }
        assert( na+nbl_a<=mesh->na );
        break;
      case 2:
        /* Tria */
        /* Skip triangles with MG_ISO refs */
        if ( !mesh->info.iso ) {
          ptt = &mesh->tria[++nt];
          fscanf((*inm),"%d %d %d ",&ptt->v[0],&ptt->v[1],&ptt->v[2]);
          ptt->ref = abs(ref);
        }
        else {
          if ( abs(ref)!= MG_ISO ) {
            ptt = &mesh->tria[++nt];
            fscanf((*inm),"%d %d %d",&ptt->v[0],&ptt->v[1],&ptt->v[2]);
            ptt->ref = abs(ref);
            ina_t[nt+nbl_t]=nt;
          }
          else {
            /* Skip this triangle but advance the file pointer */
            ptt = &mesh->tria[0];
            fscanf((*inm),"%d %d %d",&ptt->v[0],&ptt->v[1],&ptt->v[2]);
            ++nbl_t;
          }
        }
        assert( nt+nbl_t<=mesh->nt );
        break;
      case 3:
        /* Quad */
        pq1 = &mesh->quadra[++nq];
        fscanf((*inm),"%d %d %d %d ",&pq1->v[0],&pq1->v[1],&pq1->v[2],&pq1->v[3]);
        pq1->ref = ref;
        assert( nq<=mesh->nquad );
        break;
      case 4:
        /* Tetra for mmg3d */
        if ( mesh->ne ) {
          pt = &mesh->tetra[++ne];
          fscanf((*inm),"%d %d %d %d ",&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3]);
          pt->ref = abs(ref);
        } else { /*skip tetra*/
          fscanf((*inm),"%d %d %d %d ",&v[0],&v[1],&v[2],&v[3]);
        }

        if(ref < 0) {
          nref++;
        }

        assert( ne<=mesh->ne );
        break;
      case 6:
        /* Prism for mmg3d */
        if ( mesh->nprism )
        {
          pp = &mesh->prism[++npr];
          fscanf((*inm),"%d %d %d %d %d %d ",&pp->v[0],&pp->v[1],&pp->v[2],
                 &pp->v[3],&pp->v[4],&pp->v[5]);
          pp->ref = abs(ref);
        }
        if(ref < 0) {
          nref++;
        }
        assert( npr<=mesh->nprism );
        break;
      case 15:
        /* Node */
        fscanf((*inm),"%d ",&l);
        ppt = &mesh->point[l];
        ppt->ref = ref;
        assert( l<=mesh->np );
        break;
      default:
        if ( !mmgWarn ) {
          fprintf(stderr,"\n  ## Warning: %s: unexpected type for at least 1 element:"
                  " element %d, type %d\n",__func__,k,typ );
          mmgWarn = 1;
        }
      }
    }
  }
  else {
    k = 0;

    while ( k<nelts ) {
      fread(&typ,sw,1,(*inm));
      if(iswp) typ = MMG5_swapbin(typ);

      switch (typ) {
      case 1:
        /* Edge */
        fread(&num,sw,1,(*inm));
        fread(&tagNum,sw,1,(*inm));
        if(iswp) {
          num = MMG5_swapbin(num);
          tagNum = MMG5_swapbin(tagNum);
        }
        if ( tagNum < 2 ) {
          fprintf(stderr,"\n  ## Error: %s: Expected at least 2 tags per element (%d given).\n",
                  __func__,tagNum);
          if ( ina_t ) MMG5_SAFE_FREE(ina_t);
          if ( ina_a ) MMG5_SAFE_FREE(ina_a);
          fclose(*inm);
          return -1;
        }

        for ( idx=0; idx<num; ++idx ) {
          fread(&i,sw,1,(*inm));

          fread(&ref,sw,1,(*inm));
          fread(&i,sw,1,(*inm));

          for ( l=2; l<tagNum; ++l ) {
            fread(&i,sw,1,(*inm));
          }

          if(iswp) ref = MMG5_swapbin(ref);

          /* Skip edges with MG_ISO refs */
          if ( !mesh->info.iso ) {
            pa = &mesh->edge[++na];
            fread(&pa->a,sw,1,(*inm));
            fread(&pa->b,sw,1,(*inm));
            if ( iswp ) {
              pa->a = MMG5_swapbin(pa->a);
              pa->b = MMG5_swapbin(pa->b);
            }
            pa->ref = abs(ref);
            pa->tag |= MG_REF;
          }
          else {
            if( abs(ref)!= MG_ISO ) {
              pa = &mesh->edge[++na];
              fread(&pa->a,sw,1,(*inm));
              fread(&pa->b,sw,1,(*inm));
              if ( iswp ) {
                pa->a = MMG5_swapbin(pa->a);
                pa->b = MMG5_swapbin(pa->b);
              }
              pa->ref = abs(ref);
              pa->tag |= MG_REF;
              ina_a[na+nbl_a]=na;
            }
            else {
              /* Skip this edge but advance the file pointer */
              pa = &mesh->edge[0];
              fread(&pa->a,sw,1,(*inm));
              fread(&pa->b,sw,1,(*inm));
              ++nbl_a;
            }
          }
          assert( na+nbl_a<=mesh->na );
        }
        k += num;

        break;
      case 2:
        /* Tria */
        fread(&num,sw,1,(*inm));
        fread(&tagNum,sw,1,(*inm));
        if(iswp) {
          num = MMG5_swapbin(num);
          tagNum = MMG5_swapbin(tagNum);
        }
        if ( tagNum < 2 ) {
          fprintf(stderr,"\n  ## Error: %s: Expected at least 2 tags per element (%d given).\n",
                  __func__,tagNum);
          if ( ina_t ) MMG5_SAFE_FREE(ina_t);
          if ( ina_a ) MMG5_SAFE_FREE(ina_a);
          fclose(*inm);
          return -1;
        }

        for ( idx=0; idx<num; ++idx ) {
          fread(&i,sw,1,(*inm));

          fread(&ref,sw,1,(*inm));
          fread(&i,sw,1,(*inm));

          for ( l=2; l<tagNum; ++l ) {
            fread(&i,sw,1,(*inm));
          }

          if(iswp) ref = MMG5_swapbin(ref);

          /* Skip triangles with MG_ISO refs */
          if ( !mesh->info.iso ) {
            ptt = &mesh->tria[++nt];
            for ( i=0; i<3 ; ++i ) {
              fread(&l,sw,1,(*inm));
              if ( iswp ) l = MMG5_swapbin(l);
              ptt->v[i] = l;
            }
            ptt->ref = abs(ref);
          }
          else {
            if( abs(ref)!= MG_ISO ) {
              ptt = &mesh->tria[++nt];
              for ( i=0; i<3 ; ++i ) {
                fread(&l,sw,1,(*inm));
                if ( iswp ) l = MMG5_swapbin(l);
                ptt->v[i] = l;
              }
              ptt->ref = abs(ref);
              ina_t[nt+nbl_t]=nt;
            }
            else {
              /* Skip this triangle but advance the file pointer */
              for ( i=0; i<3 ; ++i ) {
                fread(&l,sw,1,(*inm));
              }
              ++nbl_t;
            }
          }
          assert( nt+nbl_t<=mesh->nt );
        }
        k += num;
        break;
      case 3:
        /* Quad */
        fread(&num,sw,1,(*inm));
        fread(&tagNum,sw,1,(*inm));
        if(iswp) {
          num = MMG5_swapbin(num);
          tagNum = MMG5_swapbin(tagNum);
        }
       if ( tagNum < 2 ) {
          fprintf(stderr,"\n  ## Error: %s: Expected at least 2 tags per element (%d given).\n",
                  __func__,tagNum);
          if ( ina_t ) MMG5_SAFE_FREE(ina_t);
          if ( ina_a ) MMG5_SAFE_FREE(ina_a);
          fclose(*inm);
          return -1;
        }

        for ( idx=0; idx<num; ++idx ) {
          fread(&i,sw,1,(*inm));

          fread(&ref,sw,1,(*inm));
          fread(&i,sw,1,(*inm));

          for ( l=2; l<tagNum; ++l ) {
            fread(&i,sw,1,(*inm));
          }

          if(iswp) ref = MMG5_swapbin(ref);

          pq1 = &mesh->quadra[++nq];
          for ( i=0; i<4 ; ++i ) {
            fread(&l,sw,1,(*inm));
            if ( iswp ) l = MMG5_swapbin(l);
            pq1->v[i] = l;
          }
          pq1->ref = ref;
          assert( nq<=mesh->nquad );
        }
        k += num;
        break;
      case 4:
        /* Tetra */
        fread(&num,sw,1,(*inm));
        fread(&tagNum,sw,1,(*inm));
        if(iswp) {
          num = MMG5_swapbin(num);
          tagNum = MMG5_swapbin(tagNum);
        }

        if ( tagNum < 2 ) {
          fprintf(stderr,"\n  ## Error: %s: Expected at least 2 tags per element (%d given).\n",
                  __func__,tagNum);
          if ( ina_t ) MMG5_SAFE_FREE(ina_t);
          if ( ina_a ) MMG5_SAFE_FREE(ina_a);
          fclose(*inm);
          return -1;
        }

        for ( idx=0; idx<num; ++idx ) {
          fread(&i,sw,1,(*inm));

          fread(&ref,sw,1,(*inm));
          fread(&i,sw,1,(*inm));

          for ( l=2; l<tagNum; ++l ) {
            fread(&i,sw,1,(*inm));
          }

          if(iswp) ref = MMG5_swapbin(ref);

          if ( mesh->ne ) {
            pt = &mesh->tetra[++ne];
            for ( i=0; i<4 ; ++i ) {
              fread(&l,sw,1,(*inm));
              if ( iswp ) l = MMG5_swapbin(l);
              pt->v[i] = l;
            }
            pt->ref = abs(ref);
            assert( ne<=mesh->ne );
          }
          else
            for ( i=0; i<4 ; ++i )
              fread(&l,sw,1,(*inm));

          if(ref < 0) {
            nref++;
          }
        }
        k += num;
        break;
      case 6:
        /* Prism */
        fread(&num,sw,1,(*inm));
        fread(&tagNum,sw,1,(*inm));
        if(iswp) {
          num = MMG5_swapbin(num);
          tagNum = MMG5_swapbin(tagNum);
        }
        if ( tagNum < 2 ) {
          fprintf(stderr,"\n  ## Error: %s: Expected at least 2 tags per element (%d given).\n",
                  __func__,tagNum);
          if ( ina_t ) MMG5_SAFE_FREE(ina_t);
          if ( ina_a ) MMG5_SAFE_FREE(ina_a);
          fclose(*inm);
          return -1;
        }

        for ( idx=0; idx<num; ++idx ) {
          fread(&i,sw,1,(*inm));

          fread(&ref,sw,1,(*inm));
          fread(&i,sw,1,(*inm));

          for ( l=2; l<tagNum; ++l ) {
            fread(&i,sw,1,(*inm));
          }

          if(iswp) ref = MMG5_swapbin(ref);

          if ( mesh->nprism ) {
            pp = &mesh->prism[++npr];
            for ( i=0; i<6 ; ++i ) {
              fread(&l,sw,1,(*inm));
              if ( iswp ) l = MMG5_swapbin(l);
              pp->v[i] = l;
            }
            pp->ref = abs(ref);
            assert( npr<=mesh->nprism );
          }
          else {
            for ( i=0; i<6 ; ++i )
              fread(&l,sw,1,(*inm));
          }

          if(ref < 0) {
            nref++;
          }
        }
        k += num;
        break;
      case 15:
        /* Node */
        fread(&num,sw,1,(*inm));
        fread(&tagNum,sw,1,(*inm));
        if(iswp) {
          num = MMG5_swapbin(num);
          tagNum = MMG5_swapbin(tagNum);
        }
        if ( tagNum < 2 ) {
          fprintf(stderr,"\n  ## Error: %s: Expected at least 2 tags per element (%d given).\n",
                  __func__,tagNum);
          if ( ina_t ) MMG5_SAFE_FREE(ina_t);
          if ( ina_a ) MMG5_SAFE_FREE(ina_a);
          fclose(*inm);
          return -1;
        }

        for ( idx=0; idx<num; ++idx ) {
          fread(&i,sw,1,(*inm));

          fread(&ref,sw,1,(*inm));
          fread(&i,sw,1,(*inm));

          for ( l=2; l<tagNum; ++l ) {
            fread(&i,sw,1,(*inm));
          }

          if(iswp) ref = MMG5_swapbin(ref);

          fread(&l,sw,1,(*inm));
          if(iswp) l = MMG5_swapbin(l);
          ppt = &mesh->point[l];
          ppt->ref = ref;
          assert( l<=mesh->np );
        }
        k += num;
        break;
      default:
        fprintf(stderr,"\n  ## Error: %s: unexpected type of element (%d)\n",
                __func__,typ);
        if ( ina_t ) MMG5_SAFE_FREE(ina_t);
        if ( ina_a ) MMG5_SAFE_FREE(ina_a);
        fclose(*inm);
        return -1;
      }
    }
  }


  if ( mesh->dim==3 && mesh->info.iso ) {
    if ( mesh->nt ) {
      if( !nt )
        MMG5_DEL_MEM(mesh,mesh->tria);

      else if ( nt < mesh->nt ) {
        MMG5_ADD_MEM(mesh,(nt-mesh->nt)*sizeof(MMG5_Tria),"triangles",
                      fprintf(stderr,"  Exit program.\n");
                      fclose(*inm);
                      if ( ina_t ) MMG5_SAFE_FREE(ina_t);
                      if ( ina_a ) MMG5_SAFE_FREE(ina_a);
                      return 0);
        MMG5_SAFE_RECALLOC(mesh->tria,mesh->nt+1,(nt+1),MMG5_Tria,"triangles",
                            return 0);
      }
      MMG5_SAFE_FREE(ina_t);
      mesh->nt = nt;
    }
    if ( mesh->na ) {
      if( !na )
        MMG5_DEL_MEM(mesh,mesh->edge);
      else if ( na < mesh->na ) {
        MMG5_ADD_MEM(mesh,(na-mesh->na)*sizeof(MMG5_Edge),"edges",
                      fprintf(stderr,"  Exit program.\n");
                      fclose(*inm);
                      if ( ina_t ) MMG5_SAFE_FREE(ina_t);
                      if ( ina_a ) MMG5_SAFE_FREE(ina_a);
                      return 0);
        MMG5_SAFE_RECALLOC(mesh->edge,mesh->na+1,(na+1),MMG5_Edge,"edges",
                            return 0);
      }
      MMG5_SAFE_FREE(ina_a);
      mesh->na = na;
    }
  }


  if(nref) {
    fprintf(stdout,"\n     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n");
    fprintf(stdout,"         WARNING : %d elements (tetra or prisms) with ref < 0.",nref);
    fprintf(stdout," We take their absolute values.\n");
    fprintf(stdout,"     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n");
  }

  /* Cross the elements, mark their vertices as used and reorient the elements
   * with negative orientation */
  mesh->xt = 0;

  if ( mesh->dim == 2 ) {
    for (k=1; k<=mesh->nt; k++) {
      ptt = &mesh->tria[k];

      for (i=0; i<3; i++) {
        ppt = &mesh->point[ ptt->v[i] ];
        ppt->tag &= ~MG_NUL;
      }

      /* Set the elements references to 0 in iso mode */
      if ( mesh->info.iso )  ptt->ref = 0;

      for(i=0 ; i<3 ; i++)
        ptt->edg[i] = 0;

      if ( MMG2D_quickarea(mesh->point[ptt->v[0]].c,mesh->point[ptt->v[1]].c,
                          mesh->point[ptt->v[2]].c) < 0.0 ) {
        /* mesh->xt temporary used to count reoriented tetra*/
        mesh->xt++;
        aux = ptt->v[2];
        ptt->v[2] = ptt->v[1];
        ptt->v[1] = aux;
      }
    }
  }
  else {
    if ( mesh->ne ) {
      /* mmg3d */
      for ( k=1; k<=mesh->ne; k++ ) {
        pt = &mesh->tetra[k];
        if ( !MG_EOK(pt) ) continue;

        for (i=0; i<4; i++) {
          ppt = &mesh->point[pt->v[i]];
          ppt->tag &= ~MG_NUL;
        }

        /* Set the elements references to 0 in iso mode */
        if ( mesh->info.iso )  pt->ref = 0;

        /* Possibly switch 2 vertices number so that each tet is positively oriented */
        if ( MMG5_orvol(mesh->point,pt->v) < 0.0 ) {
          /* mesh->xt temporary used to count reoriented tetra*/
          mesh->xt++;
          aux = pt->v[2];
          pt->v[2] = pt->v[3];
          pt->v[3] = aux;
        }
      }
    }
    else {
      /* mmgs */
      for ( k=1; k<=mesh->nt; k++ ) {
        ptt = &mesh->tria[k];
        if ( !MG_EOK(ptt) ) continue;
        for (i=0; i<3; i++) {
          ppt = &mesh->point[ptt->v[i]];
          ppt->tag &= ~MG_NUL;
        }
      }
    }
  }

  if(mesh->xt) {
    fprintf(stdout,"\n     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n");
    fprintf(stdout,"         BAD ORIENTATION : vol < 0 -- %8d element(s) reoriented\n",mesh->xt);
    fprintf(stdout,"     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n");
  }
  mesh->xt = 0;

  /* Cross the prisms and mark their vertices as used */
  for ( k=1; k<=mesh->nprism; k++ ) {
    pp = &mesh->prism[k];
    for (i=0; i<6; i++) {
      ppt = &mesh->point[pp->v[i]];
      ppt->tag &= ~MG_NUL;
    }
  }

  /* stats */
  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"     NUMBER OF VERTICES       %8d\n",mesh->np);
    if ( mesh->ne )
      fprintf(stdout,"     NUMBER OF TETRAHEDRA     %8d\n",mesh->ne);

    if ( mesh->nprism )
      fprintf(stdout,"     NUMBER OF PRISMS         %8d\n",mesh->nprism);
    if ( mesh->nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES      %8d\n",mesh->nt);
    if ( mesh->nquad )
      fprintf(stdout,"     NUMBER OF QUADRILATERALS %8d\n",mesh->nquad);
    if ( mesh->na ) {
      fprintf(stdout,"     NUMBER OF EDGES          %8d\n",mesh->na);
    }
  }

  /** Read the solution at nodes */
  /* Init (*sol)[0] for the case where nsols=0 */
  psl = *sol;
  psl->ver = mesh->ver;
  psl->dim = mesh->dim;
  psl->type = 1;

  for ( isol=0; isol < nsols; ++isol ) {
    assert ( posNodeData[isol] );

    rewind((*inm));
    fseek((*inm),posNodeData[isol],SEEK_SET);

    psl = *sol + isol;

    psl->ver = mesh->ver;
    psl->dim = mesh->dim;
    psl->type = 1;

    /* String tags: The first one stores the solution name */
    fscanf((*inm),"%d ",&tagNum);
    if ( 1 != fscanf(*inm,"\"%127s\"\n",&chaine[0]) ) {
      fscanf(*inm,"%127s\n",&chaine[0]);
    }
    ptr = NULL;
    ptr = strstr(chaine,":metric");

    if ( ptr ) {
      *ptr = '\0';
      mesh->info.inputMet = 1;
    }


    if ( !MMG5_Set_inputSolName(mesh,psl,chaine) ) {
      if ( !mmgWarn1 ) {
        mmgWarn1 = 1;
        fprintf(stderr,"\n  ## Warning: %s: unable to set solution name for"
                " at least 1 solution.\n",__func__);
      }
    }

    for ( k=1; k<tagNum; ++k ) {
      fscanf((*inm),"%*[^\n]%*c");
    }

    /* Real tags ignored */
    if ( fscanf((*inm),"%d",&tagNum) ) {
      for ( k=0; k<tagNum; ++k ) {
        if ( 1 != fscanf((*inm),"%f",&fc) ) {
          fprintf(stderr,"\n  ## Error: unable to skip real tags\n");
          fclose(*inm);
          return -1;
        }
      }
    }

    /* Integer tags : allow to recover the number of sols and their types */
    fscanf((*inm),"%d ",&tagNum);
    if ( tagNum < 3 ) {
      fprintf(stderr,"   Error: %s: node data: Expected at least 3 tags (%d given).\n",
              __func__,tagNum);
      fclose(*inm);
      return -1;
    }

    fscanf((*inm),"%d ",&i); //time step;
    fscanf((*inm),"%d ",&typ); //type of solution: 1=scalar, 3=vector, 9=tensor ;
    fscanf((*inm),"%d ",&psl->np);

    for ( k=3; k<tagNum; ++k ) {
      fscanf((*inm),"%d",&i);
    }

    if ( mesh->np != psl->np ) {
      fprintf(stderr,"  ** MISMATCHES DATA: THE NUMBER OF VERTICES IN "
              "THE MESH (%d) DIFFERS FROM THE NUMBER OF VERTICES IN "
              "THE SOLUTION (%d) \n",mesh->np,psl->np);
      fclose(*inm);
      return -1;
    }

    if ( typ == 1 ) {
        psl->size = 1;
        psl->type = 1;
    }
    else if ( typ == 3 ) {
      psl->size = psl->dim;
      psl->type = 2;
    }
    else if ( typ == 9 ) {
      psl->size = (psl->dim*(psl->dim+1))/2;
      psl->type = 3;
    }
    else {
      fprintf(stderr,"  ** DATA TYPE IGNORED %d \n",typ);
      fclose(*inm);
      return -1;
    }

    /* mem alloc */
    if ( psl->m )  MMG5_DEL_MEM(mesh,psl->m);
    psl->npmax = mesh->npmax;

    MMG5_ADD_MEM(mesh,(psl->size*(psl->npmax+1))*sizeof(double),"initial solution",
                  fprintf(stderr,"  Exit program.\n");
                  fclose(*inm);
                  return 0);
    MMG5_SAFE_CALLOC(psl->m,psl->size*(psl->npmax+1),double,return 0);

    /* isotropic solution */
    if ( psl->size == 1 ) {
      if ( psl->ver == 1 ) {
        for (k=1; k<=psl->np; k++) {
          if(!bin){
            fscanf((*inm),"%d ",&idx);
            fscanf((*inm),"%f ",&fbuf[0]);
          } else {
            fread(&idx,sw,1, (*inm));
            if(iswp) idx = MMG5_swapbin(idx);
            fread(&fbuf[0],sw,1,(*inm));
            if(iswp) fbuf[0]=MMG5_swapf(fbuf[0]);
          }
          psl->m[idx] = fbuf[0];
        }
      }
      else {
        for (k=1; k<=psl->np; k++) {
          if(!bin){
            fscanf((*inm),"%d ",&idx);
            fscanf((*inm),"%lf ",&dbuf[0]);
          } else {
            fread(&idx,sw,1, (*inm));
            if(iswp) idx = MMG5_swapbin(idx);
            fread(&dbuf[0],sd,1,(*inm));
            if(iswp) dbuf[0]=MMG5_swapd(dbuf[0]);
          }
          psl->m[idx] = dbuf[0];
        }
      }
    }
    /* vector displacement only */
    else if ( psl->size == psl->dim ) {
      if ( psl->ver == 1 ) {
        for (k=1; k<=psl->np; k++) {
          if(!bin){
            fscanf((*inm),"%d ",&idx);
            for (i=0; i<psl->dim; i++) {
              fscanf((*inm),"%f ",&fbuf[0]);
              psl->m[psl->dim*idx+i] = fbuf[0];
            }
          } else {
            fread(&idx,sw,1, (*inm));
            if(iswp) idx = MMG5_swapbin(idx);
            for (i=0; i<psl->dim; i++) {
              fread(&fbuf[0],sw,1,(*inm));
              if(iswp) fbuf[0]=MMG5_swapf(fbuf[0]);
              psl->m[psl->dim*idx+i] = fbuf[0];
            }
          }
        }
      }
      else {
        for (k=1; k<=psl->np; k++) {
          if(!bin){
            fscanf((*inm),"%d ",&idx);

            for (i=0; i<psl->dim; i++) {
              fscanf((*inm),"%lf ",&dbuf[0]);
              psl->m[psl->dim*idx+i] = dbuf[0];
            }

          } else {
            fread(&idx,sw,1, (*inm));
            if(iswp) idx = MMG5_swapbin(idx);

            for (i=0; i<psl->dim; i++) {
              fread(&dbuf[0],sd,1,(*inm));
              if(iswp) dbuf[0]=MMG5_swapd(dbuf[0]);
              psl->m[psl->dim*idx+i] = dbuf[0];
            }

          }
        }
      }
    }
    /* anisotropic sol */
    else {
      if ( psl->ver == 1 ) {
        for (k=1; k<=psl->np; k++) {
          if(!bin){
            fscanf((*inm),"%d ",&idx);
            for(i=0 ; i<9 ; i++)
              fscanf((*inm),"%f ",&fbuf[i]);
          } else {
            fread(&idx,sw,1, (*inm));
            if(iswp) idx = MMG5_swapbin(idx);
            for(i=0 ; i<9 ; i++) {
              fread(&fbuf[i],sw,1,(*inm));
              if(iswp) fbuf[i]=MMG5_swapf(fbuf[i]);
            }
          }

          assert(fbuf[1]==fbuf[3] && fbuf[2]==fbuf[6] && fbuf[5]==fbuf[7]);

          if ( psl->dim ==2 ) {
            iadr = 3*idx;
            psl->m[iadr] = fbuf[0];
            psl->m[iadr+1] = fbuf[1];
            psl->m[iadr+2] = fbuf[4];
          }
          else {
            iadr = 6*idx;
            psl->m[iadr+0] = fbuf[0];
            psl->m[iadr+1] = fbuf[1];
            psl->m[iadr+2] = fbuf[2];
            psl->m[iadr+3] = fbuf[4];
            psl->m[iadr+4] = fbuf[5];
            psl->m[iadr+5] = fbuf[8];
          }
        }
      }
      else {
        for (k=1; k<=psl->np; k++) {
          if(!bin){
            fscanf((*inm),"%d ",&idx);
            for(i=0 ; i<9 ; i++)
              fscanf((*inm),"%lf ",&dbuf[i]);
          } else {
            fread(&idx,sw,1, (*inm));
            if(iswp) idx = MMG5_swapbin(idx);
            for(i=0 ; i<9 ; i++) {
              fread(&dbuf[i],sd,1,(*inm));
              if(iswp) dbuf[i]=MMG5_swapf(dbuf[i]);
            }
          }

          assert(dbuf[1]==dbuf[3] && dbuf[2]==dbuf[6] && dbuf[5]==dbuf[7]);

          if ( psl->dim ==2 ) {
            iadr = 3*idx;
            psl->m[iadr  ] = dbuf[0];
            psl->m[iadr+1] = dbuf[1];
            psl->m[iadr+2] = dbuf[4];
          }
          else {
            iadr = 6*idx;
            psl->m[iadr  ] = dbuf[0];
            psl->m[iadr+1] = dbuf[1];
            psl->m[iadr+2] = dbuf[2];
            psl->m[iadr+3] = dbuf[4];
            psl->m[iadr+4] = dbuf[5];
            psl->m[iadr+5] = dbuf[8];
          }
        }
      }
    }

    psl->npi = psl->np;

  }

  fclose((*inm));

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward the sol structure.
 * \param index of point in which we want to build the metric
 * \param dbuf builded metric
 *
 * Build the metric at point \a ip depending with its type (ridge/not ridge).
 *
 */
static inline
void  MMG5_build3DMetric(MMG5_pMesh mesh,MMG5_pSol sol,int ip,
                         double dbuf[6]) {
  MMG5_pPoint ppt;
  double      mtmp[3],r[3][3];
  int         i;

  ppt = &mesh->point[ip];
  if ( !(MG_SIN(ppt->tag) || (ppt->tag & MG_NOM) || (ppt->tag & MG_NOSURF))
       && (ppt->tag & MG_GEO) ) {
    if ( mesh->xp ) {
      // Arbitrary, we take the metric associated to the surface ruled by n_1
      mtmp[0] = sol->m[sol->size*(ip)];
      mtmp[1] = sol->m[sol->size*(ip)+1];
      mtmp[2] = sol->m[sol->size*(ip)+3];

      // Rotation matrix.
      r[0][0] = ppt->n[0];
      r[1][0] = ppt->n[1];
      r[2][0] = ppt->n[2];
      r[0][1] = mesh->xpoint[ppt->xp].n1[1]*ppt->n[2]
        - mesh->xpoint[ppt->xp].n1[2]*ppt->n[1];
      r[1][1] = mesh->xpoint[ppt->xp].n1[2]*ppt->n[0]
        - mesh->xpoint[ppt->xp].n1[0]*ppt->n[2];
      r[2][1] = mesh->xpoint[ppt->xp].n1[0]*ppt->n[1]
        - mesh->xpoint[ppt->xp].n1[1]*ppt->n[0];
      r[0][2] = mesh->xpoint[ppt->xp].n1[0];
      r[1][2] = mesh->xpoint[ppt->xp].n1[1];
      r[2][2] = mesh->xpoint[ppt->xp].n1[2];

      // Metric in the canonic space
      dbuf[0] = mtmp[0]*r[0][0]*r[0][0] + mtmp[1]*r[0][1]*r[0][1] + mtmp[2]*r[0][2]*r[0][2];
      dbuf[1] = mtmp[0]*r[0][0]*r[1][0] + mtmp[1]*r[0][1]*r[1][1] + mtmp[2]*r[0][2]*r[1][2];
      dbuf[2] = mtmp[0]*r[0][0]*r[2][0] + mtmp[1]*r[0][1]*r[2][1] + mtmp[2]*r[0][2]*r[2][2];
      dbuf[3] = mtmp[0]*r[1][0]*r[1][0] + mtmp[1]*r[1][1]*r[1][1] + mtmp[2]*r[1][2]*r[1][2];
      dbuf[4] = mtmp[0]*r[1][0]*r[2][0] + mtmp[1]*r[1][1]*r[2][1] + mtmp[2]*r[1][2]*r[2][2];
      dbuf[5] = mtmp[0]*r[2][0]*r[2][0] + mtmp[1]*r[2][1]*r[2][1] + mtmp[2]*r[2][2]*r[2][2];
    }
    else { // Cannot recover the metric
      for (i=0; i<sol->size; i++)  dbuf[i] = 0.;
    }
  }
  else {
    for (i=0; i<sol->size; i++)  dbuf[i] = sol->m[sol->size*ip+i];
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward an array of solutions.
 * \param filename name of file.
 * \param metricData 1 if the data saved is a metric (if only 1 data)
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of solutions at MSH  file format (.msh extension).
 * Write binary file for .mshb extension.and ASCII for .msh one.
 *
 */
int MMG5_saveMshMesh(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename,
                     int metricData) {
  FILE*       inm;
  MMG5_pPoint ppt;
  MMG5_pTetra pt;
  MMG5_pPrism pp;
  MMG5_pTria  ptt;
  MMG5_pQuad  pq;
  MMG5_pEdge  pa;
  MMG5_pSol   psl;
  double      dbuf[6];
  int         bin,k,i,typ,nelts,word, header[3],iadr;
  int         nq,ne,npr,np,nt,na,isol,nsols;
  char        *ptr,*data;
  static char mmgWarn = 0;

  bin = 0;

  MMG5_SAFE_CALLOC(data,strlen(filename)+7,char,return 0);
  strcpy(data,filename);

  ptr = strstr(data,".msh");
  if ( !ptr ) {
    /* data contains the filename without extension */
    strcat(data,".mshb");
    if (!(inm = fopen(data,"wb")) ) {
      ptr  = strstr(data,".msh");
      *ptr = '\0';
      strcat(data,".msh");
      if( !(inm = fopen(data,"wb")) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        MMG5_SAFE_FREE(data);
        return 0;
      }
    }
    else bin=1;
  }
  else {
    ptr = strstr(data,".mshb");
    if ( ptr ) bin = 1;
    if( !(inm = fopen(data,"wb")) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
      MMG5_SAFE_FREE(data);
      return 0;
    }
  }

  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

  /* Entete fichier*/
  fprintf(inm,"$MeshFormat\n");
  fprintf(inm,"2.2 %d %d\n",bin,8);
  if ( bin ) {
    word = 1;
    fwrite(&word,sw,1,inm);
    fprintf(inm,"\n");
  }
  fprintf(inm,"$EndMeshFormat\n");

  /* Vertices */
  np = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      ppt->tmp  = ++np;
      if ( mesh->dim==2 ) ppt->c[2] = 0.;
    }
  }
  fprintf(inm,"$Nodes\n");
  fprintf(inm,"%d\n",np);

  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) ) {
      if(!bin) {
        fprintf(inm," %d",ppt->tmp);
        for ( i=0; i<3; ++i )
          fprintf(inm," %.15lg",ppt->c[i]);
        fprintf(inm,"\n");
      } else {
        fwrite(&ppt->tmp,sw,1,inm);
        fwrite(&ppt->c[0],sd,3,inm);
      }
    }
  }
  if ( bin )  fprintf(inm,"\n");
  fprintf(inm,"$EndNodes\n");

  /* Data reading */
  /** First step: Count the number of elements of each type */
  ne = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    ++ne;
  }
  npr =  0;
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) ) continue;
    npr++;
  }
  nq = 0;
  for (k=1; k<=mesh->nquad; k++) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) )  continue;
    nq++;
  }
  nt = 0;
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MG_EOK(ptt) )  continue;
    nt++;
  }
  na = 0;
  for (k=1; k<=mesh->na; k++) {
    pa = &mesh->edge[k];
    if ( (!pa || !pa->a) )  continue;
    na++;
  }

  fprintf(inm,"$Elements\n");
  fprintf(inm,"%d\n", np + ne + npr + nt + nq + na );

  /** Second step: save the elements at following format:
      "idx type tagNumber tag0 tag1... v0_elt v1_elt..." */
  nelts = 0;

  /* Nodes */
  if ( bin ) {
    header[0] = 15;// Node keyword
    header[1] = np;
    header[2] = 2; // 2 tags per node
    fwrite(header, sw, 3, inm);
  }

  for ( k=1; k<= mesh->np; ++k)
  {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;
    ++nelts;

    if ( !bin ) fprintf(inm,"%d 15 2 %d %d %d\n",nelts,abs(ppt->ref),
                        abs(ppt->ref),ppt->tmp);
    else {
      fwrite(&nelts,sw,1,inm);
      word = abs(ppt->ref);
      fwrite(&word,sw,1,inm);
      fwrite(&word,sw,1,inm);
      fwrite(&ppt->tmp,sw,1,inm);
    }
  }

  /* Edges */
  if ( bin && na ) {
    header[0] = 1;// Edge keyword
    header[1] = na;
    header[2] = 2; // 2 tags per edge
    fwrite(header, sw, 3, inm);
  }

  for (k=1; k<=mesh->na; ++k) {
    pa = &mesh->edge[k];

    if ( !pa || !pa->a ) continue;
    ++nelts;

    if(!bin)
      fprintf(inm,"%d 1 2 %d %d %d %d\n",nelts,pa->ref,pa->ref,
              mesh->point[pa->a].tmp,mesh->point[pa->b].tmp);
    else {
      fwrite(&nelts,sw,1,inm);
      fwrite(&pa->ref,sw,1,inm);
      fwrite(&pa->ref,sw,1,inm);
      fwrite(&mesh->point[pa->a].tmp,sw,1,inm);
      fwrite(&mesh->point[pa->b].tmp,sw,1,inm);
    }
  }

  /* Triangles */
  if ( bin && nt ) {
    header[0] = 2;// Tria keyword
    header[1] = nt;
    header[2] = 2; // 2 tags per tria
    fwrite(header, sw, 3, inm);
  }

  for (k=1; k<=mesh->nt; ++k) {
    ptt = &mesh->tria[k];

    if ( !MG_EOK(ptt) ) continue;
    ++nelts;

    if(!bin)
      fprintf(inm,"%d 2 2 %d %d %d %d %d\n",nelts,ptt->ref,ptt->ref,
              mesh->point[ptt->v[0]].tmp,mesh->point[ptt->v[1]].tmp,
              mesh->point[ptt->v[2]].tmp);
    else {
      fwrite(&nelts,sw,1,inm);
      fwrite(&ptt->ref,sw,1,inm);
      fwrite(&ptt->ref,sw,1,inm);
      fwrite(&mesh->point[ptt->v[0]].tmp,sw,1,inm);
      fwrite(&mesh->point[ptt->v[1]].tmp,sw,1,inm);
      fwrite(&mesh->point[ptt->v[2]].tmp,sw,1,inm);
    }
  }

  /* Quads */
  if ( bin && nq ) {
    header[0] = 3;// Quad keyword
    header[1] = nq;
    header[2] = 2; // 2 tags per quad
    fwrite(header, sw, 3, inm);
  }

  for (k=1; k<=mesh->nquad; ++k) {
    pq = &mesh->quadra[k];
    if ( !MG_EOK(pq) ) continue;
    ++nelts;

    if(!bin)
      fprintf(inm,"%d 3 2 %d %d %d %d %d %d\n",nelts,pq->ref,pq->ref,
              mesh->point[pq->v[0]].tmp,mesh->point[pq->v[1]].tmp,
              mesh->point[pq->v[2]].tmp,mesh->point[pq->v[3]].tmp);
    else {
      fwrite(&nelts,sw,1,inm);
      fwrite(&pq->ref,sw,1,inm);
      fwrite(&pq->ref,sw,1,inm);
      fwrite(&mesh->point[pq->v[0]].tmp,sw,1,inm);
      fwrite(&mesh->point[pq->v[1]].tmp,sw,1,inm);
      fwrite(&mesh->point[pq->v[2]].tmp,sw,1,inm);
      fwrite(&mesh->point[pq->v[3]].tmp,sw,1,inm);
    }
  }

  /* Tetra */
  if ( bin && ne ) {
    header[0] = 4;// Tetra keyword
    header[1] = ne;
    header[2] = 2; // 2 tags per quad
    fwrite(header, sw, 3, inm);
  }

  for (k=1; k<=mesh->ne; ++k) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    ++nelts;

    if(!bin)
      fprintf(inm,"%d 4 2 %d %d %d %d %d %d\n",nelts,pt->ref,pt->ref,
              mesh->point[pt->v[0]].tmp,mesh->point[pt->v[1]].tmp,
              mesh->point[pt->v[2]].tmp,mesh->point[pt->v[3]].tmp);
    else {
      fwrite(&nelts,sw,1,inm);
      fwrite(&pt->ref,sw,1,inm);
      fwrite(&pt->ref,sw,1,inm);
      fwrite(&mesh->point[pt->v[0]].tmp,sw,1,inm);
      fwrite(&mesh->point[pt->v[1]].tmp,sw,1,inm);
      fwrite(&mesh->point[pt->v[2]].tmp,sw,1,inm);
      fwrite(&mesh->point[pt->v[3]].tmp,sw,1,inm);
    }
  }

  /* Prisms */
  if ( bin && npr ) {
    header[0] = 6;// Prism keyword
    header[1] = npr;
    header[2] = 2; // 2 tags per prism
    fwrite(header, sw, 3, inm);
  }

  for (k=1; k<=mesh->nprism; ++k) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) ) continue;
    ++nelts;

    if(!bin)
      fprintf(inm,"%d 6 2 %d %d %d %d %d %d %d %d\n",nelts,pp->ref,pp->ref,
              mesh->point[pp->v[0]].tmp,mesh->point[pp->v[1]].tmp,
              mesh->point[pp->v[2]].tmp,mesh->point[pp->v[3]].tmp,
              mesh->point[pp->v[4]].tmp,mesh->point[pp->v[5]].tmp);
    else {
      fwrite(&nelts,sw,1,inm);
      fwrite(&pp->ref,sw,1,inm);
      fwrite(&pp->ref,sw,1,inm);
      fwrite(&mesh->point[pp->v[0]].tmp,sw,1,inm);
      fwrite(&mesh->point[pp->v[1]].tmp,sw,1,inm);
      fwrite(&mesh->point[pp->v[2]].tmp,sw,1,inm);
      fwrite(&mesh->point[pp->v[3]].tmp,sw,1,inm);
      fwrite(&mesh->point[pp->v[4]].tmp,sw,1,inm);
      fwrite(&mesh->point[pp->v[5]].tmp,sw,1,inm);
    }
  }
  if ( bin )  fprintf(inm,"\n");
  fprintf(inm,"$EndElements\n");

  /* stats */
  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"     NUMBER OF VERTICES       %8d\n",np);
    fprintf(stdout,"     NUMBER OF TETRAHEDRA     %8d\n",ne);
    if ( npr )
      fprintf(stdout,"     NUMBER OF PRISMS         %8d\n",npr);
    if ( nt )
      fprintf(stdout,"     NUMBER OF TRIANGLES      %8d\n",nt);
    if ( nq )
      fprintf(stdout,"     NUMBER OF QUADRILATERALS %8d\n",nq);
    if ( na ) {
      fprintf(stdout,"     NUMBER OF EDGES          %8d\n",na);
    }

  }

  /** Write solution */
  nsols = (metricData==1)? 1 : mesh->nsols;

  for ( isol=0; isol<nsols; ++isol) {
    psl = *sol + isol;

    if ( !psl->m ) {
      if ( !mmgWarn ) {
        mmgWarn = 1;
        fprintf(stderr, "  ## Warning: %s: missing data for at least 1 solution."
                " Skipped.\n",__func__);
      }
      continue;
    }

    fprintf(inm,"$NodeData\n");

    /* One string tag saying the type of solution saved.
       Add the "metric" keyword if we save a metric */
    fprintf(inm,"1\n");

    if ( psl->size == 1 ) {
      typ = 1;
    }
    else if ( psl->size == psl->dim ) {
      typ = 3;
    }
    else {
      typ = 9;
    }
    if ( metricData ) {
      fprintf(inm,"\"%s:metric\"\n",psl->namein);
    }
    else {
      fprintf(inm,"\"%s\"\n",psl->namein);
    }

    /* One real tag unused */
    fprintf(inm,"1\n");
    fprintf(inm,"0.0\n"); // time value: unused

    /* Three integer tags */
    fprintf(inm,"3\n");
    fprintf(inm,"0\n"); // Time step: unused
    fprintf(inm,"%d\n",typ);
    fprintf(inm,"%d\n",np);

    /** Save the solution at following format:
        "idx sol" */
    if ( psl->size!= (psl->dim*(psl->dim+1))/2 ) {
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) ) continue;

        iadr = k*psl->size;
        if ( psl->dim == 2 )
          dbuf[2] = 0; // z-component for a vector field

        for ( i=0; i<psl->size; ++i )
          dbuf[i] = psl->m[iadr+i];

        if ( !bin ) {
          fprintf(inm,"%d",ppt->tmp);
          for ( i=0; i<typ; ++i )
            fprintf(inm," %lg",dbuf[i]);
          fprintf(inm,"\n");
        }
        else {
          fwrite(&ppt->tmp,sw,1,inm);
          fwrite(&dbuf,sd,typ,inm);
        }
      }
    }
    else {
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) ) continue;

        if ( psl->dim == 3 ) {
          if ( metricData ) {
            assert(!mesh->nsols);
            MMG5_build3DMetric(mesh,psl,k,dbuf);
          }
          else {
            for (i=0; i<psl->size; i++)  dbuf[i] = psl->m[psl->size*k+i];
          }
        }

        if(!bin) {
          fprintf(inm,"%d",ppt->tmp);
          if ( psl->dim==2 ) {
            iadr = k*psl->size;
            fprintf(inm," %.15lg %.15lg %.15lg %.15lg %.15lg %.15lg"
                    " %.15lg %.15lg %.15lg \n",
                    psl->m[iadr],psl->m[iadr+1],0.,psl->m[iadr+1],psl->m[iadr+2],0.,0.,0.,1.);
          }
          else {
            fprintf(inm," %.15lg %.15lg %.15lg %.15lg %.15lg %.15lg"
                    " %.15lg %.15lg %.15lg \n", dbuf[0],dbuf[1],dbuf[2],
                    dbuf[1],dbuf[3],dbuf[4],dbuf[2],dbuf[4],dbuf[5]);
          }
        }
        else {
          fwrite(&ppt->tmp,sw,1,inm);
          if ( psl->dim==2 ) {
            iadr = k*psl->size;
            fwrite(&psl->m[iadr],sd,2,inm);
            dbuf[0] = dbuf[1] = dbuf[2] = 0.;
            dbuf[3] = 1.;
            fwrite(&dbuf,sd,1,inm);
            fwrite(&psl->m[iadr+1],sd,2,inm);
            fwrite(&dbuf,sd,4,inm);
          }
          else {
            fwrite(&dbuf[0],sd,3,inm);
            fwrite(&dbuf[1],sd,1,inm);
            fwrite(&dbuf[3],sd,2,inm);
            fwrite(&dbuf[2],sd,1,inm);
            fwrite(&dbuf[4],sd,2,inm);
          }
        }
      }
    }
    if ( bin ) fprintf(inm,"\n");
    fprintf(inm,"$EndNodeData\n");
  }
  fclose(inm);

  return 1;
}

/**
 * \param filename name of file.
 * \param meshDim mesh dimenson.
 * \param inm allocatable pointer toward the FILE structure
 * \param ver file version (1=simple precision, 2=double)
 * \param bin 1 if the file is a binary
 * \param iswp 1 or 0 depending on the endianness (binary only)
 * \param np number of solutions of each type
 * \param dim solution dimension
 * \param nsols number of solutions of different types in the file
 * \param type type of solutions
 * \param posnp pointer toward the position of the point list in the file
 * \param imprim verbosity
 *
 * \return -1 data invalid or we fail, 0 no file, 1 ok.
 *
 * Open the "filename" solution file and read the file header.
 *
 */
int MMG5_loadSolHeader( const char *filename,int meshDim,FILE **inm,int *ver,
                        int *bin,int *iswp,int *np,int *dim,int *nsols,int **type,
                        long *posnp, int imprim) {
  int         binch,bdim;
  int         bpos,i;
  char        *ptr,*data,chaine[128];

  *posnp = 0;
  *bin   = 0;
  *iswp  = 0;
  *ver   = 0;

  MMG5_SAFE_CALLOC(data,strlen(filename)+6,char,return -1);
  strcpy(data,filename);

  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';

  ptr = strstr(data,".sol");

  if ( !ptr ) {
    /* data contains the filename without extension */
    strcat(data,".solb");
    if (!(*inm = fopen(data,"rb"))  ) {
      /* our file is not a .solb file, try with .sol ext */
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      if (!(*inm = fopen(data,"rb"))  ) {
        if ( imprim >= 0 )
          fprintf(stderr,"  ** %s  NOT FOUND. USE DEFAULT METRIC.\n",data);
        MMG5_SAFE_FREE(data);
        return 0;
      }
    } else {
      *bin = 1;
    }
  }
  else {
    ptr = strstr(data,".solb");
    if ( ptr )  *bin = 1;

    if (!(*inm = fopen(data,"rb")) ) {
      if ( imprim >= 0 )
        fprintf(stderr,"  ** %s  NOT FOUND. USE DEFAULT METRIC.\n",data);
      MMG5_SAFE_FREE(data);
      return 0;
    }
  }
  if ( imprim >= 0 )
    fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

  /* read solution or metric */
  if(!*bin) {
    strcpy(chaine,"DDD");
    while(fscanf(*inm,"%127s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if(!strncmp(chaine,"Dimension",strlen("Dimension"))) {
        fscanf(*inm,"%d",dim);
        if ( *dim!=meshDim ) {
          fprintf(stderr,"BAD SOL DIMENSION: %d\n",*dim);
          fclose(*inm);
          return -1;
        }
        continue;
      } else if(!strncmp(chaine,"SolAtVertices",strlen("SolAtVertices"))) {
        fscanf(*inm,"%d",np);
        fscanf(*inm,"%d",nsols);
        MMG5_SAFE_CALLOC(*type,*nsols,int,return -1);
        for ( i=0; i<*nsols; ++i ) {
          fscanf(*inm,"%d",&((*type)[i]));
        }
        *posnp = ftell(*inm);
        break;
      }
    }
  } else {
    fread(&binch,sw,1,*inm);
    if(binch==16777216) (*iswp)=1;
    else if(binch!=1) {
      fprintf(stdout,"BAD FILE ENCODING\n");
      fclose(*inm);
      return -1;
    }
    fread(ver,sw,1,*inm);
    if ( *iswp ) *ver = MMG5_swapbin(*ver);
    while(fread(&binch,sw,1,*inm)!=EOF && binch!=54 ) {
      if ( *iswp ) binch=MMG5_swapbin(binch);
      if(binch==54) break;
      if(binch==3) {  //Dimension
        fread(&bdim,sw,1,*inm);  //NulPos=>20
        if ( *iswp ) bdim=MMG5_swapbin(bdim);
        fread(dim,sw,1,*inm);
        if ( *iswp ) *dim=MMG5_swapbin(*dim);
        if ( *dim!=meshDim ) {
          fprintf(stderr,"BAD SOL DIMENSION: %d\n",*dim);
          printf("  Exit program.\n");
          fclose(*inm);
          return -1;
        }
        continue;
      } else if(binch==62) {  //SolAtVertices
        fread(&binch,sw,1,*inm); //Pos
        if ( *iswp ) binch=MMG5_swapbin(binch);
        fread(np,sw,1,*inm);
        if ( *iswp ) *np=MMG5_swapbin(*np);
        fread(nsols,sw,1,*inm); //nb sol
        if ( *iswp ) *nsols =MMG5_swapbin(*nsols);

        MMG5_SAFE_CALLOC(*type,*nsols,int,return -1); //typSol
        for ( i=0; i<*nsols; ++i ) {
          fread(&((*type)[i]),sw,1,*inm);
          if ( *iswp ) (*type)[i]=MMG5_swapbin((*type)[i]);
        }
        *posnp = ftell(*inm);
        break;
      } else {
        fread(&bpos,sw,1,*inm); //Pos
        if ( *iswp ) bpos=MMG5_swapbin(bpos);
        rewind(*inm);
        fseek(*inm,bpos,SEEK_SET);
      }
    }
  }

  return 1;
}

/**
 * \param sol pointer toward an allocatable sol structure.
 * \param inm pointer toward the solution file
 * \param bin 1 if binary file
 * \param iswp Endianess
 * \param index of the readed solution
 *
 * Read the solution value for vertex of index pos in floating precision.
 *
 */
void MMG5_readFloatSol3D(MMG5_pSol sol,FILE *inm,int bin,int iswp,int pos) {
  float       fbuf[6],tmpf;
  int         i;

  switch ( sol->size ) {
  case 1: case 3:
    /* scalar or vector solution */
    for (i=0; i<sol->size; i++) {
      if(!bin){
        fscanf(inm,"%f",&fbuf[0]);
      } else {
        fread(&fbuf[0],sw,1,inm);
        if(iswp) fbuf[0]=MMG5_swapf(fbuf[0]);
      }
      sol->m[sol->size*pos+i] = fbuf[0];
    }
    break;
  case 6 :
    /* Tensor solution */
    if(!bin){
      for(i=0 ; i<sol->size ; i++)
        fscanf(inm,"%f",&fbuf[i]);
    } else {
      for(i=0 ; i<sol->size ; i++) {
        fread(&fbuf[i],sw,1,inm);
        if(iswp) fbuf[i]=MMG5_swapf(fbuf[i]);
      }
    }
    tmpf    = fbuf[2];
    fbuf[2] = fbuf[3];
    fbuf[3] = tmpf;
    for (i=0; i<6; i++)  sol->m[6*pos+i] = fbuf[i];
    break;
  }
}

/**
 * \param sol pointer toward an allocatable sol structure.
 * \param inm pointer toward the solution file
 * \param bin 1 if binary file
 * \param iswp Endianess
 * \param index of the readed solution
 *
 * Read the solution value for vertex of index pos in double precision.
 *
 */
void MMG5_readDoubleSol3D(MMG5_pSol sol,FILE *inm,int bin,int iswp,int pos) {
  double      dbuf[6],tmpd;
  int         i;

  switch ( sol->size ) {
  case 1: case 3:
    /* scalar or vector solution */
    for (i=0; i<sol->size; i++) {
      if(!bin){
        fscanf(inm,"%lf",&dbuf[i]);
      } else {
        fread(&dbuf[i],sd,1,inm);
        if(iswp) dbuf[i]=MMG5_swapd(dbuf[i]);
      }
      sol->m[sol->size*pos+i] = dbuf[i];
    }
    break;

  case 6 :
    /* tensor solution */
    if(!bin){
      for(i=0 ; i<sol->size ; i++)
        fscanf(inm,"%lf",&dbuf[i]);
    } else {
      for(i=0 ; i<sol->size ; i++) {
        fread(&dbuf[i],sd,1,inm);
        if(iswp) dbuf[i]=MMG5_swapf(dbuf[i]);
      }
    }
    tmpd    = dbuf[2];
    dbuf[2] = dbuf[3];
    dbuf[3] = tmpd;
    for (i=0; i<sol->size; i++)  sol->m[6*pos+i] = dbuf[i];
    break;
  }
}

/**
 * \param mesh pointer toward the mesh structure
 * \param sol pointer toward an allocatable sol structure.
 * \param inm pointer toward the solution file
 * \param bin 1 if binary file
 * \param pos of the writted solution
 * \param metricData 1 if the data saved is a metric (if only 1 data)
 *
 * Write the solution value for vertex of index pos in double precision.
 *
 */
void MMG5_writeDoubleSol3D(MMG5_pMesh mesh,MMG5_pSol sol,FILE *inm,int bin,
                           int pos,int metricData) {
  double      dbuf[6],tmp;
  int         i;

  switch ( sol->size ) {
  case 1: case 3:
    /* scalar or vector solution */
    for (i=0; i<sol->size; i++) {
      for (i=0; i<sol->size; i++) dbuf[i] = sol->m[sol->size*pos+i];
      if(!bin){
        for (i=0; i<sol->size; i++)
          fprintf(inm," %.15lg",dbuf[i]);
      } else {
        for(i=0; i<sol->size; i++)
          fwrite((unsigned char*)&dbuf[i],sd,1,inm);
      }
    }
    break;

  case 6 :
    /* tensor solution */
    if ( metricData )
    {
      MMG5_build3DMetric(mesh,sol,pos,dbuf);
    }
    else
    {
      for (i=0; i<sol->size; i++)
      {
        dbuf[i] = sol->m[6*pos+i];
      }
    }

    tmp = dbuf[2];
    dbuf[2] = dbuf[3];
    dbuf[3] = tmp;

    if(!bin) {
      for(i=0; i<sol->size; i++)
        fprintf(inm," %.15lg",dbuf[i]);
    } else {
      for(i=0; i<sol->size; i++)
        fwrite((unsigned char*)&dbuf[i],sd,1,inm);
    }
    break;
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \param inm allocatable pointer toward the FILE structure.
 * \param ver file version (1=simple precision, 2=double).
 * \param bin 1 if the file is a binary.
 * \param np number of solutions of each type.
 * \param dim solution dimension.
 * \param nsols number of solutions of different types in the file.
 * \param type type of solutions.
 * \param size size of solutions.
 *
 * \return 0 if unable to open the file, 1 if success.
 *
 * Open the "filename" solution file and read the file header.
 *
 */
int MMG5_saveSolHeader( MMG5_pMesh mesh,const char *filename,
                        FILE **inm,int ver,int *bin,int np,int dim,
                        int nsols,int *type,int *size) {
  MMG5_pPoint ppt;
  int         binch,bpos;
  int         k;
  char        *ptr,*data,chaine[128];

  *bin = 0;

  MMG5_SAFE_CALLOC(data,strlen(filename)+6,char,return 0);
  strcpy(data,filename);
  ptr = strstr(data,".sol");
  if ( ptr ) {
    // filename contains the solution extension
    ptr = strstr(data,".solb");

    if ( ptr )  *bin = 1;

    if( !(*inm = fopen(data,"wb")) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
      MMG5_SAFE_FREE(data);
      return 0;
    }
  }
  else
  {
    // filename don't contains the solution extension
    ptr = strstr(data,".mesh");
    if ( ptr ) *ptr = '\0';

    strcat(data,".sol");
    if (!(*inm = fopen(data,"wb")) ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      if (!(*inm = fopen(data,"wb")) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        MMG5_SAFE_FREE(data);
        return 0;
      }
      else *bin = 1;
    }
  }

  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

  /*entete fichier*/
  binch=bpos=0;
  if(!*bin) {
    strcpy(&chaine[0],"MeshVersionFormatted\n");
    fprintf(*inm,"%s %d",chaine,ver);
    strcpy(&chaine[0],"\n\nDimension\n");
    fprintf(*inm,"%s %d",chaine,dim);
  } else {
    binch = 1; //MeshVersionFormatted
    fwrite(&binch,sw,1,*inm);
    binch = ver; //version
    fwrite(&binch,sw,1,*inm);
    binch = 3; //Dimension
    fwrite(&binch,sw,1,*inm);
    bpos = 20; //Pos
    fwrite(&bpos,sw,1,*inm);
    binch = dim;
    fwrite(&binch,sw,1,*inm);
  }

  np = 0;
  for (k=1; k<=mesh->np; k++) {
    ppt = &mesh->point[k];
    if ( MG_VOK(ppt) )  np++;
  }

  if(!*bin) {
    strcpy(&chaine[0],"\n\nSolAtVertices\n");
    fprintf(*inm,"%s",chaine);
    fprintf(*inm,"%d\n",np);
    fprintf(*inm,"%d",nsols);
    for (k=0; k<nsols; ++k )
      fprintf(*inm," %d",type[k]);
    fprintf(*inm,"\n");
  } else {
    binch = 62; //Vertices
    fwrite(&binch,sw,1,*inm);
    bpos += 16;

    for (k=0; k<nsols; ++k )
      bpos += 4 + (size[k]*ver)*4*np; //Pos
    fwrite(&bpos,sw,1,*inm);

    fwrite(&np,sw,1,*inm);
    fwrite(&nsols,sw,1,*inm);
    for (k=0; k<nsols; ++k )
      fwrite(&type[k],sw,1,*inm);
  }

  return 1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param type type of the metric
 * \param inm metric file
 * \return 1 if success, -1 if fail
 *
 * Check if the type of the metric is compatible with the remeshing mode.
 * If not, deallocate the type array and close the metric file.
 *
 */
int MMG5_chkMetricType(MMG5_pMesh mesh,int *type, FILE *inm) {

  /* 1: scalar solution (isotropic metric or ls function,
     2: vector field (displacement in Lagrangian mode),
     3: anisotropic metric */
  if ( mesh->info.lag == -1 ) {
    if ( type[0]!=1 && type[0]!=3) {
      fprintf(stderr,"  ** DATA TYPE IGNORED %d \n",type[0]);
      MMG5_SAFE_FREE(type);
      fclose(inm);
      return -1;
    }
  }
  else {
    if ( type[0] != 2 ) {
      fprintf(stderr,"  ** MISMATCH DATA TYPE FOR LAGRANGIAN MODE %d \n",
              type[0]);
      MMG5_SAFE_FREE(type);
      fclose(inm);
      return -1;
    }
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 *
 * print metric statistics
 *
 */
void MMG5_printMetStats(MMG5_pMesh mesh,MMG5_pSol met) {
  if ( abs(mesh->info.imprim) > 3 ) {
    if ( met->size == 1 )
      fprintf(stdout,"     NUMBER OF SCALAR VALUES %8d\n",met->np);
    else if ( met->size == 3 )
      fprintf(stdout,"     NUMBER OF VECTOR VALUES %8d\n",met->np);
    else
      fprintf(stdout,"     NUMBER OF TENSOR VALUES %8d\n",met->np);
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the solutions array.
 *
 * print solutions statistics
 *
 */
void MMG5_printSolStats(MMG5_pMesh mesh,MMG5_pSol *sol) {
  int j;

  if ( abs(mesh->info.imprim) > 3 ) {
    fprintf(stdout,"     NUMBER OF SOLUTIONS PER ENTITY %8d\n",mesh->nsols);
    fprintf(stdout,"     TYPE OF SOLUTIONS:\n          ");
    for ( j=0; j<mesh->nsols; ++j ) {
      if ( (*sol)[j].size == 1 )
        fprintf(stdout," SCALAR");
      else if ( (*sol)[j].size == 3 )
        fprintf(stdout," VECTOR");
      else
        fprintf(stdout," TENSOR");
    }
    fprintf(stdout,"\n");
  }
}
