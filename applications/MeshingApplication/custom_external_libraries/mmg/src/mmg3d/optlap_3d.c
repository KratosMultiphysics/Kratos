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
 * \file mmg3d/optlap_3d.c
 * \brief Functions for the optimization with laplacian/anti-laplacian.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "inlined_functions_3d.h"
/**
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the sol structure
 * \return 0 if fail, 1 otherwise.
 *
 **/

int _MMG3D_optlap(MMG5_pMesh mesh,MMG5_pSol sol) {
  MMG5_pTetra    pt,pt1;
  MMG5_pPoint    ppt,pptb,ppta;

  int       it,i,k,lon,l,iel,ipt,list[MMG3D_LMAX+2];
  int       maxiter,ipta,iptb,ipt0,ipt1,ipt2,ipt3,*compt;
  double    vol,ax,ay,az,bx,by,bz;
  double    *nv,*pos,res,dd,ox,oy,oz,declic;
  double LLAMBDA  = 0.33;
  double LMU      = 0.331;

  maxiter = 3;
  _MMG5_ADD_MEM(mesh,(3*mesh->np+1)*sizeof(double),"nv",
                return(0));
  _MMG5_ADD_MEM(mesh,(3*mesh->np+1)*sizeof(double),"pos",
                return(0));
  _MMG5_ADD_MEM(mesh,(mesh->np+1)*sizeof(double),"compt",
                return(0));
  _MMG5_SAFE_CALLOC(nv, 3*mesh->np+1, double);
  _MMG5_SAFE_CALLOC(pos, 3*mesh->np+1, double);
  _MMG5_SAFE_CALLOC(compt, mesh->np+1, int);

  it  = 1;
  declic = 3./_MMG5_ALPHAD;
  do {

    /*initialisation*/
    for(i = 1 ; i<=mesh->np ; i++) {
      ppt          = &mesh->point[i];
      compt[i] = 0;
      pos[3*(i-1) + 1 + 0] = 0.;
      pos[3*(i-1) + 1 + 1] = 0.;
      pos[3*(i-1) + 1 + 2] = 0.;
    }

    /*1st stage : laplacian*/
    for(k = 1 ; k<=mesh->ne ; k++) {
      pt = &mesh->tetra[k];
      if (!pt->v[0]) continue;
      if (pt->qual > declic) continue;

      for(i=0 ; i<6 ; i++) {
        ipta   = pt->v[_MMG5_iare[i][0]];
        ppta   = &mesh->point[ipta];
        
        iptb   = pt->v[_MMG5_iare[i][1]];
        pptb   = &mesh->point[iptb];
        
        if(!(ppta->tag & MG_BDY)) {
          pos[3*(ipta-1) + 1 + 0] += pptb->c[0];
          pos[3*(ipta-1) + 1 + 1] += pptb->c[1];
          pos[3*(ipta-1) + 1 + 2] += pptb->c[2];
          compt[ipta]++;
        }
        if(!(pptb->tag & MG_BDY)) {
          pos[3*(iptb-1) + 1 + 0] += ppta->c[0];
          pos[3*(iptb-1) + 1 + 1] += ppta->c[1];
          pos[3*(iptb-1) + 1 + 2] += ppta->c[2];
          compt[iptb]++;
        }
      }
    }
    
    for(i=1 ; i<=mesh->np ; i++) {
      ppt           = &mesh->point[i];
      if(compt[i]) {
        dd            = 1./(double) compt[i];
        pos[3*(i-1) + 1 + 0] *= dd;
        pos[3*(i-1) + 1 + 1] *= dd;
        pos[3*(i-1) + 1 + 2] *= dd;
        nv[3*(i-1) + 1] = ppt->c[0] + LLAMBDA * (ppt->c[0] - pos[3*(i-1) + 1 + 0]);
        nv[3*(i-1) + 2] = ppt->c[1] + LLAMBDA * (ppt->c[1] - pos[3*(i-1) + 1 + 1]);
        nv[3*(i-1) + 3] = ppt->c[2] + LLAMBDA * (ppt->c[2] - pos[3*(i-1) + 1 + 2]);
      } else {
        nv[3*(i-1) + 1] = ppt->c[0];
        nv[3*(i-1) + 2] = ppt->c[1];
        nv[3*(i-1) + 3] = ppt->c[2];
        
      }
      compt[i] = 0;
      pos[3*(i-1) + 1 + 0] = 0.;
      pos[3*(i-1) + 1 + 1] = 0.;
      pos[3*(i-1) + 1 + 2] = 0.;
      
    }
    
    /*2nd stage : anti-laplacian*/
    for(k = 1 ; k<=mesh->ne ; k++) {
      pt = &mesh->tetra[k];
      if (!pt->v[0]) continue;
      if (pt->qual > declic) continue;
      
      for(i=0 ; i<6 ; i++) {
        ipta   = pt->v[_MMG5_iare[i][0]];
        ppta   = &mesh->point[ipta];
        
        iptb   = pt->v[_MMG5_iare[i][1]];
        pptb   = &mesh->point[iptb];
        
        if(!(ppta->tag & MG_BDY)) {
          pos[3*(ipta-1) + 1 + 0] += nv[3*(iptb-1) + 1];
          pos[3*(ipta-1) + 1 + 1] += nv[3*(iptb-1) + 2];
          pos[3*(ipta-1) + 1 + 2] += nv[3*(iptb-1) + 3];
          compt[ipta]++;
        }
        if(!(pptb->tag & MG_BDY)) {
          pos[3*(iptb-1) + 1 + 0] += nv[3*(ipta-1) + 1];
          pos[3*(iptb-1) + 1 + 1] += nv[3*(ipta-1) + 2];
          pos[3*(iptb-1) + 1 + 2] += nv[3*(ipta-1) + 3];
          compt[iptb]++;
        }
      }
    }
    
    res= 0.;
    for(i=1 ; i<=mesh->np ; i++) {
      ppt           = &mesh->point[i];
      if(compt[i]) {
        dd            = 1./(double) compt[i];
        pos[3*(i-1) + 1 + 0] *= dd;
        pos[3*(i-1) + 1 + 1] *= dd;
        pos[3*(i-1) + 1 + 2] *= dd;
        ox = nv[3*(i-1) + 1];
        oy = nv[3*(i-1) + 2];
        oz = nv[3*(i-1) + 3];
        nv[3*(i-1) + 1] = nv[3*(i-1) + 1] - LMU * (nv[3*(i-1) + 1] - pos[3*(i-1) + 1 + 0]);
        nv[3*(i-1) + 2] = nv[3*(i-1) + 2] - LMU * (nv[3*(i-1) + 2] - pos[3*(i-1) + 1 + 1]);
        nv[3*(i-1) + 3] = nv[3*(i-1) + 3] - LMU * (nv[3*(i-1) + 3] - pos[3*(i-1) + 1 + 2]);
        
        dd = (nv[3*(i-1) + 1]-ox)*(nv[3*(i-1) + 1]-ox)
          + (nv[3*(i-1) + 2]-oy)*(nv[3*(i-1) + 2]-oy)
          + (nv[3*(i-1) + 3]-oz)*(nv[3*(i-1) + 3]-oz);
        res +=dd;
        
      }
      
      
      compt[i] = 0;
      pos[3*(i-1) + 1 + 0] = 0.;
      pos[3*(i-1) + 1 + 1] = 0.;
      pos[3*(i-1) + 1 + 2] = 0.;
    }

    /*check new coor*/
    for(k = 1 ; k<=mesh->ne ; k++) {
      pt = &mesh->tetra[k];
      if(!pt->v[0]) continue;
      
      for(i=0 ; i<4 ; i++) {
        ipt   = pt->v[i];
        ppt   = &mesh->point[ipt];
        if(ppt->tag & MG_BDY) continue;
        //if(ppt->tmp) continue;
        //ppt->tmp = 1;
        lon =_MMG5_boulevolp(mesh,k,i,&list[0]);
        for (l=0; l<lon; l++) {
          iel    = list[l] /4;
          pt1    = &mesh->tetra[iel];
          ipt0   = 3*(pt1->v[0] - 1);
          ipt1   = 3*(pt1->v[1] - 1);
          ipt2   = 3*(pt1->v[2] - 1);
          ipt3   = 3*(pt1->v[3] - 1);
          
          ax = nv[ipt2 + 1] - nv[ipt0 + 1];
          ay = nv[ipt2 + 2] - nv[ipt0 + 2];
          az = nv[ipt2 + 3] - nv[ipt0 + 3];
          
          bx = nv[ipt3 + 1] - nv[ipt0 + 1];
          by = nv[ipt3 + 2] - nv[ipt0 + 2];
          bz = nv[ipt3 + 3] - nv[ipt0 + 3];
          
          vol = (nv[ipt1 + 1] - nv[ipt0 + 1]) * (ay*bz - az*by) \
            + (nv[ipt1 + 2] - nv[ipt0 + 2]) * (az*bx - ax*bz)   \
            + (nv[ipt1 + 3] - nv[ipt0 + 3]) * (ax*by - ay*bx);
          if(vol < 0) {/*printf("reject1 %e\n",vol);*/break;}
        }
        if(l<=lon) {
          memcpy(&pos[3*(ipt-1) + 1],ppt->c,3*sizeof(double));
          for (l=0; l<lon; l++) {
            iel    = list[l] / 4;
            pt1    = &mesh->tetra[iel];
            ipt0   = 3*(pt1->v[0] - 1);
            ipt1   = 3*(pt1->v[1] - 1);
            ipt2   = 3*(pt1->v[2] - 1);
            ipt3   = 3*(pt1->v[3] - 1);
            
            ax = nv[ipt2 + 1] - nv[ipt0 + 1];
            ay = nv[ipt2 + 2] - nv[ipt0 + 2];
            az = nv[ipt2 + 3] - nv[ipt0 + 3];
            
            bx = nv[ipt3 + 1] - nv[ipt0 + 1];
            by = nv[ipt3 + 2] - nv[ipt0 + 2];
            bz = nv[ipt3 + 3] - nv[ipt0 + 3];
            
            vol = (nv[ipt1 + 1] - nv[ipt0 + 1]) * (ay*bz - az*by) \
              + (nv[ipt1 + 2] - nv[ipt0 + 2]) * (az*bx - ax*bz)   \
              + (nv[ipt1 + 3] - nv[ipt0 + 3]) * (ax*by - ay*bx);
            if(vol < 0) {/*printf("reject %e\n",vol);*/break;}
          }
          if(l<lon) break;
        }
      }
      if(i<4) break;
    }
    if(k > mesh->ne) {
      /*update coor*/
      for(i=1 ; i<=mesh->np ; i++) {
        ppt   = &mesh->point[i];
        ppt->c[0] = nv[3*(i-1) + 1];
        ppt->c[1] = nv[3*(i-1) + 2];
        ppt->c[2] = nv[3*(i-1) + 3];
      }
      for(k=1 ; k<=mesh->ne ; k++) {
        pt = &mesh->tetra[k];
        if(!pt->v[0]) continue;
        pt->qual = _MMG5_caltet(mesh,sol,pt);
      }
      if( mesh->info.imprim > 5) fprintf(stdout,"              LAPLACIAN : %8f\n",res);
    } else {
      if( mesh->info.imprim > 5) fprintf(stdout,"              NO LAPLACIAN\n");
      break;
    }
    if(res<1e-5) break;
    
  } while(it++ < maxiter);
  
  _MMG5_DEL_MEM(mesh,nv,(3*mesh->np+1)*sizeof(double));
  _MMG5_DEL_MEM(mesh,pos,(3*mesh->np+1)*sizeof(double));
  _MMG5_DEL_MEM(mesh,compt,(mesh->np+1)*sizeof(double));
  return(1);
}
