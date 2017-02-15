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

#define SHORT_MAX    0x7f //127


int optlen_iso_bar(MMG5_pMesh mesh,MMG5_pSol sol,double declic,int base);

//return le facteur rendant le deplacement valide
// si > 1 --> on bouge sans soucis
// si < 1 --> on ne reussira pas a bouger sans croiser!
double MMG_maxdep(MMG5_pMesh mesh,MMG5_pSol sol) {
  printf("comment because of the merge needs\n");
  exit(EXIT_FAILURE);
  /* MMG5_pTria pt; */
  /* MMG5_pPoint ppt; */
  /* Displ pd; */
  /* double    c[3][2],aire,grad11,grad12,grad21,grad22,degrad; */
  /* double    u1x,u1y,u2x,u2y,u3x,u3y,lambda,lambdamin,mod,lambda1,lambda2;    */
  /* double    lambdamax,airenew,airemin; */
  /* double    a,b,ce,delta,xba,xca,yba,yca,dd;    */
  /* int k,i,imin;   */
  /* airemin = 1e20; */
  /* imin    = 0; */
  /* for (k=1; k<=mesh->nt; k++) { */
  /*   pt = &mesh->tria[k]; */
  /*   if ( !pt->v[0] )  continue; */
  /*   for (i=0; i<3; i++) {   */
  /*     ppt = &mesh->point[pt->v[i]]; */
  /*     memcpy(c[i],ppt->c,2*sizeof(double));   */
  /*   } */
  /*   aire = MMG2_quickarea(c[0],c[1],c[2]);; */
  /*   airemin = (airemin<aire)?airemin:aire; */
  /*   imin = (airemin<aire)?imin:k; */
  /* } */
  /* degrad = 0.9;//0.95; */
  /* printf("aire min : %e %e : %d %d %d\n",airemin,airemin*degrad,mesh->tria[imin].v[0], */
  /*   mesh->tria[imin].v[1],mesh->tria[imin].v[2]); */
  /* airemin*=degrad;  */
  /* pd = mesh->disp; */
  /* lambdamax = 1e20; */
  /* lambdamin = -1e20; */
  /* for (k=1; k<=mesh->nt; k++) { */
  /*   pt = &mesh->tria[k]; */
  /*   if ( !pt->v[0] )  continue; */
  /*   //if(k==5) exit(EXIT_FAILURE); */
  /*   //printf("triangle %d %d %d\n",pt->v[0],pt->v[1],pt->v[2]); */
  /*   mod = 0; */
  /*   for (i=0; i<3; i++) {   */
  /*     ppt = &mesh->point[pt->v[i]]; */
  /*     memcpy(c[i],ppt->c,2*sizeof(double));   */
  /*     mod += pd.mv[2*(pt->v[i]-1) + 1 + 0]*pd.mv[2*(pt->v[i]-1) + 1 + 0] + */
  /*  pd.mv[2*(pt->v[i]-1) + 1 + 1]*pd.mv[2*(pt->v[i]-1) + 1 + 1];  */
  /*   } */
  /*   // printf("mod %e\n",mod);   */
  /*   if(sqrt(mod) < 1e-24) continue; */
  /*   aire = MMG2_quickarea(c[0],c[1],c[2]);  */
  /*   dd = 1./sqrt(aire);    */
  /*   //printf("pour %d aire intiale %e -- normalisation by %e\n",k,aire,dd);  */
  /*   xba = dd*(c[1][0]-c[0][0]); */
  /*   yba = dd*(c[1][1]-c[0][1]); */
  /*   xca = dd*(c[2][0]-c[0][0]); */
  /*   yca = dd*(c[2][1]-c[0][1]); */
  /*   u1x = dd*pd.mv[2*(pt->v[0]-1) + 1 + 0]; */
  /*   u1y = dd*pd.mv[2*(pt->v[0]-1) + 1 + 1]; */
  /*   u2x = dd*pd.mv[2*(pt->v[1]-1) + 1 + 0]; */
  /*   u2y = dd*pd.mv[2*(pt->v[1]-1) + 1 + 1]; */
  /*   u3x = dd*pd.mv[2*(pt->v[2]-1) + 1 + 0]; */
  /*   u3y = dd*pd.mv[2*(pt->v[2]-1) + 1 + 1];    */

  /*   grad11 = (u2x-u1x); */
  /*   grad21 = (u3x-u1x); */
  /*   grad12 = (u2y-u1y);              */
  /*   grad22 = (u3y-u1y);     */
  /*   // printf("grad %e %e %e %e\n",grad11,grad21,grad12,grad22);  */

  /*   b = xba*grad22+yca*grad11-(yba*grad21+xca*grad12);   */
  /*   a = grad11*grad22-grad21*grad12;          */
  /*   ce = xba*yca-yba*xca-airemin*dd*dd; */
  /*   if(fabs(a)<1e-20) { */
  /*     //printf("poly de degre 1!!!!! derive %e -- %e %e\n",b,a,ce);    */
  /*     //printf("si lambda = 1 : %e\n",a+b+ce); */
  /*     lambda = -1./b;  */
  /*     if(b > 0) { */
  /*  //printf("positif si lambda > %e -- lambdamin %e\n",lambda,lambdamin);  */
  /*  assert(lambda < 0); */
  /*  lambda = lambdamax;   */
  /*     } else { */
  /*  //printf("positif si lambda < %e -- lambdamax %e\n",lambda,lambdamax);  */
  /*  if(lambdamax < lambda) { */
  /*    lambda = lambdamax; */
  /*  }        */
  /*     } */
  /*   } else{ */
  /*     delta = b*b-4*a*ce; */
  /*     //printf("resolution : %e\n",delta); */
  /*     if(delta < 0) {   */
  /*  //printf("pour %d aire %e -- ",k,aire); */
  /*  //printf("pas de solution!!!!!!!!!!! %e\n",delta); */
  /*  //printf("a b %e %e %e\n",a,b,ce);  */
  /*  //printf("grad %e %e %e %e --- lmax %e\n",grad11,grad21,grad12,grad22,lambdamax);  */
  /*  airenew = 1./(dd*dd)*(a+b+ce); */
  /*  //printf("avec notre methode on trouve %e\n",airenew); */
  /*  for (i=0; i<3; i++) {   */
  /*    ppt = &mesh->point[pt->v[i]];   */
  /*    c[i][0] = ppt->c[0] + pd.mv[2*(pt->v[i]-1) + 1 + 0]; */
  /*    c[i][1] = ppt->c[1] + pd.mv[2*(pt->v[i]-1) + 1 + 1];  */
  /*  }      */
  /*  //printf("le triangle marche ?? %e\n",MMG2_quickarea(c[0],c[1],c[2]));      */
  /*  if(airenew > 0) continue;  */
  /*  printf("delta negatif on sait pas quoi faire!!!!!\n"); */
  /*  exit(EXIT_FAILURE); */
  /*     }           */
  /*     delta = sqrt(delta); */

  /*     lambda1 = (-b - delta)/(2.*a); */
  /*     lambda2 = (-b + delta)/(2.*a); */
  /*     //exit(EXIT_FAILURE);  */
  /*     lambda = 0.5*(lambda1+lambda2);  */
  /*     // printf("lambdamilieu %e -- %e et %e\n",lambda,1+lambda*b+lambda*lambda*a,1+b+a); */
  /*     if((ce+lambda*b+lambda*lambda*a) > 0) { */
  /*  lambda = ((lambda1) > (lambda2)) ? (lambda1):(lambda2); */
  /*     } else { */
  /*  // printf("eh eh autre cas!!!!\n"); */
  /*  // printf("  on trouve lambda %e %e ---- lambdamax %e\n",aire*lambda1,aire*lambda2,lambdamax);  */
  /*  lambda = ((lambda1) > (lambda2)) ? (lambda1):(lambda2);  */
  /*  if(lambda < 0) { */
  /*    lambda = lambdamax; */
  /*  } else { */
  /*    lambdamin = (lambdamin > lambda)? lambdamin : lambda; */
  /*    lambda1 = ((lambda1) < (lambda2)) ? (lambda1):(lambda2); */
  /*    if(!(lambdamax<lambda1)) { */
  /*      lambda = lambda1; */
  /*    } */
  /*  } */
  /*     }   */
  /*   } */
  /*   //     if(lambda < 1) {   */
  /*   //   airenew = ce+b+a; */
  /*   //   //printf("lambda < 1!!! avec notre methode on trouve %e\n",airenew); */
  /*   //   airenew = (xba+u2x-u1x)*(yca+u3y-u1y)-((yba+u2y-u1y)*(xca+u3x-u1x)); */
  /*   //   //printf("2e methode --  %e\n",airenew); */
  /*   //   for (i=0; i<3; i++) {   */
  /*   //     ppt = &mesh->point[pt->v[i]];   */
  /*   //     c[i][0] = ppt->c[0] + pd.mv[2*(pt->v[i]-1) + 1 + 0]; */
  /*   //         c[i][1] = ppt->c[1] + pd.mv[2*(pt->v[i]-1) + 1 + 1];  */
  /*   //       }      */
  /*   //   //printf("le triangle marche ?? %e\n",MMG2_quickarea(c[0],c[1],c[2]));   */
  /*   //   //printf("  on trouve lambda %e %e ---- lambdamax %e\n",lambda1,lambda2,lambdamax);  */
  /*   //   //exit(EXIT_FAILURE);   */
  /*   // } */
  /*   lambdamax = (lambdamax < (lambda)) ? lambdamax:(lambda);     */
  /*   if(lambdamax<0) {printf("heu lambdamax < 0!!! %e tr %d -- dep %e\n",lambdamax,k,mod);exit(EXIT_FAILURE);   } */
  /*   //printf("  on trouve lambda %e %e ---- lambdamax %e\n",lambda1,lambda2,lambdamax);  */
  /* }   */
  /* printf("%e < dep or dep < %e\n",lambdamin,lambdamax); */
  /* //assert(lambdamin > lambdamax); */
  //  return(lambdamax);
}
/* dichotomy: check if nodes can move */
int MMG_dikomv(MMG5_pMesh mesh,MMG5_pSol sol,short t) {
  printf("comment because of the merge needs\n");
  exit(EXIT_FAILURE);
  /*  MMG5_pTria     pt; */
/*   MMG5_pPoint    ppt; */
/*   Displ     pd; */
/*   double    c[3][2],alpha,aire; */
/*   int       k,i,nm; */

/*   pd = mesh->disp; */

/*   alpha = (double) t / SHORT_MAX;   */
/*   for (k=1; k<=mesh->nt; k++) { */
/*     pt = &mesh->tria[k]; */
/*     if ( !pt->v[0] )  continue; */
/*     for (i=0; i<3; i++) { */
/*       ppt      = &mesh->point[ pt->v[i] ]; */
/*       ppt->tmp = k; */
/*       if ( ppt->tag & M_MOVE ) { */
/*         c[i][0] = ppt->c[0] + alpha * pd.mv[2*(pt->v[i]-1) + 1 + 0]; */
/*         c[i][1] = ppt->c[1] + alpha * pd.mv[2*(pt->v[i]-1) + 1 + 1];  */
/*       } */
/*       else */
/*         memcpy(c[i],ppt->c,2*sizeof(double)); */
/*     } */

/*     aire = MMG2_quickarea(c[0],c[1],c[2]); */
/*     if ( aire < 1e-24 )  { */
/*       // if(mesh->info.imprim < 0)  */
/* printf("aire reject %d %e %e\n",k,aire,pt->qual * ALPHA); */
/*       return(0); */
/*     } */
/*   } */

/*   /\* update metrics *\/ */

/*   /\* update point coords *\/ */
/*   nm = 0; */
/*   for (k=1; k<=mesh->np; k++) { */
/*     ppt = &mesh->point[k]; */
/*     if ( ppt->tag & M_MOVE ) { */
/*       ppt->c[0] += alpha * pd.mv[2*(k-1) + 1 + 0]; */
/*       ppt->c[1] += alpha * pd.mv[2*(k-1) + 1 + 1]; */
/*       pd.alpha[k]  = t; */
/*       if ( t == SHORT_MAX )  ppt->tag &= ~M_MOVE;        */
/*       nm++; */
/*     } */
/*   } */

/*   /\*MAJ qual*\/ */
/*   for (k=1; k<=mesh->nt; k++) { */
/*     pt = &mesh->tria[k]; */
/*  if ( !pt->v[0] ) continue;  */
/*     pt->qual = MMG2_caltri_in(mesh,sol,pt); */
/*   } */

/*   if ( mesh->info.imprim < 0 )  fprintf(stdout,"     %7d NODES UPDATED\n",nm); */
/*   return(nm); */
}


/* check if displacement ok */
int MMG_chkmov(MMG5_pMesh mesh,char level) {
  printf("comment because of the merge needs\n");
  exit(EXIT_FAILURE);
  /* MMG5_pTria      pt; */
  /* MMG5_pPoint     ppt; */
  /* // Displ      pd; */
  /* int        k,nc; */

  /* pd  = mesh->disp; */

  /* nc = 0; */
  /* for (k=1; k<=mesh->np; k++) { */
  /*   ppt = &mesh->point[k]; */
  /*   if ( ppt->tag & M_MOVE ) { */
  /*     if ( pd.alpha[k] )  return(pd.alpha[k]); */
  /*     ppt->tag &= ~M_MOVE; */
  /*     nc++; */
  /*   } */
  /* } */

  /* /\* check element validity *\/ */
  /* /\*if ( level > 0 ) { */
  /*   for (k=1 ; k<=mesh->ne; k++) { */
  /*     pt = &mesh->tetra[k]; */
  /*     if ( !pt->v[0] )  continue; */
  /*     vol = MMG_voltet(mesh,k); */
  /*     if ( vol < 0.0 )  return(0); */
  /*   } */
  /* }*\/ */

  /* return(0); */
}





/*rigid bodies moving*/
int MMG2_mmg2d9(MMG5_pMesh mesh,MMG5_pSol sol) {
  printf("comment because of the merge needs\n");
  exit(EXIT_FAILURE);
  /* MMG5_pPoint   ppt;    */
  /* MMG5_pTria    pt; */
  /* Displ    pd; */
  /* double   d1,declic,dd,qworstbef,qworst,qavg,qavgbef; */
  /* int      iter,maxiter,it,maxtou,ntreal; */
  /* int      iold,k,base,ns,nsiter,nm,nmiter,nmbar;   */
  /*  double lambda; */
  /* short    t,i,alpha; */

  /* if ( mesh->info.imprim < 0 ) { */
  /*   MMG2_outqua(mesh,sol); */
  /*   MMG2_prilen(mesh,sol); */
  /* } */
  /* /\* normalize coordinates *\/ */
  /* dd = 1;//PRECI / (mesh->info).delta;  */
  /* printf("dd %e\n",dd); */
  /* pd  = mesh->disp; */
  /* for (k=1; k<=mesh->np; k++) { */
  /*   ppt = &mesh->point[k]; */
  /*   if ( !M_VOK(ppt) )  continue; */
  /*   pd.mv[2*(k-1) + 1 + 0] *= /\*20**\/dd; */
  /*   pd.mv[2*(k-1) + 1 + 1] *= /\*20**\/dd; */
  /*   d1 = pd.mv[2*(k-1) + 1 + 0]*pd.mv[2*(k-1) + 1 + 0] */
  /*     + pd.mv[2*(k-1) + 1 + 1]*pd.mv[2*(k-1) + 1 + 1];    */
  /*   if(k==140561) printf("dep %e %e %e\n",d1,pd.mv[2*(k-1) + 1 + 0],pd.mv[2*(k-1) + 1 + 1]); */
  /*   if ( d1 > 1e-24 )  ppt->tag  |= M_MOVE; */
  /* } */
  /* lambda =  MMG_maxdep(mesh,sol); */
  /* printf("le dep max est %e\n",lambda);  */
  /* if((mesh->info.option)==99) {       */
  /*   for (k=1; k<=mesh->np; k++) { */
  /*     ppt = &mesh->point[k]; */
  /*     if ( !M_VOK(ppt) )  continue; */
  /*     pd.mv[2*(k-1) + 1 + 0] /= dd; */
  /*     pd.mv[2*(k-1) + 1 + 1] /= dd; */
  /*   } */
  /*   if(!MMG2D_saveVect(mesh,sol,mesh->namedep,lambda)) return(0); */
  /*   return(1); */
  /* } */
  /* /\*seuil declenchement du post-traitement : qdegradbef qworstbef*\/ */
  /* qworstbef = 1.;  */
  /* qavgbef   = 0.;  */
  /* ntreal = 0; */
  /* for(k=1 ; k<=mesh->nt ; k++) { */
  /*   pt = &mesh->tria[k];  */
  /*   if(!pt->v[0] )  continue; */
  /*   if(pt->qual > qworstbef) qworstbef = pt->qual;     */
  /*   ntreal++; */
  /*   qavgbef += pt->qual; */

  /* }                        */
  /* qavgbef /= (double) ntreal; */



  /* /\* move grid nodes *\/ */
  /* t       = SHORT_MAX; */
  /* alpha   = 0; */
  /* iter    = 0; */
  /* maxiter = 110; */
  /* iold    = 1; */
  /* ns      = 0; */

  /* /\* move grid nodes *\/ */
  /* t = SHORT_MAX; */
  /* if (  MMG_dikomv(mesh,sol,t) ) { */
  /*   if ( mesh->info.imprim )  fprintf(stdout,"     %7d NODES MOVED\n",mesh->np); */
  /* } */
  /* else { */
  /*   if ( mesh->info.imprim < 0) fprintf(stdout,"     TRYING DICHO\n"); */
  /*   while (t && alpha < SHORT_MAX && iter++ < maxiter) {   */
  /*     if ( mesh->info.imprim < 0)  */
  /*  printf("ITER %d  : alpha %d / %d\n",iter,alpha,SHORT_MAX);   */

  /*     /\*optim*\/  */
  /*     ns     = 0; */
  /*     nm     = 0; */
  /*     nmiter = 0; */
  /*     nsiter = 0;  */
  /*     it     = 0; */
  /*     maxtou = 0; */
  /*     do { */
  /*       /\*adaptation : insertion/collapse*\/ */
  /*       if(!mesh->info.noinsert) { */
  /*         fprintf(stdout,"\n  -- LENGTH ANALYSIS\n"); */
  /*         if ( !MMG2_mmg2d1(mesh,sol) )  { */
  /*      return(0);  */
  /*    }   */
  /*       } */
  /*       /\*edge flip*\/   */
  /*       if(!mesh->info.noswap) { */
  /*         declic = 1.1 / ALPHA; */
  /*         nsiter = MMG2_cendel(mesh,sol,declic,-1); */
  /*         if ( nsiter && mesh->info.imprim < 0) */
  /*      fprintf(stdout,"     %7d SWAPPED\n",nsiter);  */

  /*    ns+=nsiter;     */
  /*       } */
  /*       /\*point relocation*\/ */
  /*       if(!mesh->info.nomove) {  */
  /*         declic = 1.1 / ALPHA; */
  /*         base   = mesh->flag; */
  /*         nmiter = MMG2_optlen(mesh,sol,declic,-1); */
  /*         nmbar  = optlen_iso_bar(mesh,sol,declic,-1); */
  /*         nm += nmiter + nmbar; */
  /*         if ( mesh->info.imprim < 0) */
  /*           fprintf(stdout,"     %7d + %7d MOVED\n",nmiter,nmbar); */
  /*       }      */
  /*     } while( (nmiter+nsiter > 0) && (++it <= maxtou) ); */
  /*     if ( mesh->info.imprim ) */
  /*       fprintf(stdout,"     %7d SWAPPED %7d MOVED\n",ns,nm);     */


  /*     t = SHORT_MAX - alpha; */
  /*     i = 0; */
  /*     do {          */
  /*       nm = MMG_dikomv(mesh,sol,t); */
  /*       if ( nm )  { */
  /*         fprintf(stdout,"              ---- MOVE %d (%d)\n",t,i);                  */
  /*    alpha += t; */
  /*    break; */
  /*  } */
  /*  t = t >> 1;    */

  /*     } while (t && (++i < 8)); //en combien est-ce qu'on accepte de decouper */
  /*   } */
  /* } */

  /* /\* check mesh *\/  */
  /* if(alpha < SHORT_MAX) { */
  /*   if ( MMG_chkmov(mesh,1) ) { */
  /*     fprintf(stdout,"  ## UNCOMPLETE DISPLACEMENT (%d / %d)\n",alpha,SHORT_MAX); */
  /*     return(0); */
  /*   } */
  /* }  */

  /* /\*seuil declenchement du post-traitement : qdegrad qworst*\/ */
  /* qworst = 1.;  */
  /* qavg   = 0.;  */
  /* ntreal = 0; */
  /* for(k=1 ; k<=mesh->nt ; k++) { */
  /*   pt = &mesh->tria[k];  */
  /*   if(!pt->v[0] )  continue; */
  /*   if(pt->qual > qworst) qworst = pt->qual;     */
  /*   ntreal++; */
  /*   qavg += pt->qual; */
  /* }                        */
  /* qavg /= (double) ntreal; */

  /* if(abs(mesh->info.imprim) > 3) { */
  /*   fprintf(stdout,"\n     AVERAGE QUALITY %8f  ---> %8f  (%8f >? %8f)\n",qavgbef*ALPHA,qavg*ALPHA,qavg/qavgbef,mesh->info.qdegrad[1]);   */
  /*   fprintf(stdout,"     WORST QUALITY   %8f ---> %8f >? %8f\n",qworstbef*ALPHA,qworst*ALPHA,mesh->info.qdegrad[0]*ALPHA);   */
  /* } */

  /* if( (qworst < mesh->info.qdegrad[0]) && (qavg < qavgbef*mesh->info.qdegrad[1]) ) { */
  /*   /\*point relocation*\/ */
  /*   if(!mesh->info.nomove) {  */
  /*     fprintf(stdout,"\n  -- MESH OPTIMISATION\n"); */
  /*     nm     = 0; */
  /*     nmiter = 0; */
  /*     it     = 0; */
  /*     maxtou = 10; */
  /*     do { */
  /*       declic = 1.1 / ALPHA; */
  /*       base   = mesh->flag; */
  /*       nmiter = MMG2_optlen(mesh,sol,declic,-1); */
  /*       nmbar =  optlen_iso_bar(mesh,sol,declic,-1); */
  /*       nm += nmiter+nmbar; */
  /*       if ( mesh->info.imprim < 0) */
  /*         fprintf(stdout,"     %7d + %7d MOVED \n",nmiter,nmbar);   */

  /*     } while((nmiter+nsiter > 0) && (++it <= maxtou)); */

  /*     if ( mesh->info.imprim ) */
  /*  fprintf(stdout,"     %7d SWAPPED %7d MOVED\n",ns,nm); */
  /*   } */
  /*   return(1); */
  /* } */

  /* /\*adaptation : insertion/collapse*\/ */
  /* if(!mesh->info.noinsert) { */
  /*   fprintf(stdout,"\n  -- LENGTH ANALYSIS\n"); */
  /*   if ( !MMG2_mmg2d1(mesh,sol) )  { */

  /*     return(0); */
  /*   }    */
  /* } */

  /* /\*optim*\/  */
  /* fprintf(stdout,"\n  -- MESH OPTIMISATION\n"); */
  /* ns     = 0; */
  /* nm     = 0; */
  /* nmiter = 0; */
  /* nsiter = 0;  */
  /* it     = 0; */
  /* maxtou = 10; */
  /* do { */
  /*   /\*edge flip*\/   */
  /*   if(!mesh->info.noswap) { */
  /*     declic = 1.1 / ALPHA; */
  /*     nsiter = MMG2_cendel(mesh,sol,declic,-1); */
  /*     if ( nsiter && mesh->info.imprim < 0) */
  /*       fprintf(stdout,"     %7d SWAPPED\n",nsiter);  */

  /*     ns+=nsiter;     */
  /*   } */
  /*   /\*point relocation*\/ */
  /*   if(!mesh->info.nomove) {  */
  /*     declic = 1.1 / ALPHA; */
  /*     base   = mesh->flag; */
  /*     nmiter = MMG2_optlen(mesh,sol,declic,-1); */
  /*     nmbar =  optlen_iso_bar(mesh,sol,declic,-1); */
  /*     nm += nmiter+nmbar; */
  /*     if ( mesh->info.imprim < 0) */
  /*       fprintf(stdout,"     %7d + %7d MOVED \n",nmiter,nmbar); */
  /*   } */


  /* } while((nmiter+nsiter > 0) && (++it <= maxtou)); */
  /* if ( mesh->info.imprim ) */
  /*   fprintf(stdout,"     %7d SWAPPED %7d MOVED\n",ns,nm);  */


  return(1);
}
