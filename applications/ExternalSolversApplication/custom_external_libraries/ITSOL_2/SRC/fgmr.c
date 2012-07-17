#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globheads.h"
#include "protos.h"

#define  epsmac  1.0e-16

int fgmr(SMatptr Amat, SPreptr lu, double *rhs, double *sol, 
         double tol, int im, int *itmax, FILE *fits) { 
/*----------------------------------------------------------------------
|                 *** Preconditioned FGMRES ***                  
+-----------------------------------------------------------------------
| This is a simple version of the ARMS preconditioned FGMRES algorithm. 
+-----------------------------------------------------------------------
| Y. S. Dec. 2000. -- Apr. 2008   -- Jul. 2009 
+-----------------------------------------------------------------------
| on entry:
|---------- 
|
|(Amat)   = matrix struct. the matvec operation is Amat->matvec.
|(lu)     = preconditioner struct.. the preconditioner is lu->precon
|           if (lu == NULL) the no-preconditioning option is invoked.
| rhs     = real vector of length n containing the right hand side.
| sol     = real vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| im      = Krylov subspace dimension 
| (itmax) = max number of iterations allowed. 
| fits    = NULL: no output
|        != NULL: file handle to output " resid vs time and its" 
|
| on return:
|---------- 
| fgmr      int =  0 --> successful return.
|           int =  1 --> convergence not achieved in itmax iterations.
| sol     = contains an approximate solution (upon successful return).
| itmax   = has changed. It now contains the number of steps required
|           to converge -- 
+-----------------------------------------------------------------------
| internal work arrays:
|----------       
| vv      = work array of length [im+1][n] (used to store the Arnoldi
|           basis)
| hh      = work array of length [im][im+1] (Arnoldi matrix)
| z       = work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
|     matvec and
|     preconditionning operation 
+---------------------------------------------------------------------*/
  int n=Amat->n, maxits = *itmax; 
  int i, i1, ii, j, k, k1, its, im1, pti, pti1, ptih=0, retval, one = 1;
  double *hh, *c, *s, *rs, t;
  double negt, beta, eps1=0, gam, *vv, *z; 
  im1 = im+1;
  vv = (double *)Malloc(im1*n*sizeof(double), "fgmres:vv");
  z  = (double *)Malloc(im*n*sizeof(double), "fgmres:z");
  im1 = im+1;
  hh = (double *)Malloc((im1*(im+3))*sizeof(double), "fgmres:hh");
  c  = hh+im1*im ; s  = c+im1;  rs = s+im1;
/*-------------------- outer loop starts here */
  retval = 0;
  its = 0;
/*-------------------- Outer loop */
  while (its < maxits) {
/*-------------------- compute initial residual vector */
    Amat->matvec(Amat, sol, vv); 
    for (j=0; j<n; j++)
      vv[j] = rhs[j] - vv[j];    /*  vv[0]= initial residual */
    beta = DNRM2(n, vv, one);
/*-------------------- print info if fits != null */
    if (fits != NULL && its == 0) 
      fprintf(fits, "%8d   %10.2e\n",its, beta) ;     
    if (beta == 0.0) 
      break;
    t = 1.0 / beta;
/*--------------------   normalize:  vv    =  vv   / beta */
    DSCAL(n, t, vv, one);
    if (its == 0) 
      eps1 = tol*beta;
/*--------------------initialize 1-st term  of rhs of hessenberg mtx */
    rs[0] = beta;
    i = 0;
/*-------------------- Krylov loop*/
    i = -1;
    pti=pti1=0;
    while((i < im-1) && (beta > eps1) && (its++ < maxits))  {
      i++;
      i1   = i+1; 
      pti  = i*n;
      pti1 = i1*n;
/*------------------------------------------------------------
|  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
+-----------------------------------------------------------*/
      if (lu == NULL) 
        memcpy(z+pti, vv+pti, n*sizeof(double));
      else
        lu->precon(vv+pti, z+pti, lu) ;
/*-------------------- matvec operation w = A z_{j} = A M^{-1} v_{j} */
      Amat->matvec(Amat, &z[pti], &vv[pti1]); 
/*-------------------- modified gram - schmidt...
|     h_{i,j} = (w,v_{i});  
|     w  = w - h_{i,j} v_{i}
+------------------------------------------------------------*/
      ptih=i*im1;
      for (j=0; j<=i; j++) {
	t = DDOT(n, &vv[j*n], one, &vv[pti1], one);
	hh[ptih+j] = t;
	negt = -t;
	DAXPY(n, negt, &vv[j*n], one, &vv[pti1], one);
      }
/*-------------------- h_{j+1,j} = ||w||_{2}    */
      t = DNRM2(n, &vv[pti1], one);
      hh[ptih+i1] = t;
      if (t == 0.0) 
	return(1);
      t = 1.0/t;
/*-------------------- v_{j+1} = w / h_{j+1,j}  */
      DSCAL(n, t, &vv[pti1], one);
/*-------- done with modified gram schimdt/arnoldi step
| now  update factorization of hh.
| perform previous transformations  on i-th column of h
+-------------------------------------------------------*/
      for (k=1; k<=i; k++) {
	k1 = k-1;
	t = hh[ptih+k1];
	hh[ptih+k1] = c[k1]*t + s[k1]*hh[ptih+k];
	hh[ptih+k] = -s[k1]*t + c[k1]*hh[ptih+k];
      }
      gam = sqrt(pow(hh[ptih+i],2) + pow(hh[ptih+i1],2) );
/*-------------------- check if gamma is zero */
      if (gam == 0.0) gam = epsmac;
/*-------------------- get  next plane rotation    */
      c[i] = hh[ptih+i]/gam;
      s[i] = hh[ptih+i1]/gam;
      rs[i1] = -s[i]*rs[i];
      rs[i] =  c[i]*rs[i];
/*-------------------- get residual norm + test convergence*/
      hh[ptih+i] = c[i]*hh[ptih+i] + s[i]*hh[ptih+i1];
      beta = fabs(rs[i1]);
      if( fits != NULL ) 
	fprintf(fits,"%8d   %10.2e\n", its, beta) ;          
/*-------------------- end [inner] while loop [Arnoldi] */
    }
/*-------------------- now compute solution. 1st, solve upper 
                       triangular system*/
    rs[i] = rs[i]/hh[ptih+i];
    for (ii=i-1; ii>=0; ii--) {
      t=rs[ii];
      for (j=ii+1; j<=i; j++)
        t -= hh[j*im1+ii]*rs[j];
      rs[ii] = t/hh[ii*im1+ii];
    }
/*---------- linear combination of z_j's to get sol. */
    for (j=0; j<= i; j++) 
      DAXPY(n, rs[j], &z[j*n], one, sol, one);
/*--------------------  restart outer loop if needed */
    if (beta < eps1) 
      break;
    else 
      if (its >= maxits) 
	retval = 1; 
/*---------- end main [outer] while loop */
  } 
/*-------------------- prepare to return */
  *itmax = its; 
  free(vv);
  free(z);
  free(hh);
  return (retval); 
}
/*-----------------end of fgmr -----------------------------------*/
