#include <stdio.h> 
#include <stdlib.h> 

#include "feast.h"
#include "feast_dense.h"
int main() { 
  /*  eigenvalue system  */
  int   N=2,LDA=2;
  char  UPLO='F'; // 'L' and 'U' also fine
  double A[4]={2.0,-1.0,-1.0,2.0}; 
  double Emin=-5.0, Emax=5.0;
  int M0=2; //size initial subspace
  /* input parameters for FEAST */
  int feastparam[64];
  /* output variables for FEAST */
  double *E, *res, *X; 
  double  epsout,trace; 
  int i,loop,info,M;

  /* Allocate memory for eigenvalues.eigenvectors/residual */
  E=calloc(M0,sizeof(double));  //eigenvalues
  res=calloc(M0,sizeof(double));//eigenvectors
  X=calloc(N*M0,sizeof(double));//residual 

  /* !!!!!!!!!! FEAST !!!!!!!!!*/
  feastinit(feastparam); 
  feastparam[0]=1;  /*change from default value */
  dfeast_syev(&UPLO,&N,A,&LDA,feastparam,&epsout,&loop,&Emin,&Emax,&M0,E,X,&M,res,&info);

  /*!!!!!!!!!! REPORT !!!!!!!!!*/
  printf("FEAST OUTPUT INFO %d\n",info);
  if (info==0) {
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("# Search interval [Emin,Emax] %.15e %.15e\n",Emin,Emax);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace=0.0;
    for (i=0;i<=M-1;i=i+1){
      trace=trace+*(E+i);
    }	  
    printf("TRACE %.15e\n", trace);
    printf("Relative error on the Trace %.15e\n",epsout );
    printf("Eigenvalues/Residuals\n");
    for (i=0;i<=M-1;i=i+1){
      printf("   %d %.15e %.15e\n",i,*(E+i),*(res+i));
    }
    printf("Eigenvectors\n");
    for (i=0;i<=M-1;i=i+1){
      printf("   %d (%.15e, %.15e)\n",i,*(X+i*M),*(X+1+i*M));
    }
  }
free(E);
free(res);
free(X);
  return 0;
}




