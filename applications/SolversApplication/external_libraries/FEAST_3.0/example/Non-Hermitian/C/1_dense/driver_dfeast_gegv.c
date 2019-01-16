/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! FEAST Driver Example  - Dense Storage
  !!!!!!! solving Ax=ex with A real non-symmetric 
  !!!!!!! James Kestyn, Eric Polizzi 2015
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#include <stdio.h> 
#include <stdlib.h> 
#include <sys/time.h>


#include "feast.h"
#include "feast_dense.h"
int main() {
  /*!!!!!!!!!!!!!!!!! Feast declaration variable */
  int  feastparam[64]; 
  int loop;
  char UPLO='F'; // ! 'L' or 'U' also fine

  /*!!!!!!!!!!!!!!!!! Matrix declaration variable */
  FILE *fp;
  char name[]="../../system3";
  int  N,LDA,LDB,nnz;
  double *A,*B;
  double *Zedge,*Zne,*Wne;
  int ccN, *Nedge, *Tedge;
  double *Emid,epsout,*trace;
  double r;
  /*!!!!!!!!!!!!!!!!! Others */
  struct timeval t1, t2;
  int  i,j,k,n2,err;
  int  M0,M,info;
  double *E,*X; //! eigenvectors
  double *res; //! eigenvalue+resridual

  Emid = calloc(2,sizeof(double));
  trace = calloc(2,sizeof(double));
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!Read Coordinate format and convert to dense format
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  fp = fopen (name, "r");
  err=fscanf (fp, "%d%d\n",&N,&nnz);
  printf ("n/nnx = %d %d\n",N,nnz);
  LDA=N;
  n2=N*N;  // factor 2 because of complex number
  A=calloc(n2,sizeof(double));
  B=calloc(n2,sizeof(double));
  for (i=0;i<=n2-1;i++){
    *(A+i)=(double) 0.0;
    *(B+i)=(double) 0.0;
  };
  for (k=0;k<=nnz-1;k++){
    err=fscanf(fp,"%d %d",&i,&j);
    err=fscanf(fp,"%lf%lf\n",A+(j-1)*N+i-1,B+(j-1)*N+i-1);
  };
  fclose(fp);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  printf("dense matrix -system3- size %.d\n",N);
  
  gettimeofday(&t1,NULL);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FEAST in dense format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  feastinit(feastparam);
  feastparam[0]=1;  /*change from default value */
  /*!!! search interval [Emin,Emax] including M eigenpairs*/
  M0=40; // !! M0>=M
  Emid[0] = 0.590e0;
  Emid[1] = 0.0e0;
  r = 0.410e0;

  /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
  E=calloc(2*M0,sizeof(double));  // eigenvalues
  X=calloc(2*N*M0*2,sizeof(double));// right and left eigenvector // factor 2 because of complex number
  res=calloc(M0*2,sizeof(double));// right and left eigenvectors residual


  printf("Enter FEAST\n");
  /*!!!!!!!!!!!!  FEAST */
  dfeast_gegv(&N,A,&N,B,&N,feastparam,&epsout,&loop,Emid,&r,&M0,E,X,&M,res,&info);
  gettimeofday(&t2,NULL);
  /*!!!!!!!!!! REPORT !!!!!!!!!*/
  info=0;
  printf("FEAST OUTPUT INFO %d\n",info);
  if (info==0) {
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("SIMULATION TIME %f\n",(t2.tv_sec-t1.tv_sec)*1.0+(t2.tv_usec-t1.tv_usec)*0.000001);
    printf("# Search interval [Emid,r] (%.15e+i*%.15e) %.15e\n",Emid[0],Emid[1],r);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace[0]=(double) 0.0;
    trace[1]=(double) 0.0;
    for (i=0;i<=M-1;i=i+1){
      trace[0]=trace[0]+*(E+2*i);
      trace[1]=trace[1]+*(E+2*i+1);
    }	  
    printf("TRACE %.15e + i*%.15e \n", trace[0],trace[1]);
    printf("Relative error on the Trace %.15e\n",epsout);
    printf("Eigenvalues/Residuals\n");
    for (i=0;i<=M-1;i=i+1){
      printf("   %d %.15e+i*%.15e \t %.15e\n",i,*(E+2*i),*(E+2*i+1),*(res+i));
    }
  }
  return 0;
}





