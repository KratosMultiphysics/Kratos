/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! FEAST Driver dense example !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! solving Ax=eBx with A real and B spd --- A and B dense matrix!!!!
  !!!!!!! by Eric Polizzi- 2009-2012!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#include <stdio.h> 
#include <stdlib.h> 
#include <sys/time.h>

#include "feast.h"
#include "feast_dense.h"
int main() {
  /*!!!!!!!!!!!!!!!!! Feast declaration variable */
  int  feastparam[64]; 
  double epsout;
  int loop;
  char UPLO='F'; // ! 'L' or 'U' also fine

  /*!!!!!!!!!!!!!!!!! Matrix declaration variable */
  FILE *fp;
  char name[]="../../system1";
  int  N,nnz;
  double *A,*B;

  /*!!!!!!!!!!!!!!!!! Others */
  struct timeval t1, t2;
  int  i,j,k,n2,err;
  int  M0,M,info;
  double Emin,Emax,trace;
  double *X; //! eigenvectors
  double *E,*res; //! eigenvalue+residual


  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!Read Coordinate format and convert to dense format
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  fp = fopen (name, "r");
  err=fscanf (fp, "%d%d\n",&N,&nnz);
  n2=N*N;
  A=calloc(n2,sizeof(double));
  B=calloc(n2,sizeof(double));
  for (i=0;i<=n2-1;i++){
    *(A+i)=(double) 0.0;
    *(B+i)=(double) 0.0;
  };
  for (k=0;k<=nnz-1;k++){
    err=fscanf(fp,"%d%d",&i,&j);
    err=fscanf(fp,"%lf%lf\n",A+(j-1)*N+i-1,B+(j-1)*N+i-1);
  };
  fclose(fp);
    
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  printf("dense matrix -system1- size %.d\n",N);


  gettimeofday(&t1,NULL);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FEAST in dense format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*!!! search interval [Emin,Emax] including M eigenpairs*/
  Emin=(double) 0.18;
  Emax=(double) 1.0;
  M0=25; // !! M0>=M

  /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
  E=calloc(M0,sizeof(double));  // eigenvalues
  res=calloc(M0,sizeof(double));// eigenvectors 
  X=calloc(N*M0,sizeof(double));// residual 


  /*!!!!!!!!!!!!  FEAST */
  feastinit(feastparam);
  feastparam[0]=1;  /*change from default value */
  dfeast_sygv(&UPLO,&N,A,&N,B,&N,feastparam,&epsout,&loop,&Emin,&Emax,&M0,E,X,&M,res,&info);

  gettimeofday(&t2,NULL);
  /*!!!!!!!!!! REPORT !!!!!!!!!*/
  printf("FEAST OUTPUT INFO %d\n",info);
  if (info==0) {
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("SIMULATION TIME %f\n",(t2.tv_sec-t1.tv_sec)*1.0+(t2.tv_usec-t1.tv_usec)*0.000001);
    printf("# Search interval [Emin,Emax] %.15e %.15e\n",Emin,Emax);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace=(double) 0.0;
    for (i=0;i<=M-1;i=i+1){
      trace=trace+*(E+i);
    }	  
    printf("TRACE %.15e\n", trace);
    printf("Relative error on the Trace %.15e\n",epsout );
    printf("Eigenvalues/Residuals\n");
    for (i=0;i<=M-1;i=i+1){
      printf("   %d %.15e %.15e\n",i,*(E+i),*(res+i));
    }
  }
  return 0;
}






