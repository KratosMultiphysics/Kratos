/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! FEAST Driver sparse example !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! solving Ax=eBx with A real and B spd --- A and B sparse matrix!!!
  !!!!!!! by Eric Polizzi- 2009-2012!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#include <stdio.h> 
#include <stdlib.h> 
#include <sys/time.h>

#include "feast.h"
#include "feast_sparse.h"
int main() {
  /*!!!!!!!!!!!!!!!!! Feast declaration variable */
  int  feastparam[64]; 
  double epsout;
  int loop;
  char UPLO='F'; 

  /*!!!!!!!!!!!!!!!!! Matrix declaration variable */
  FILE *fp;
  char name[]="../../system1";
  int  N,nnz;
  double *sa,*sb;
  int *isa,*jsa;
  /*!!!!!!!!!!!!!!!!! Others */
  struct timeval t1, t2;
  int  i,k,err;
  int  M0,M,info;
  double Emin,Emax,trace;
  double *X; //! eigenvectors
  double *E,*res; //! eigenvalue+residual


  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!! read input file in csr format!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  // !!!!!!!!!! form CSR arrays isa,jsa,sa 
  fp = fopen (name, "r");
  err=fscanf (fp, "%d%d\n",&N,&nnz);
  sa=calloc(nnz,sizeof(double));
  sb=calloc(nnz,sizeof(double));
  isa=calloc(N+1,sizeof(int));
  jsa=calloc(nnz,sizeof(int));

  for (i=0;i<=N;i++){
    *(isa+i)=0;
  };
  *(isa)=1;
  for (k=0;k<=nnz-1;k++){
    err=fscanf(fp,"%d%d%lf%lf\n",&i,jsa+k,sa+k,sb+k);
    *(isa+i)=*(isa+i)+1;
  };
  fclose(fp);
  for (i=1;i<=N;i++){
    *(isa+i)=*(isa+i)+*(isa+i-1);
  };

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  printf("sparse matrix -system1- size %.d\n",N);
  printf("nnz %d \n",nnz);

  gettimeofday(&t1,NULL);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FEAST in sparse format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  dfeast_scsrgv(&UPLO,&N,sa,isa,jsa,sb,isa,jsa,feastparam,&epsout,&loop,&Emin,&Emax,&M0,E,X,&M,res,&info);

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






