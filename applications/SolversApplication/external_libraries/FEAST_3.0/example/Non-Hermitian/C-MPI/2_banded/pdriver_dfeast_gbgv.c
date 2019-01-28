/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! FEAST Driver Example - Banded Storage 
  !!!!!!! solving Ax=ex with A real non-symmetric (non-Hermitian)
  !!!!!!! James Kestyn, Eric Polizzi 2015
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#include <stdio.h> 
#include <stdlib.h> 
#include <sys/time.h>
#include <mpi.h> 

#include "feast.h"
#include "feast_banded.h"

int main(int argc, char **argv) {
  /*!!!!!!!!!!!!!!!!! Feast declaration variable */
  int  feastparam[64]; 
  int loop,LDA,LDB;

  /*!!!!!!!!!!!!!!!!! Matrix declaration variable */
  FILE *fp;
  char name[]="../../system3";
  int  N,nnz,kl,ku;
  double *A,*B;

  /*!!!!!!!!!!!!!!!!! Others */
  struct timeval t1, t2;
  int  i,j,k,n2,err;
  int  M0,M,info;
  double r,trace[2],epsout,dum;
  double Emid[2];
  double *X; //! eigenvectors
  double *E,*res; //! eigenvalue+residual
/*********** MPI *****************************/
int rank,numprocs;
MPI_Init(&argc,&argv); 
MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
/*********************************************/


  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!! read input file in coordinate format!!!!!!!
    !!!!!!!!!!!!!!!!transform  in banded format directly !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  // !! find kl,ku (kl,ku here)
  kl=0;
  ku=0;
  fp = fopen (name, "r");
  err=fscanf (fp, "%d%d\n",&N,&nnz);
  for (k=0;k<=nnz-1;k++){

    err=fscanf(fp,"%d%d%lf%lf\n",&i,&j,&dum,&dum);
    if( i > j ){
        if (abs(i-j)>kl) kl=abs(i-j); 
    }else{
        if (abs(i-j)>ku) ku=abs(i-j); 
    }
 
  }
  fclose(fp);
  // !! form the banded matrices A and B
  LDA=(kl+ku+1);
  LDB=LDA;
  n2=2*N*(kl+ku+1); 
  A=calloc(n2,sizeof(double));
  B=calloc(n2,sizeof(double));
  for (i=0;i<=n2-1;i++){
    *(A+i)=(double) 0.0;
    *(B+i)=(double) 0.0;
  };

  fp = fopen (name, "r");
  err=fscanf (fp, "%d%d\n",&N,&nnz); 
  for (k=0;k<=nnz-1;k++){
    err=fscanf(fp,"%d%d",&i,&j);
    err=fscanf(fp,"%lf%lf\n",A+(j-1)*LDA+ku+(i-j),B+(j-1)*LDB+ku+(i-j));
  };
  fclose(fp);

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  if(rank==0) printf("dense matrix -system3- size %.d\n",N);
  if(rank==0) printf("bandwidth %d \n",kl+ku+1);

  gettimeofday(&t1,NULL);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FEAST in banded format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*!!! search interval [Emin,Emax] including M eigenpairs*/
  M0=40; // !! M0>=M
  Emid[0] = 0.59e0;
  Emid[1] = 0.41e0;
  r = 0.51e0;

  /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
  E=calloc(M0*2,sizeof(double));  // eigenvalues  // factor 2 for complex
  res=calloc(M0*2,sizeof(double));// eigenvectors // factor 2 for L and R res
  X=calloc(2*N*M0*2,sizeof(double));// residual (if needed) // factor 2 for complex //factor 2 for L/R vectors


  /*!!!!!!!!!!!!  FEAST */
  feastinit(feastparam);
  feastparam[0]=1;  /*change from default value */
  dfeast_gbgv(&N,&kl,&ku,A,&LDA,&kl,&ku,B,&LDB,feastparam,&epsout,&loop,Emid,&r,&M0,E,X,&M,res,&info);
  gettimeofday(&t2,NULL);
  /*!!!!!!!!!! REPORT !!!!!!!!!*/
  if(rank==0) printf("FEAST OUTPUT INFO %d\n",info);
  if (info==0 && rank==0) {
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n"); 
    printf("SIMULATION TIME %f\n",(t2.tv_sec-t1.tv_sec)*1.0+(t2.tv_usec-t1.tv_usec)*0.000001);
    printf("# Search interval [Emid; r] %.15e+i*%.15e;%.15e\n",Emid[0],Emid[1],r);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace[0]=(double) 0.0;
    trace[1]=(double) 0.0;
    for (i=0;i<=M-1;i=i+1){
      trace[0]=trace[0]+*(E+2*i);
      trace[1]=trace[1]+*(E+2*i+1);
    }	  
    printf("TRACE %.15e + i*%.15e \n", trace[0],trace[1]);
    printf("Relative error on the Trace %.15e\n",epsout );
    printf("Eigenvalues/Residuals\n");
    for (i=0;i<=M-1;i=i+1){
      printf("   %d %.15e+i*%.15e \t %.15e\n",i,*(E+2*i),*(E+2*i+1),*(res+i));
    }
  }
  MPI_Finalize();
  return 0;
}






